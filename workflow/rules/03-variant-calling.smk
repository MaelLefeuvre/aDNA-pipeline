localrules: generate_bam_list
# ------------------------------------------------------------------------------------------------------------------- #
# ---- 01. Generate a list of input BAM list for cohort variant calling.

def define_rescale_input_bam(wildcards):
    rescaler = config['preprocess']['pmd-rescaling']['rescaler']
    match rescaler:
        case  "mapdamage":
            return "results/01-preprocess/07-rescale/mapdamage/{sample}/{sample}.srt.rmdup.rescaled.bam" 
        case "pmdtools":
            return "results/01-preprocess/07-rescale/pmdtools/{sample}/{sample}.srt.rmdup.filtercontam.bam"
        case other:
            raise RuntimeError(f'Invalid rescaler value "{rescaler}')


def get_pileup_input_bams(wildcards):
    """
    Define the appropriate input bam for the variant caller, based on which 
    PMD-rescaling method was requested by the user.
    """
    # Run through the config file and extract sample ids 
    samples = get_sample_names(wildcards)

    # if masking is required, delegate input definition to the appropraite rule.
    apply_masking = config['preprocess']['pmd-rescaling']['apply-masking']
    if apply_masking:
        print("Applying pmd-mask for variant calling.", file=sys.stderr)
        return expand(rules.run_pmd_mask.output.bam, sample = samples)


    # Return a list of input bam files for pileup
    rescaler = config['preprocess']['pmd-rescaling']['rescaler']
    if rescaler is None:
        print("WARNING: Skipping PMD Rescaling for variant calling!", file=sys.stderr)
        return expand(define_dedup_input_bam(wildcards), sample = samples)
    else:
        print("Note: Applying {rescaler} for variant calling.", file=sys.stderr)
        return expand(define_rescale_input_bam(wildcards), sample = samples)

    raise RuntimeError(f'Invalid rescaler value "{rescaler}')


rule generate_bam_list:
    """
    Create an ordered list of bam files using the original samples file
    """
    input:
        bamlist = get_pileup_input_bams
    output:
        bamlist = "results/02-variant-calling/01-pileup/samples.bam.list"
    log:   "logs/02-variant-calling/generate_bam_list.log"
    priority: 15
    shell: """
        ls {input.bamlist} > {output.bamlist} 2> {log}
    """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 02. Generate pileup files using samtools.

def get_samtools_optargs(wildcards):
    """
    Inject the appropriate flag to disable samtools' Base alignment quality,
    if the user requested it.
    """
    optargs = ""
    if config['variant-calling']['pileup']['disable-BAQ']:
        optargs += "-B "
    return optargs


rule samtools_pileup:
    """
    Generate a legacy pileup file using samtools mpileup.
    """
    input:
        bamlist   = rules.generate_bam_list.output.bamlist,
        targets   = os.path.splitext(config['variant-calling']['targets'])[0] + '.ucscbed',
        reference = config["reference"],
        fai       = config["reference"] + ".fai",
    output:
        pileup    = "results/02-variant-calling/01-pileup/samples.pileup"
    params:
        min_MQ    = config['variant-calling']['pileup']['min-MQ'],
        min_BQ    = config['variant-calling']['pileup']['min-BQ'],
        optargs   = get_samtools_optargs
    log:       "logs/02-variant-calling/samtools_pileup.log"
    benchmark: "benchmarks/02-variant-calling/samtools_pileup.tsv"
    conda:     "../envs/samtools-1.15.yml"
    shell: """
        samtools mpileup \
        -R {params.optargs} \
        -q {params.min_MQ} \
        -Q {params.min_BQ} \
        -l {input.targets} \
        -f {input.reference} \
        -b {input.bamlist} | LC_ALL=C sort -k1,1n -k2,2n > {output.pileup} 2> {log}
        """


rule get_bamlist_panel_coverage:
    """
    Obtain a raw estimate of the SNP target panel coverage for each sample.
    """
    input:
        bamlist  = rules.generate_bam_list.output.bamlist,
        targets  = rules.samtools_pileup.input.targets
    output:
        coverage = "results/02-variant-calling/01-pileup/bamlist-panel-coverage.tsv"
    log:       "logs/02-variant-calling/01-pileup/get_bamlist_panel_coverage.log"
    benchmark: "benchmarks/02-variant-calling/01-pileup/get_bamlist_panel_coverage.tsv"
    conda:     "../envs/samtools-1.15.yml"
    threads:   4
    shell: """
        panel_size=$(cat {input.targets} | wc -l)
        for bam in $(cat {input.bamlist}); do
            depth=$(samtools depth -@ {threads} -b {input.targets} $bam | awk '($3>0)' | wc -l);
            coverage=$(python -c "print(${{depth}}/${{panel_size}})")
            echo -e ${{bam}}"\t"${{depth}}\t${{coverage}};
        done | sort -h -k2 | column -t > {output.coverage} 2> {log}
    """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 03-A. Perform pseudo-haploid random variant calling with SequenceTools' PileupCaller

def parse_pileup_caller_flags(wildcards):
    args=""
    if config['variant-calling']['pileupCaller']["skip-transitions"]:
        args+="--skipTransitions "

    seed = config['variant-calling']['pileupCaller']['seed']
    if seed is None:
        with open(rules.meta.output.metadata) as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']

    args += f"--seed {seed} "

    match config['variant-calling']['pileupCaller']["mode"]:
        case "randomHaploid":
            args += "--randomHaploid "
        case "majorityCall":
            args += "--majorityCall "
        case _:
            raise RuntimeError("Incorrect pileupCaller mode selected.")
    return args


rule pileup_caller:
    """
    Perform random sampling variant calling uing SequenceTools' pileupCaller
    """
    input:
        pileup            = rules.samtools_pileup.output.pileup,
        bamlist           = rules.generate_bam_list.output.bamlist,
        targets           = os.path.splitext(config["variant-calling"]["targets"])[0] + ".snp",
        metadata          = "results/meta/pipeline-metadata.yml"
    output:
        plink             = multiext("results/02-variant-calling/02-pileupCaller/samples", ".bed", ".bim", ".fam"),
        sample_names_file = "results/02-variant-calling/02-pileupCaller/samples-names.txt"
    params:
        basename          = "results/02-variant-calling/02-pileupCaller/samples",
        optargs           = parse_pileup_caller_flags,
        min_depth         = config['variant-calling']["pileupCaller"]["min-depth"],
        seed              = config['variant-calling']['pileupCaller']["seed"],
        sample_pop_name   = config['variant-calling']['sample-pop-name']
    log:       "logs/02-variant-calling/pileup_caller.log"
    benchmark: "benchmarks/02-variant-calling/pileup_caller.tsv"
    conda:     "../envs/sequencetools-1.5.2.yml"
    shell: """
        cat {input.bamlist} | rev | cut -d'/' -f1 | rev | cut -d'.' -f1  > {output.sample_names_file}

        pileupCaller \
        {params.optargs} \
        --snpFile <( awk '($2<=22)' {input.targets}) \
        --sampleNameFile {output.sample_names_file} \
        --plinkOut {params.basename} \
        --minDepth {params.min_depth} \
        --samplePopName {params.sample_pop_name} < {input.pileup} > {log} 2>&1 
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 03-B. Perform pseudo-haploid random variant calling with ANGSD

rule ANGSD_random_haploid:
    """
    Perform random pseudo-haploid variant calling using ANGSD
    """
    input:
        pileup    = rules.samtools_pileup.output.pileup,
        bamlist   = rules.generate_bam_list.output.bamlist,
        reference = config['reference'],
        fai       = config["reference"] + ".fai",
    output:
        haplos    = "results/02-variant-calling/02-ANGSD/samples.haplo.gz"
    params:
        out       = "results/02-variant-calling/02-ANGSD/samples"
    log:       "logs/02-variant-calling/ANGSD_random_haploid.log"
    benchmark: "benchmarks/02-variant-calling/ANGSD_random_haploid.tsv"
    conda:     "../envs/angsd-0.939.yml"
    threads:   4
    priority:  15
    shell: """
        angsd -pileup {input.pileup} \
        -fai {input.fai} \
        -nInd $(cat {input.bamlist} | wc -l) \
        -rf <(seq 1 22) \
        -dohaplocall 1 \
        -doCounts 1 \
        -nthreads {threads} \
        -out {params.out} 2> {log}
    """


rule ANGSD_haplo_to_plink:
    """
    Convert a .haplo.gz ANGSD file to a set of PLINK .tped / .tfam files 
    """
    input:
        haplos  = rules.ANGSD_random_haploid.output.haplos,
        bamlist = rules.generate_bam_list.output.bamlist 
    output:
        tped = "results/02-variant-calling/02-ANGSD/samples.tped",
        tfam = "results/02-variant-calling/02-ANGSD/samples.tfam",
    params:
        outputname      = "results/02-variant-calling/02-ANGSD/samples",
        sample_pop_name = config['variant-calling']['sample-pop-name']

    log:       "logs/02-variant-calling/ANGSD_haplo_to_plink.log"
    benchmark: "benchmarks/02-variant-calling/ANGSD_haplo_to_plink.tsv"
    conda:     "../envs/angsd-0.939.yml"
    priority:  15
    shell: """
        haploToPlink {input.haplos} {params.outputname} 2>  {log}
        sed -i 's/N/0/g' {output.tped}                  2>> {log}
        cat {input.bamlist} \
        | xargs basename -a \
        | grep -oP '^[^.]+(?=(\.[^.]+)*(\.bam$))' \
        | awk 'BEGIN{{OFS="\t"}}{{print "{params.sample_pop_name}", $1, 0, 0, 0, 0}}' \
        > {output.tfam} 2>> {log}
    """
