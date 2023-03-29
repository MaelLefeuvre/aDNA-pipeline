configfile: "./config/config.yml"

# ------------------------------------------------------------------------------------------------------------------- #
# ---- Generic / Utility rules
rule samtools_index:
    """
    Index a generic .bam file
    """
    input:
        bam = "{bam}"
    output:
        bai = "{bam}.bai"
    threads: 8
    conda: "../envs/samtools-1.15.yml"
    log: "logs/generics/{bam}/samtools_index.log"
    shell: """
        samtools index -@ {threads} {input.bam}
    """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 00. Resolve the required aligner algorithm.

def assign_aligner_algorithm(wildcards):
    """
    Decide on the appropriate bwa algorithm (aln or mem), based on the user input.
    """
    aligner        = config['preprocess']['bwa']['aligner']
    collapsed_only = config['preprocess']['bwa']['collapsed-only']
    protocol       = get_illumina_protocol(wildcards)

    if protocol == "single":
        if aligner == "mem":
            return expand(rules.bwa_mem_se.output.sam, sample="{sample}", run="{run}", extension="truncated")
        elif aligner == "aln":
            return expand(rules.bwa_samse.output.sam, sample="{sample}", run="{run}", extension="truncated")
            
    elif protocol == "paired":
        if aligner == "mem":
            if collapsed_only:
                return expand(rules.bwa_mem_se.output.sam, sample="{sample}", run="{run}", extension="collapsed")
            else:
                return rules.samtools_merge_mem.output.merged
        elif aligner == "aln":
            if collapsed_only:
                return expand(rules.bwa_samse.output.sam, sample="{sample}", run="{run}", extension="collapsed")
            else:
                return rules.samtools_merge_aln.output.merged


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 01. Perform read-length and BQ quality filtration.

rule samtools_filter_unmapped:
    """
    Apply QC filtration after raw alignment. Three filters are applied concurrently:
    - (user-defined) minimum base quality
    - (user-defined) minimum length
    - (constant)     remove unmapped
    """
    input:
        sam        = assign_aligner_algorithm,
        reference  = config["reference"],
        bwt        = rules.index_reference_genome.output.bwt
    output:
        bam        = temp("results/01-preprocess/03-filter/{sample}/{run}/{sample}.bam")
    params:
        min_MQ     = config["preprocess"]["filter"]["min-MQ"],
        min_length = config["preprocess"]["filter"]["min-length"],
    log:       "logs/01-preprocess/samtools_filter_unmapped/{sample}/{run}.log"
    benchmark: "benchmarks/01-preprocess/samtools_filter_unmapped/{sample}/{run}.tsv"
    conda:     "../envs/samtools-1.15.yml"
    priority:  6
    threads:   8
    shell: """
        samtools view \
        --threads {threads} \
        --reference {input.reference} \
        -q {params.min_MQ} \
        -e 'length(seq)>{params.min_length}' \
        -F4 -Sb \
        {input.sam} > {output.bam} 2> {log}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 02. Merge BAM files across runs. 

def define_merged_bams(wildcards):
    """
    Define the appropriate input bam for the variant caller, based on which 
    PMD-rescaling method was requested by the user.
    """
    # Run through the initial samples files and extract pedigree ids 
    # run_ids = os.listdir(f"original-data/samples/{wildcards.sample}")
    run_ids = get_requested_sample_runs(wildcards)
    
    #print(get_requested_sample_runs(wildcards))

    # Return a list of input bam files for pileup
    return expand(
        rules.samtools_filter_unmapped.output.bam,
        sample = "{sample}",
        run    = run_ids
    )
    

rule samtools_merge_runs:
    input:
        bams = define_merged_bams
    output:
        merged = "results/01-preprocess/04-merge-runs/{sample}/{sample}.merged.bam"
    params:
        RG     = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina' # @TODO ADD it LATER
    log:       "logs/01-preprocess/samtools_merge_runs/{sample}.log"
    benchmark: "benchmarks/01-preprocess/samtools_merge_runs/{sample}.tsv"
    conda:     "../envs/samtools-1.15.yml"
    threads:   8
    shell: """
        samtools merge -@ {threads} -o - {input.bams} \
        | samtools addreplacerg -r '{params.RG}' -w -OBAM -o {output.merged} - 2> {log}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 03. Sort BAM file. 

rule samtools_sort:
    """
    Apply coordinate sorting on a BAM file.
    """
    input:
        bam       = rules.samtools_merge_runs.output.merged,
        reference = config["reference"]
    output:
        bam       = "results/01-preprocess/05-sort/{sample}/{sample}.srt.bam"
    log:       "logs/01-preprocess/samtools_sort/{sample}.log"
    benchmark: "benchmarks/01-preprocess/samtools_sort/{sample}.tsv"
    conda:     "../envs/samtools-1.15.yml"
    priority: 7
    threads:  8
    shell: """
        samtools sort -@ {threads} --reference {input.reference} --output-fmt BAM -o {output.bam} {input.bam} 2> {log}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 04. Remove duplicates using Picard, Dedup or Samtools.

rule picard_rmdup:
    """
    Remove PCR duplicates from a BAM file using picard.
    """
    input:
        bam     = rules.samtools_sort.output.bam,
    output:
        bam     = "results/01-preprocess/06-dedup/picard/{sample}/{sample}.srt.rmdup.bam",
        metrics = "results/01-preprocess/06-dedup/picard/{sample}/{sample}.rmdup.metrics.txt"
    params:
        tmpdir  = config["tempdir"]
    log:       "logs/01-preprocess/picard_rmdup/{sample}.log"
    benchmark: "benchmarks/01-preprocess/picard_rmdup/{sample}.tsv"
    conda:     "../envs/picard-2.27.4.yml"
    threads:   8
    shell: """
        picard MarkDuplicates \
        -I {input.bam} \
        -O {output.bam} \
        -M {output.metrics} \
        --ASSUME_SORT_ORDER coordinate \
        --REMOVE_DUPLICATES true \
        --VALIDATION_STRINGENCY LENIENT \
        --TMP_DIR {params.tmpdir} 2> {log}
    """

rule apeltzer_dedup:
    """
    Remove PCR duplicates from a BAM file using Alexander Peltzer's 'dedup'
    tool (i.e. from nf-core EAGER pipeline)
    @ TODO: FIXME: not yet fully implemented: AdapterRemovalFixPrefix
            should be applied on fastq files, right after AdapterRemoval not bams.
    """
    input:
        bam = rules.samtools_sort.output.bam
    output:
        bam = temp("results/01-preprocess/06-dedup/dedup/{sample}/{sample}.srt.rmdup.bam")
    log: "logs/01-preprocess/apeltzer_dedup/{sample}.log"
    conda: "../envs/dedup-0.12.8.yml"
    shell: """
        AdapterRemovalFixPrefix {input.bam} | dedup --output {output.bam} 2> {log}
    """


rule samtools_rmdup:
    """
    Remove PCR duplicates from a BAM file using samtool's rmdup
    """
    input:
        bam = rules.samtools_sort.output.bam
    output:
        bam = temp("results/01-preprocess/06-dedup/samtools/{sample}/{sample}.srt.rmdup.bam")
    params:
        tmpdir = lambda wildcards, output: splitext(output.bam)[0]
    log:       "logs/01-preprocess/samtools_rmdup/{sample}.log"
    benchmark: "benchmarks/01-preprocess/samtools_rmdup/{sample}.tsv"
    conda:     "../envs/samtools-1.15.yml"
    threads: 8
    shell: """
    samtools index -@ {threads} {input.bam}
    samtools collate -O {input.bam} \
    | samtools fixmate -m - - \
    | samtools sort - \
    | samtools markdup -r -T {params.tmpdir} -s - {output.bam} 2> {log}
    """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 05. Perform PMD Base Quality rescaling using MapDamage 

def define_dedup_input_bam(wildcards):
    """
    Choose the appropriate input for run_pmdtools / run_mapdamage, based on
    which duplicate removal tool was requested by the user.
    """
    match config['preprocess']['dedup']['method']:
        case "picard":
            return rules.picard_rmdup.output.bam
        case "dedup":
            return rules.apeltzer_dedup.output.bam
        case "samtools":
            return rules.samtools_rmdup.output.bam

    raise RuntimeError(f"Invalid config value for dedup method '{config['preprocess']['dedup']['method']}'.")


rule run_pmdtools:
    input:
        bam       = define_dedup_input_bam,
        bai       = lambda wildcards: define_dedup_input_bam(wildcards) + ".bai",
        reference = config['reference'],
    output:
        bam       = "results/01-preprocess/07-rescale/pmdtools/{sample}/{sample}.srt.rmdup.filtercontam.bam"
    params:
        threshold  = config['preprocess']['pmd-rescaling']['pmdtools']['threshold'],
        mask_deams = config['preprocess']['pmd-rescaling']['pmdtools']['mask-terminal-deams'],
    log:       "logs/01-preprocess/run_pmdtools/{sample}.log"
    benchmark: "benchmarks/01-preprocess/run_pmdtools/{sample}.tsv"
    conda:     "../envs/pmdtools-0.60.yml"
    threads:   3
    shell: """
        ### --threshold {params.threshold} 
        samtools calmd {input.bam} {input.reference} \
        | pmdtools \
          --header \
          --adjustbaseq \
          --stats \
          --threshold -9999 \
          --upperthreshold 999999 \
          --maskterminaldeaminations {params.mask_deams} \
        | samtools view -OBAM > {output.bam} 2> {log}
    """


def define_mapdamage_seed(wildcards):
    """
    Add additional seeding arguments for MapDamage if the user-requested for 
    downsampling. If the downsample seed was provided within the config file,
    provide the program with said value. If not, fetch the randomly generated
    backup seed from our metadata file.
    """
    downsample_n = config['preprocess']['pmd-rescaling']['map-damage']['downsample']

    if downsample_n is None:
        return ""

    seed = config['preprocess']['pmd-rescaling']['map-damage']['downsample-seed']
    if seed is None:
        with open(rules.meta.output.metadata) as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']

    return f"--downsample {downsample_n} --downsample-seed {seed}"


rule run_mapdamage:
    """
    Apply PMD base quality score recalibration on a bam file using MapDamage.
    """
    input:
        bam       = define_dedup_input_bam,
        bai       = lambda wildcards: define_dedup_input_bam(wildcards) + ".bai",
        reference = config["reference"],
        metadata     = "results/meta/pipeline-metadata.yml"
    output:
        bam = "results/01-preprocess/07-rescale/mapdamage/{sample}/{sample}.srt.rmdup.rescaled.bam",
        misincorporation = "results/01-preprocess/07-rescale/mapdamage/{sample}/misincorporation.txt"
    params:
        downsample_seed = define_mapdamage_seed
    log:       "logs/01-preprocess/run_mapdamage/{sample}.log"
    benchmark: "benchmarks/01-preprocess/run_mapdamage/{sample}.tsv"
    conda:     "../envs/mapdamage-2.2.1.yml"
    threads: 1
    priority: 10
    shell: """
        mapDamage \
        -i {input.bam}                   \
        -r {input.reference}             \
        --rescale                        \
        --folder $(dirname {output.bam}) \
        --rescale-out {output.bam}       \
        {params.downsample_seed}         \
        --verbose > {log} 2>&1 
    """


def define_masking_input_bam(wildcards):
    rescaler = config['preprocess']['pmd-rescaling']['rescaler']
    if rescaler is None:
        return define_dedup_input_bam(wildcards)
    else:
        return define_rescale_input_bam(wildcards)
    raise RuntimeError("Failed to define a proper input bam file for pmd-mask") 

rule run_pmd_mask:
    input:
        bam              = define_masking_input_bam,
        bai              = lambda wildcards: define_masking_input_bam(wildcards) + ".bai",
        reference        = config["reference"],
        misincorporation = rules.run_mapdamage.output.misincorporation
    output:
        bam       = "results/01-preprocess/07-rescale/pmd-mask/{sample}/{sample}.pmd_masked.bam",
        metrics   = "results/01-preprocess/07-rescale/pmd-mask/{sample}/{sample}.pmd_masked.metrics"
    params:
        threshold = config['preprocess']['pmd-rescaling']['pmd-mask']['threshold']
    log:   "logs/01-preprocess/run_pmd_mask/{sample}.log"
    benchmark: "benchmarks/01-preprocess/run_pmd_mask/{sample}.tsv"
    conda: "../envs/pmd-mask.yml"
    threads: 8
    shell: """
        pmd-mask -@ {threads} -b {input.bam} -f {input.reference} -m {input.misincorporation} --threshold {params.threshold} -M {output.metrics} -Ob -o {output.bam} --verbose 2> {log}
    """
