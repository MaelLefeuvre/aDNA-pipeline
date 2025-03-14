import itertools
from os.path import splitext, basename
from functools import partial

configfile: "config/netrules.yml"

localrules: merge_TKGWV2_results


# ------------------------------------------------------------------------------------------------------------------- #



def get_TKGWV2_input_bams(wildcards):
    """
    Define the appropriate input bam files for TKGWV2, based on which PMD
    rescaler was requested by the user
    """
    # ---- If masking is required, delegate input definition to the appropraite rule.
    apply_masking = config['preprocess']['pmd-rescaling']['apply-masking']
    rescaler = config['preprocess']['pmd-rescaling']['rescaler']
    if apply_masking:
        msg  = "[NOTE]: Using pmd-masked files for TKGWV2"
        root = "results/01-preprocess/07-rescale/pmd-mask/{sample}/{sample}.pmd_masked.{ext}"
    elif rescaler is None:
        msg  = "[NOTE]: Skipping PMD rescaling for TKGWV2"
        root = define_dedup_input_bam(wildcards)
    else:
        msg  = "[Note]: Using {rescaler} rescaled bams for TKGWV2"
        root = define_rescale_input_bam(wildcards)

    # ---- Print message on first call only. 
    if getattr(get_TKGWV2_input_bams, 'first_call', None) is None:
        print(msg, file=sys.stderr)
        get_TKGWV2_input_bams.first_call = False

    return expand(root, sample = ["{pairA}", "{pairB}"], ext = ["bam", "bam.bai"] )


def TKGWV2_downsample_seed(wildcards):
    seed = config['kinship']['TKGWV2']['downsample-seed']
    if seed is None:
        with checkpoints.meta.get().output.metadata.open() as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']
    return seed


rule TKGWV2_downsample_bam:
    """
    Downsample the .bam file if there are more than 1_500_000 reads.
     - As is, this helper script is applied to all files ending with the ".bam" suffix within the working directory.
     - output files are identified with the "_subsampled.bam" suffix.
    """
    input:
        pairs    = get_TKGWV2_input_bams,
        metadata = rules.meta.output
    output:
        pairA = "results/03-kinship/TKGWV2/{pairA}_{pairB}/{pairA}.subsampled.bam",
        pairB = "results/03-kinship/TKGWV2/{pairA}_{pairB}/{pairB}.subsampled.bam"
    params:
        workdir     = lambda wildcards, output: dirname(output.pairA),
        downsampleN = config['kinship']['TKGWV2']['downsample-N'],
        seed        = TKGWV2_downsample_seed,
    log:       "logs/03-kinship/TKGWV2/TKGWV2_downsample_bam/{pairA}_{pairB}.log"
    benchmark: "benchmarks/03-kinship/TKGWV2/TKGWV2_downsample_bam/{pairA}_{pairB}.tsv"
    conda:     "../envs/TKGWV2.yml"
    shell: """
        root_dir=`pwd`                                                   # Keep current dir in memory.
        ln -sfrt {params.workdir} {input.pairs} >  $root_dir/{log}       # temporary symlink
        cd {params.workdir}                     >> $root_dir/{log}       # go to output workdir
        TK-helpers.py downsampleBam \
        --downsampleN {params.downsampleN} \
        --downsampleSeed {params.seed}          >> $root_dir/{log} 2>&1  # run TK-helpers
        find . -maxdepth 1 -type l -delete      >> $root_dir/{log}       # delete symlinks.
        
        mv {wildcards.pairA}.*_subsampled.bam $(basename {output.pairA}) >> $root_dir/{log} 2>&1
        mv {wildcards.pairB}.*_subsampled.bam $(basename {output.pairB}) >> $root_dir/{log} 2>&1
    """


def define_TKGWV2_input(wildcards):
    """
    Define the input for TKGWV2: either a downsampled bam file if the user
    requested it, or skip downsampling entirely.
    """
    if config['kinship']['TKGWV2']['downsample']:
        return rules.TKGWV2_downsample_bam.output
    else:
        return get_TKGWV2_input_bams(wildcards)


rule run_TKGWV2:
    """
    Run TKGWV2 on a single pair of individuals.
    @ TODO: Keep track of the 'frq' and 'tped' output files (TKGWV2 reorders pairA and pairB in lexical order.......)
    """
    input:
        bams          = define_TKGWV2_input,
        reference     = ReferenceGenome.get_path(),
        bed_targets   = "data/TKGWV2/genomeWideVariants_hg19/1000GP3_22M_noFixed_noChr.bed",
        plink_targets = multiext("data/TKGWV2/genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed", ".bed", ".bim", ".fam"),
        frequencies   = config['kinship']['TKGWV2']['target-frequencies']
    output:
        results = "results/03-kinship/TKGWV2/{pairA}_{pairB}/TKGWV2_Results.txt",
        #frq     = "results/03-kinship/TKGWV2/ped{gen}/{pairA}_{pairB}/commped{gen}_{A}____ped{gen}_{B}.frq",
        #tped    = "results/03-kinship/TKGWV2/ped{gen}/{pairA}_{pairB}/ped{gen}_{pairA}____ped{gen}_{pairB}.tped",
        peds    = [
            "results/03-kinship/TKGWV2/{pairA}_{pairB}/{pairA}.ped",
            "results/03-kinship/TKGWV2/{pairA}_{pairB}/{pairB}.ped"
        ],
        maps    = [
            "results/03-kinship/TKGWV2/{pairA}_{pairB}/{pairA}.map",
            "results/03-kinship/TKGWV2/{pairA}_{pairB}/{pairB}.map"
        ],
        pileups = [
            "results/03-kinship/TKGWV2/{pairA}_{pairB}/{pairA}.pileupsamtools.gwv.txt",
            "results/03-kinship/TKGWV2/{pairA}_{pairB}/{pairB}.pileupsamtools.gwv.txt"
        ],
    params:
        plink_basename = lambda wildcards, input: splitext(input.plink_targets[0])[0],
        min_MQ     = config["kinship"]["TKGWV2"]["min-MQ"],
        min_BQ     = config["kinship"]["TKGWV2"]["min-BQ"],
        min_depth  = config["kinship"]["TKGWV2"]["min-depth"],
        bam_ext    = lambda wildcards, input: basename(input.bams[0]).split(".",1)[1]
    log:       "logs/03-kinship/TKGWV2/run_TKGWV2/{pairA}_{pairB}.log"
    benchmark: "benchmarks/03-kinship/TKGWV2/run_TKGWV2/{pairA}_{pairB}.tsv"
    conda:     "../envs/TKGWV2.yml"
    shell: """
        base_dir=`pwd`                                     # Keep a record of the base directory
        cd $(dirname {output.results}) 2> $base_dir/{log}  # Go into the results directory

        # If the file is not present (i.e. no downsample has been made, create symlink.)
        for bam in {input.bams}; do
            [[ -f $(basename $bam) ]] || ln -sr $base_dir/$bam
        done 2>> $base_dir/{log}

        # Run TKGWV2
        TKGWV2.py bam2plink \
        --referenceGenome $base_dir/{input.reference} \
        --gwvList $base_dir/{input.bed_targets} \
        --gwvPlink $base_dir/{params.plink_basename} \
        --minMQ {params.min_MQ} \
        --minBQ {params.min_BQ} \
        --bamExtension .{params.bam_ext} \
        plink2tkrelated \
        --freqFile $base_dir/{input.frequencies} \
        --ignoreThresh {params.min_depth} \
        --verbose >> $base_dir/{log} 2>&1
    """

def define_TKGWV2_requested_dyads(wildcards):
    """
    Extract all relevant samples pairwise comparisons for TKGWV2
    """
    # Extract the relevant comparisons.
    pair_A, pair_B =  zip(*itertools.combinations(get_sample_names(), 2))
    relevant_comparisons = expand(
        rules.run_TKGWV2.output.results,
        zip,
        pairA=pair_A,
        pairB=pair_B
    )
    return relevant_comparisons



rule merge_TKGWV2_results:
    """
    Merge the pair-specific results of run_TKGWV2 into a single results file.
    (while removing header duplicates.)
    """
    input:
        results = define_TKGWV2_requested_dyads
    output:
        result = "results/03-kinship/TKGWV2/TKGWV2_Results.txt"
    log: "logs/03-kinship/TKGWV2/merge_TKGWV2_results.log"
    shell: """
        awk 'FNR>1 || NR==1' {input.results} > {output.result} 2> {log}
    """