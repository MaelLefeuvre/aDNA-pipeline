import yaml
from os.path import dirname
# ------------------------------------------------------------------------------------------------------------------- #

rule GRUPS_generate_fst_set:
    """
    Generate a .fst and .fst.frq dataset from a set of VCF files. Allows for
    generally faster IO processing when performing multiple runs.
    """
    input:
        data    = expand(rules.download_1000_genomes.output.vcf, chr=range(1,23)),
        panel   = rules.fetch_samples_panel.output.panel
    output:
        fst     = protected(expand(
            "data/grups/fst/g1k-phase3-v5/{ped_pop}-{cont_pop}/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes-{ped_pop}-{cont_pop}.{ext}",
            ped_pop  = config["kinship"]["GRUPS"]["pedigree-pop"],
            cont_pop = config["kinship"]["GRUPS"]["contam-pop"],
            chrom    = range(1,23),
            ext      = ["fst", "fst.frq"]
        ))
    params:
        pedigree_pop = config["kinship"]["GRUPS"]["pedigree-pop"],
        contam_pop   = config["kinship"]["GRUPS"]["contam-pop"]
    log:       "logs/03-kinship/GRUPS/GRUPS_generate_fst_set/{params.pedigree_pop}-{params.contam_pop}-GRUPS_generate_fst_set.log"
    benchmark: "benchmarks/03-kinship/GRUPS/GRUPS_generate_fst_set/{params.pedigree_pop}-{params.contam_pop}-GRUPS_generate_fst_set.tsv"
    conda:     "../envs/grups-rs.yml"
    threads:   22
    shell: """
        grups fst \
        --vcf-dir $(dirname {input.data} | uniq) \
        --output-dir $(dirname {output.fst} | uniq) \
        --pop-subset {params.pedigree_pop} {params.contam_pop} \
        --panel {input.panel} \
        --threads {threads} \
        --verbose > {log} 2>&1
    """

def format_grups_optargs(wildcards):
    """
    Return a user-defined seed from the config file if it was set. Else, fetch
    the randomly generated backup seed from our metadata file.
    """
    optargs        = ""
    seed           = config['kinship']['GRUPS']['seed']
    seq_error_rate = config['kinship']['GRUPS']['seq-error-rate']
    if seed is None:
        with open(rules.meta.output.metadata) as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']

    optargs += f"--seed {seed} "

    if seq_error_rate is not None:
        optargs += f"--seq-error-rate {seq_error_rate}"
    
    return optargs


rule run_GRUPS:
    """
    Run grups on an entire batch of samples.
    @ TODO: VCF mode is not yet implemented within the pipeline.
            A simple function should do the trick.
    """
    input:
        pileup       = rules.samtools_pileup.output.pileup,
        data         = rules.GRUPS_generate_fst_set.output.fst,
        bamlist      = rules.generate_bam_list.output.bamlist,
        recomb_map   = rules.download_HapMapII_recombination_map.output.map,
        targets      = config["variant-calling"]["targets"],
        metadata     = "results/meta/pipeline-metadata.yml"
    output:
        output_dir   = directory("results/03-kinship/GRUPS"),
        results      = multiext("results/03-kinship/GRUPS/samples", ".pwd", ".result")
    params:
        data_dir     = lambda wildcards, input: dirname(input.data[0]),
        recomb_dir   = lambda wildcards, input: dirname(input.recomb_map[0]),
        pedigree     = config["kinship"]["GRUPS"]["pedigree"],
        pedigree_pop = config["kinship"]["GRUPS"]["pedigree-pop"],
        contam_pop   = config["kinship"]["GRUPS"]["contam-pop"],
        reps         = config["kinship"]["GRUPS"]["reps"],
        mode         = config["kinship"]["GRUPS"]["mode"],
        min_depth    = config["kinship"]["GRUPS"]["min-depth"],
        min_quality  = config["kinship"]["GRUPS"]["min-qual"], 
        maf          = config["kinship"]["GRUPS"]["maf"],
        optargs      = format_grups_optargs
    log:       "logs/03-kinship/GRUPS/run_GRUPS.log"
    benchmark: "benchmarks/03-kinship/GRUPS/run_GRUPS.tsv"
    conda:     "../envs/grups-rs.yml"
    shell: """
        grups pedigree-sims \
        --pileup {input.pileup} \
        --data-dir {params.data_dir} \
        --recomb-dir {params.recomb_dir} \
        --pedigree {params.pedigree} \
        --pedigree-pop {params.pedigree_pop} \
        --contam-pop {params.contam_pop} \
        --min-depth {params.min_depth} \
        --samples 0-$(($(wc -l {input.bamlist} | cut -f1 -d" ")-1)) \
        --sample-names $(cat {input.bamlist} | rev | cut -d'/' -f1 | rev | cut -d'.' -f1 | tr '\n' ' ') \
        --reps {params.reps} \
        --mode {params.mode} \
        --output-dir {output.output_dir} \
        --maf {params.maf} \
        --min-qual {params.min_quality} \
        {params.optargs} \
        --print-blocks \
        --ignore-dels \
        --verbose > {log} 2>&1
    """