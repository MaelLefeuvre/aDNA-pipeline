from os.path import splitext
localrules: render_parfile_eigenstrat_to_packedped, make_plink_keepfile
# ------------------------------------------------------------------------------------------------------------------- #
# ---- Generic / Utility rules and functions.

def assign_plink_optargs(wildcards):
    """
    Return a user-defined seed from the config file if it was set. Else, fetch
    the randomly generated backup seed from our metadata file.
    """
    seed = config['variant-calling']['plink']['seed']
    if seed is None:
        with checkpoints.meta.get().output.metadata.open() as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']

    return f"--allow-no-sex --keep-allele-order --seed {seed}"

rule eigenstrat_to_UCSC_BED:
    """
    Convert an eigensoft's .snp file to a generic bed file 
    (termed 'ucscbed', to avoid collisions with plink's file format).
    """
    input:
        snp   = "{directory}/{eigenstrat}.snp"
    output:
        bed   = "{directory}/{eigenstrat}.ucscbed"
    resources:
        cores = lambda w, threads: threads
    log:       "logs/generics/{directory}/eigenstrat_to_UCSC_BED-{eigenstrat}.log"
    benchmark: "benchmarks/generics/{directory}/eigenstrat_to_UCSC_BED-{eigenstrat}.tsv"
    conda:     "../envs/coreutils-9.1.yml"
    threads:   1
    shell: """
        awk 'BEGIN{{OFS="\t"}}{{print $2, $4-1, $4, $5, $6}}' {input.snp} > {output.bed}
    """


rule plink_bfile_to_tped:
    """
    Convert a binarized PLINK fileset into a human readable transposed set.
    """
    input:
        metadata = rules.meta.output,
        bfile    = multiext("{directory}/{file}", ".bed", ".bim", ".fam"),
    output:
        tplink   = multiext("{directory}/{file}", ".tped", ".tfam") 
    params:
        basename = "{directory}/{file}",
        optargs  = assign_plink_optargs
    log:       "logs/{directory}/{file}/plink_bfile_to_tped.log"
    benchmark: "benchmarks/{directory}/{file}/plink_bfile_to_tped.tsv"
    conda:     "../envs/plink-1.9.yml"
    threads:   16
    shell: """
        plink --threads {threads} --bfile {params.basename} {params.optargs} --recode transpose tab --out {params.basename} > {log} 2>&1
    """


# ------------------------------------------------------------------------------------------------------------------- #

rule render_parfile_eigenstrat_to_packedped:
    """
    Generate an eigenstrat parameter file for convertf conversion file.
    """
    input:
        template = "resources/templates/par.eigenstrat.packedped.template"
    output:
        parfile  = "{directory}/par.{filestem}.EIGENSTRAT.PED"
    params:
        basename = "{directory}/{filestem}"
    log: "logs/{directory}/{filestem}/render_parfile_eigenstrat_to_packedped.log"
    template_engine: "jinja2"


rule convertf_eigenstrat_to_plink:
    """
    Convert an eigenstrat [snp|geno|ind] fileset into a binarized plink fileset
    [bed|bim|fam] using convertf.
    """
    input:
        eigenstrat = multiext("{directory}/{filestem}", ".snp", ".geno", ".ind"),
        parfile    = rules.render_parfile_eigenstrat_to_packedped.output.parfile,
    output:
        plink      = multiext("{directory}/{filestem}", ".bed", ".bim", ".fam")
    params:
        basename   = "{directory}/{filestem}"
    log:       "logs/{directory}/{filestem}/convertf_eigenstrat_to_plink.log"
    benchmark: "benchmarks/{directory}/{filestem}/convertf_eigenstrat_to_plink.tsv"
    conda:     "../envs/eigensoft-8.0.0.yml"
    shell: """
        convertf -p {input.parfile} > {log} 2>&1 
    """

rule plink_merge_reich:
    """
    Merge our samples variant callset with the 1240K compendium dataset.
    Note that this operation uses plink, and not mergeit.
    """
    input:
        metadata         = rules.meta.output,
        reich            = multiext(splitext(rules.download_reich_1240K.output.eigenstrat[0])[0].format(**resolve_aadr_version()), ".bed", ".bim", ".fam"),
        samples          = multiext("results/02-variant-calling/02-pileupCaller/samples", ".bed", ".bim", ".fam"),
    output:
        merged           = multiext("results/02-variant-calling/03-merged-reich/v{major}.{minor}_1240K_public-merged-samples".format(**resolve_aadr_version()), ".bed", ".bim", ".fam")
    params:
        reich_basename   = lambda wildcards, input: splitext(input.reich[0])[0],
        samples_basename = lambda wildcards, input: splitext(input.samples[0])[0],
        merged_basename  = lambda wildcards, output: splitext(output.merged[0])[0],
        optargs          = assign_plink_optargs
    log:       "logs/03-kinship/READ/plink_merge_reich.log"
    benchmark: "benchmarks/03-kinship/READ/plink_merge_reich.tsv"
    conda:     "../envs/plink-1.9.yml"
    threads: 16
    shell: """
        plink --threads {threads} --bfile {params.reich_basename} {params.optargs} --bmerge {params.samples_basename} --merge-mode 2 --make-bed --out {params.merged_basename} > {log} 2>&1
    """

rule plink_filter_autosomes:
    """
    Filter out non-autossomal chromosomes from a generic bfile set [bed|bim|fam].
    """
    input:
        metadata        = rules.meta.output,
        bfile           = multiext("{directory}/{filestem}", ".bed", ".bim", ".fam")
    output:
        bfile           = multiext("{directory}/{filestem}-autosomes", ".bed", ".bim", ".fam")
    params:
        input_basename  = "{directory}/{filestem}",
        output_basename = "{directory}/{filestem}-autosomes",
        optargs         = assign_plink_optargs
    log:       "logs/{directory}/{filestem}/plink_filter_autosomes.log"
    benchmark: "benchmarks/{directory}/{filestem}/plink_filter_autosomes.tsv"
    conda:     "../envs/plink-1.9.yml"
    threads:   16
    shell: """
        plink  --threads {threads} --bfile {params.input_basename} {params.optargs} --chr 1-22 --make-bed --out {params.output_basename} > {log} 2>&1
    """

rule plink_filter_maf:
    """
    Perform minor allele filtration on a generic bfile set [bed|bim|fam].
    """
    input:
        metadata        = rules.meta.output,
        bfile           = multiext("{directory}/{filestem}", ".bed", ".bim", ".fam")
    output:
        bfile           = multiext("{directory}/{filestem}-maf{maf}", ".bed", ".bim", ".fam")
    params:
        input_basename  = lambda wildcards, input : splitext(input.bfile[0] )[0],
        output_basename = lambda wildcards, output: splitext(output.bfile[0])[0],
        optargs         = assign_plink_optargs
    wildcard_constraints:
        maf = "[.0-9]+"
    log:       "logs/{directory}/{filestem}/plink_filter_maf-maf{maf}.log"
    benchmark: "benchmarks/{directory}/{filestem}/plink_filter_maf-maf{maf}.tsv"
    conda:     "../envs/plink-1.9.yml"
    threads:   16
    shell: """
        plink  --threads {threads} --bfile {params.input_basename} {params.optargs} --maf {wildcards.maf} --make-bed --out {params.output_basename} > {log} 2>&1
    """

def get_requested_samples(wildcards):
    """
    Return a list of samples to extract from the 1240K-merge dataset.
    - user-defined samples are always extracted
    - the list of 'proxies' are only extracted if the user requested it
      using the 'add-proxies' config file directive.
    
    @TODO: implement ability to work directly w/ bam samples. 
    """
    proxies = get_sample_names()
    match config['kinship']['READ']['add-proxies']:
        case None:
            return proxies
        case "Reich":    
            proxies += config['kinship']['READ']['proxies']
            return proxies
        case other:
            raise RuntimeError(
                f"Attempting to extract proxies from Reich dataset, \
                but config file has 'use-proxies' set to {other}"
            )
        

rule make_plink_keepfile:
    """
    Generate a list of samples to 'extract' from the reich dataset
    """
    input:
        fam      = "{directory}/{filestem}.fam"
    output:
        keepfile = "{directory}/{filestem}-keep-samples.txt"
    params:
        samples  = get_requested_samples
    log:       "logs/{directory}/{filestem}/plink_extract_requested_samples.log"
    shell: """
        grep -f <(echo {params.samples} | tr " " "\n") {input.fam} | cut -f1,2 -d" " > {output.keepfile}
    """



rule plink_extract_requested_samples:
    """
    Extract a list of samples from a generic bfile set into a separate
    [bed|bim|fam] fileset
    """
    input:
        metadata        = rules.meta.output,
        keepfile        = rules.make_plink_keepfile.output.keepfile,
        bfile           = multiext("{directory}/{filestem}", ".bed", ".bim", ".fam")
    output:
        bfile           = multiext("{directory}/{filestem}-reqsamples", ".bed", ".bim", ".fam")
    params:
        input_basename  = lambda wildcards, input : splitext(input.bfile[0] )[0],
        output_basename = lambda wildcards, output: splitext(output.bfile[0])[0],
        optargs         = assign_plink_optargs
    log:       "logs/{directory}/{filestem}/plink_extract_samples.log"
    benchmark: "benchmarks/{directory}/{filestem}/plink_extract_samples.tsv"
    conda:     "../envs/plink-1.9.yml"
    threads:   16
    shell: """
        plink --threads {threads} --bfile {params.input_basename} {params.optargs} --keep {input.keepfile} --make-bed --out {params.output_basename} > {log} 2>&1
    """
