def READv2_define_random_haploid_caller(wildcards):
    """
    Define the input for READ, based on which variant-calling method was 
    requested by the user.
    """

    variant_caller = config['variant-calling']['caller']
    use_proxies    = config['kinship']['READv2']['add-proxies']
    maf            = config['kinship']['READv2']['maf']

    if use_proxies == "Reich" and variant_caller != 'pileupCaller':
        raise RuntimeError("[ERROR]: Cannot use proxy individuals from Reich compendium if variant-caller is not set to pileupCaller.")

    match config['variant-calling']['caller']:
        case "ANGSD":
            basefile = "results/02-variant-calling/02-ANGSD/samples"
        case "pileupCaller":
            if use_proxies == "Reich":
                if maf is not None:
                    basefile = "results/02-variant-calling/03-merged-reich/v{major}.{minor}_1240K_public-merged-samples-autosomes-maf{maf}-reqsamples".format(
                        maf=str(maf),
                        **resolve_aadr_version()
                        )
                else:
                    basefile = "results/02-variant-calling/03-merged-reich/v{major}.{minor}_1240K_public-merged-samples-autosomes-reqsamples".format(
                        **resolve_aadr_version()
                    )
            else:
                basefile = "results/02-variant-calling/02-pileupCaller/samples"
        case other:
            raise RuntimeError(f"Incorrect pileupCaller mode selected: '{other}'")

    return multiext(basefile, ".bed", ".bim", ".fam")

def get_coverage_threshold(wildcards):
    threshold = config['kinship']['READ']['coverage-threshold']
    if str(threshold).endswith('X'):
        threshold = threshold[:-1]
    
    return threshold


def get_coverage_column(wildcards):
    threshold = str(config['kinship']['READ']['coverage-threshold'])
    return 3 if threshold.endswith('X') else 2


rule READv2_filter_low_coverage_samples:
    input:
        meta     = rules.meta.output,
        coverage = rules.get_bamlist_panel_coverage.output.coverage,
        tplink   = READv2_define_random_haploid_caller
    output:
        excluded_samples = "results/03-kinship/READv2/excluded-samples.txt",
        bfile           = multiext("results/03-kinship/READv2/READv2-input-filtered", ".bed", ".bim", ".fam")
    params:
        coverage_threshold = get_coverage_threshold,
        coverage_column    = get_coverage_column,
        sample_names       = get_sample_names(),
        input_basename     = lambda wildcards, input: splitext(input.tplink[0])[0],
        output_basename    = lambda wildcards, output: splitext(output.bfile[0])[0],
        sample_pop_name   = config['variant-calling']['sample-pop-name'],
        optargs           = assign_plink_optargs

    log:       "logs/03-kinship/READ/filter_low_coverage_samples.log"
    benchmark: "benchmarks/03-kinship/READ/filter_low_coverage_samples.tsv"
    conda:     "../envs/plink-1.9.yml"
    threads: 1
    shell: """
        # Get a list of individuals not meeting the threshold
        # Note that 'cut -f1 -d. ' implies sample_names containing '.' characters is undefined behavior...
        # trailing grep commands with test $? = 1 ensures non-match is not interpreted as an error.
        LC_ALL="C" awk '${params.coverage_column}<{params.coverage_threshold}' {input.coverage} \
        | {{ grep -P "$(echo {params.sample_names} | tr ' ' '|')" || test $? = 1; }} \
        | cut -f1 -d' ' \
        | {{ grep -o '[^/]*$' || test $? = 1; }} \
        | cut -f1 -d"." \
        | sort -u \
        > {output.excluded_samples} 2>> {log}

        # Add population names
        sed -i -E 's/^(.+)/{params.sample_pop_name} \\1/' {output.excluded_samples} >> {log} 2>&1

        if [ -s {output.excluded_samples} ]; then
            plink --threads {threads} {params.optargs} \
            --tfile {params.input_basename} \
            --remove {output.excluded_samples} \
            --make-bed \
            --out {params.output_basename} \
            > {log} 2>&1
        else
            plink --threads {threads} {params.optargs} \
            --bfile {params.input_basename} \
            --make-bed \
            --out {params.output_basename} \
            > {log} 2>&1
        fi
    """




def get_READv2_norm_value(wildcards):
    """
    Inject the user-defined normalization value within READ's command line 
    arguments (if the user requested it). Or add a placeholder "-" character
    in place if the user wishes to compute that value from the cohort.
    """
    if config["kinship"]["READ"]["norm-method"] == "value":
        return config["kinship"]["READ"]["norm-value"]
    else:
        return "-"

# ------------------------------------------------------------------------------------------------------------------- #
# ---- READ

def parse_READv2_norm_method(wildcards):
    """
    Get user-defined normalization value within READv2' command line argument
    If the user requested it.
    """
    norm_method = config['kinship']['READv2']['norm-method']
    norm_value  = config['kinship']['READv2']['norm-value']
    match norm_method:
        case None:
            return ""
        case "median" | "mean" | "max":
            return f"--norm_method {norm_method}"
        case "value":
            if norm_value is None:
                raise RuntimeError("[READv2]: requested a set normalisation value, but none provided")
            return f"--norm_method value --norm_value {config['kinship']['READv2']['norm-value']}"
        case other:
            raise RuntimeError(f"[READv2]: Invalid norm-method: {other}")

def parse_READv2_window_estimate(wildcards):
    """
    Add user-defined window estimate strategy within READv2' command line argument
    if the user requested it.
    """
    output_optargs            = ""
    window_size               = config['kinship']['READv2']['window-size']
    window_estimate_requested = config['kinship']['READv2']['window-est']
    if window_estimate_requested:
        output_optargs += "--window_est"
        if window_size is not None:
            try:
                output_optargs += f" --window_size {int(window_size)}"
            except:
                raise RuntimeError(f"[READv2]: Invalid window-size provided: {window_size}")
    return output_optargs

def parse_READv2_alternate_thresholds(wildcards):
    """
    Inject the user-requested window threshold stragtegy within READv2's command line argument

    - When False: [0.96875,0.90625,0.8125,0.625]
    - When True:  [1-1/(2**4.5),1-1/(2**3.5),1-1/(2**2.5),1-1/(2**1.5)]
    """
    alternate_threshold_requested = config['kinship']['READv2']['2pow']
    return "--2pow" if alternate_threshold_requested else ""


def count_comparisons(wildcards, include_self = False):
    with open(rules.generate_bam_list.output.bamlist) as f:
        n = len([line for line in f.readlines()])
        combinations = (n * (n-1)) / 2
    return combinations + n if include_self else combinations

# ------------------------------------------------------------------------------------------------ #

rule run_READv2:
    """
    Run the updated version of READ
    Citation: Alaçamlı, E., Naidoo, T., Güler, M.N. et al. READv2: advanced and user-friendly
              detection of biological relatedness in archaeogenomics. Genome Biol 25, 216 (2024).
              https://doi.org/10.1186/s13059-024-03350-3
    Github:   https://github.com/GuntherLab/READv2
    """
    input:
        bfile   = rules.READv2_filter_low_coverage_samples.output.bfile
    output:
        results = "results/03-kinship/READv2/Read_Results.tsv",
        means   = "results/03-kinship/READv2/meansP0_AncientDNA_normalized_READv2",
        #plot    = "results/03-kinship/READv2/READ.pdf",
    params:
        norm_method   = parse_READv2_norm_method,
        window_size   = parse_READv2_window_estimate,
        alt_threshold = parse_READv2_alternate_thresholds,
        basename      = lambda wildcards, input: splitext(input.bfile[0])[0],
    log:       "logs/03-kinship/READv2/run_READv2.log"
    benchmark: "benchmarks/03-kinship/READv2/run_READv2.tsv"
    conda:     "../envs/READv2.yml"
    shell: r"""
        cwd=$(pwd)
        cd $(dirname {output.results}) > $cwd/{log} 2>&1
        READ2 --input $cwd/{params.basename} \
        {params.norm_method} {params.window_size} {params.alt_threshold} \
        >> $cwd/{log} 2>&1
    """
