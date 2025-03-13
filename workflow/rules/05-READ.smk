localrules: READ_get_pairwise_snp_overlap, READ_get_relatedness_coefficient
# ------------------------------------------------------------------------------------------------------------------- #

def READ_define_random_haploid_caller(wildcards):
    """
    Define the input for READ, based on which variant-calling method was 
    requested by the user.
    """

    variant_caller = config['variant-calling']['caller']
    use_proxies    = config['kinship']['READ']['add-proxies']
    maf            = config['kinship']['READ']['maf']

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

    return multiext(basefile, ".tped", ".tfam")

def get_coverage_threshold(wildcards):
    threshold = config['kinship']['READ']['coverage-threshold']
    if str(threshold).endswith('X'):
        threshold = threshold[:-1]
    
    return threshold


def get_coverage_column(wildcards):
    threshold = str(config['kinship']['READ']['coverage-threshold'])
    return 3 if threshold.endswith('X') else 2


rule filter_low_coverage_samples:
    input:
        meta     = rules.meta.output,
        coverage = rules.get_bamlist_panel_coverage.output.coverage,
        tplink   = READ_define_random_haploid_caller
    output:
        excluded_samples = "results/03-kinship/READ/excluded-samples.txt",
        tplink           = multiext("results/03-kinship/READ/READ-input-filtered", ".tped", ".tfam")
    params:
        coverage_threshold = get_coverage_threshold,
        coverage_column    = get_coverage_column,
        sample_names       = get_sample_names(),
        input_basename     = lambda wildcards, input: splitext(input.tplink[0])[0],
        output_basename    = lambda wildcards, output: splitext(output.tplink[0])[0],
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
            --recode transpose tab \
            --out {params.output_basename} \
            > {log} 2>&1
        else
            ln -sr {input.tplink[0]} {output.tplink[0]} > {log} 2>&1;
            ln -sr {input.tplink[1]} {output.tplink[1]} > {log} 2>&1;
        fi
    """




def get_READ_norm_value(wildcards):
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
# ---- Serial READ
rule run_READ:
    """
    Run the original, serialized version of READ
    See: Monroy Kuhn JM, Jakobsson M, Günther T (2018) Estimating genetic kin relationships in prehistoric populations.
         PLoS ONE 13(4): e0195491. https://doi.org/10.1371/journal.pone.0195491
         https://bitbucket.org/tguenther/read.git
    """
    input:
        tplink  = rules.filter_low_coverage_samples.output.tplink
    output:
        results = "results/03-kinship/READ/READ_results",
        means   = "results/03-kinship/READ/meansP0_AncientDNA_normalized",
        raw     = "results/03-kinship/READ/READ_output_ordered",
        plot    = "results/03-kinship/READ/READ_results_plot.pdf",
        tempraw = temp("results/03-kinship/READ/Read_intermediate_output"),
    params:
        window_size = config["kinship"]["READ"]["window-size"],
        norm_method = config["kinship"]["READ"]["norm-method"],
        norm_value  = get_READ_norm_value,
        basename    = lambda wildcards, input: splitext(input.tplink[0])[0],
    log:       "logs/03-kinship/READ/run_READ.log"
    benchmark: "benchmarks/03-kinship/READ/run_READ.tsv"
    conda:     "../envs/READ-1.0.yml"
    shell: r"""
        cwd=$(pwd)
        cd $(dirname {output.results})               2>  $cwd/{log}
        ln -srf $(which READscript.R) READscript.R   2>> $cwd/{log}
        touch meansP0_AncientDNA_normalized          2>> $cwd/{log}
        python2 $(which READ.py) $cwd/{params.basename} {params.norm_method} {params.norm_value} \
        --window_size {params.window_size} >> $cwd/{log} 2>&1

        # Remove symlink
        find . -type l -exec rm {{}} \;
    """

rule READ_get_pairwise_snp_overlap:
    input:
        raw = rules.run_READ.output.raw
    output:
        overlap = "results/03-kinship/READ/READ_overlapping_snps.txt"
    log: "logs/03-kinship/READ/READ_get_pairwise_snp_overlap.log"
    shell: """
        tail -n+2 {input.raw} \
        | awk '{{a[$1]+=$10}}END{{for (i in a) print i, a[i]}}' \
        | sort \
        | column -t > {output.overlap} 2> {log} 
    """

rule READ_get_relatedness_coefficient:
    input:
        means   = rules.run_READ.output.means,
        results = rules.run_READ.output.results
    output:
        relatedness = "results/03-kinship/READ/READ_relatedness_coefficients.txt"
    log: "logs/03-kinship/READ/READ_get_relatedness_coefficient.log"
    shell: r"""
        export LC_NUMERIC="en_US.UTF-8"
        sed 's/ /\t/g' {input.means} \
        | join -j1 -t$'\t' - {input.results} \
        | awk 'BEGIN{{FS="\t"; OFS=","}}{{print $1, 2*(1-$2), $6}}' \
        | column -s, -t \
        > {output.relatedness} 2> {log}
    """
