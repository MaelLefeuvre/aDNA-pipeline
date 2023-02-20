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
                    basefile = "results/02-variant-calling/03-merged-reich/v52.2_1240K_public-merged-samples-autosomes-maf{maf}-reqsamples".format(maf=str(maf))
                else:
                    basefile = f"results/02-variant-calling/03-merged-reich/v52.2_1240K_public-merged-samples-autosomes-reqsamples"
            else:
                basefile = "results/02-variant-calling/02-pileupCaller/samples"
        case other:
            raise RuntimeError(f"Incorrect pileupCaller mode selected: '{other}'")

    return multiext(basefile, ".tped", ".tfam")


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
        tplink = READ_define_random_haploid_caller
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
    conda:     "../envs/READ.yml"
    shell: """
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
        tail -n+2 {input.raw} | awk '{{a[$1]+=$10}}END{{for (i in a) print i, a[i]}}' | sort | column -t > {output.overlap} 2> {log} 
    """

rule READ_get_relatedness_coefficient:
    input:
        means   = rules.run_READ.output.means,
        results = rules.run_READ.output.results
    output:
        relatedness = "results/03-kinship/READ/READ_relatedness_coefficients.txt"
    log: "logs/03-kinship/READ/READ_get_relatedness_coefficient.log"
    shell: """
    export LC_NUMERIC="en_US.UTF-8"
    sed 's/ /\t/g' {input.means} \
    | join -j1 -t$'\t' - {input.results} \
    | awk 'BEGIN{{FS="\t"; OFS=","}}{{print $1, 2*(1-$2), $6}}' \
    | column -s, -t \
    > {output.relatedness} 2> {log}
    """
