localrules: symlink_KINgaroo_input_bams



rule symlink_KINgaroo_input_bams:
    input:
        bamlist = get_pileup_input_bams,
    output:
        bamlist  = "results/03-kinship/KIN/samples.bamlist",
        linkdir  = directory("results/03-kinship/KIN/symbams/")
    resources:
        cores    = lambda w, threads: threads
    log:      "logs/03-kinship/KIN/symlink_KINgaroo_input_bams.log"
    conda:    "../envs/coreutils-9.1.yml"
    threads: 1
    shell: """
        mkdir -p {output.linkdir} && ln -srft {output.linkdir} {input.bamlist} 2>  {log}
        basename -a {input.bamlist} | sed 's/.bam$//' > {output.bamlist}       2>> {log}
    """

def parse_KINgaroo_optargs(wildcards):
    correct_contamination = False
    diversity_parameter = config['kinship']['KIN']['diversity-parameter']
    noisy_windows       = config['kinship']['KIN']['noisy-windows']
    optargs = ""
    if diversity_parameter is not None:
        optargs += f"--diversity_parameter_p_0 {diversity_parameter} "
    if noisy_windows is not None:
        optargs += f"--noisy_wins {noisy_windows} "

    return optargs



rule run_KINgaroo:
    """
    Run KINgaroo on a set of pedigree individuals (generation-wise.)
    # @ TODO: Add contamination estimate. 

    # Benchmarks: 
    | depth | max h:m:s | max_rss |
    | ----- | --------- | ------- |
    | 0.01X |           |         |
    | 0.05X | 0:11:35   | 4922.52 |
    | 0.10X | 0:45:00   | 3539    |
    """
    input:
        bamlist        = rules.symlink_KINgaroo_input_bams.output.bamlist,
        linkdir        = rules.symlink_KINgaroo_input_bams.output.linkdir,
        targets        = lambda w: os.path.splitext(parse_snp_targets_file(w))[0] + ".ucscbed",
    output:
        kingaroo_dir   =  directory("results/03-kinship/KIN/kingaroo"),
        bedfiles       = directory("results/03-kinship/KIN/kingaroo/bedfiles"),
        hapProbs       = directory("results/03-kinship/KIN/kingaroo/hapProbs"),
        hbd_results    = directory("results/03-kinship/KIN/kingaroo/hbd_results"),
        hmm_parameters = directory("results/03-kinship/KIN/kingaroo/hmm_parameters"),
        splitbams      = directory("results/03-kinship/KIN/kingaroo/splitbams"),
        goodpairs      = "results/03-kinship/KIN/kingaroo/goodpairs.csv",
        overlap        = "results/03-kinship/KIN/kingaroo/identical_overlap.csv",
        diffs_hmm      = "results/03-kinship/KIN/kingaroo/input_diffs_hmm.csv",
        hbd_hmm_diffs  = "results/03-kinship/KIN/kingaroo/input_hbd_hmm_diffs.csv",
        hbd_hmm_total  = "results/03-kinship/KIN/kingaroo/input_hbd_hmm_total.csv",
        total_hmm      = "results/03-kinship/KIN/kingaroo/input_total_hmm.csv"
    params:
        interval       = config['kinship']['KIN']['interval'],
        threshold      = config['kinship']['KIN']['p0-threshold'],
        optargs        = parse_KINgaroo_optargs
    resources:
        #runtime        = 60,
        mem_mb         = 6000,
    resources:
        cores          = lambda w, threads: threads
    log:       "logs/03-kinship/KIN/run_KINgaroo.log"
    benchmark: "benchmarks/03-kinship/KIN/run_KINgaroo.tsv"
    conda:      "../envs/kin-3.1.3.yml"
    threads:    16
    shell: """
        CWD=`pwd`
        cd {output.kingaroo_dir}
        KINgaroo \
        --cores {threads} \
        --bamfiles_location $CWD/{input.linkdir} \
        --target_location $CWD/{input.bamlist} \
        --bedfile $CWD/{input.targets} \
        --contam_parameter 0 \
        --interval {params.interval} \
        --threshold {params.threshold} \
        {params.optargs} \
        > $CWD/{log} 2>&1
    """


def parse_KIN_optargs(wildcards):
    optargs = ""
    roh_threshold       = config['kinship']['KIN']['roh-threshold']
    diversity_parameter = config['kinship']['KIN']['diversity-parameter']
    if roh_threshold is not None:
        optargs += f"--threshold {roh_threshold} "
    if diversity_parameter is not None:
        optargs += f"--diversity_parameter_p_0 {diversity_parameter} "

    return optargs


rule run_KIN:
    """
    Perform generation-wise kinship estimation using KIN, and the output of KINgaroo.
    See: Popli, D., Peyrégne, S. & Peter, B.M. KIN: a method to infer relatedness from low-coverage
         ancient DNA. Genome Biol 24, 10 (2023). https://doi.org/10.1186/s13059-023-02847-7

    Repo: https://github.com/DivyaratanPopli/Kinship_Inference.git

    # @TODO: Add contamination estimate. 

    # Benchmarks: 
    | depth | max h:m:s | max_rss |
    | ----- | --------- | ------- |
    | 0.01X |           |         |
    | 0.05X | 0:00:08   | 2464    |
    | 0.10X | 0:00:17   | 2121    |
    """
    input:
        bamlist      = rules.symlink_KINgaroo_input_bams.output.bamlist,
        linkdir      = rules.symlink_KINgaroo_input_bams.output.linkdir,
        kingaroo_dir = rules.run_KINgaroo.output.kingaroo_dir
    output:
        kin_results  = directory("results/03-kinship/KIN/KIN-results"),
        gamma_files  = directory("results/03-kinship/KIN/KIN-results/gammafiles"),
        lik_files    = directory("results/03-kinship/KIN/KIN-results/likfiles"),
        res_files    = directory("results/03-kinship/KIN/KIN-results/resfiles"),
        IBDadded     = "results/03-kinship/KIN/KIN-results/IBDadded.csv",
        KIN_results  = "results/03-kinship/KIN/KIN-results/KIN_results.csv",
        likelihoods  = "results/03-kinship/KIN/KIN-results/relatable_allLikelihoods_fil0.csv"
    params:
        workdir      = lambda w, output: os.path.dirname(output.kin_results),
        interval     = config['kinship']['KIN']['interval'],
        optargs      = parse_KIN_optargs
    resources:
        #runtime      = 10,
        mem_mb       = 3000,
        cores        = lambda w, threads: threads
    log:       "logs/03-kinship/KIN/run_KIN.log"
    benchmark: "benchmarks/03-kinship/KIN/run_KIN.tsv"
    conda:      "../envs/kin-3.1.3.yml"
    threads:    16
    shell: """
        CWD=`pwd`; cd {params.workdir}
        KIN \
        --cores {threads} \
        --input_location $CWD/{input.kingaroo_dir}/ \
        --output_location $CWD/{output.kin_results}/ \
        --target_location $CWD/{input.bamlist} \
        --interval {params.interval} \
        {params.optargs} \
        > $CWD/{log} 2>&1
    """