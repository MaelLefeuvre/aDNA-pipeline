import yaml

configfile: "./config/config.yml"


# ----------------------------------------------------------------------------------------------- #
# ---- 00. Clip Adapters.

def assign_adapter_removal_seed(wildcards):
    """
    Return a user-defined seed from the config file if it was set. Else, fetch
    the randomly generated backup seed from our metadata file.
    """
    seed = config['preprocess']['trimming']['seed']
    if seed is None:
        with open(rules.meta.output.metadata) as f:
            metadata = yaml.load(f, Loader=yaml.loader.SafeLoader)
            seed     = metadata['seed']

    return seed


rule adapter_removal_pe:
    """
    Perform Adapter Trimming for Illumina Paired-End sequuencing data. Contrary to some other
    workflows, we don't output a combined fq file, to allow specific mapping using bwa aln.
    """
    input:
        forwd       = "original-data/samples/{sample}/{run}/{sample}_R1.fastq.gz",
        revrs       = "original-data/samples/{sample}/{run}/{sample}_R2.fastq.gz",
        metadata    = "results/meta/pipeline-metadata.yml",
    output:
        trimmed     = temp("results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.collapsed.gz"),
        truncated   = temp("results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.collapsed.truncated.gz"),
        discarded   = temp("results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.discarded.gz"),
        pair1       = temp("results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.pair1.truncated.gz"),
        pair2       = temp("results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.pair2.truncated.gz"),
        singleton   = temp("results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.singleton.truncated.gz")
    params:
        base_name   = "results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}",
        min_overlap = config["preprocess"]["trimming"]["min-overlap"],
        min_length  = config["preprocess"]["trimming"]["min-overlap"],
        min_quality = config["preprocess"]["trimming"]["min-quality"],
        quality_max = config["preprocess"]["trimming"]["qualitymax"],
        seed        = assign_adapter_removal_seed
    log: "logs/01-preprocess/01-adapter_removal/adapter_removal_pe/{sample}/{run}.log"
    benchmark: "benchmarks/01-preprocess/01-adapter_removal/adapter_removal_pe/{sample}/{run}-bench.tsv"
    conda: "../envs/adapterremoval-2.3.3.yml"
    priority: 4
    threads:  4
    shell: """
        AdapterRemoval \
        --threads {threads} \
        --file1 {input.forwd} \
        --file2 {input.revrs} \
        --basename {params.base_name} \
        --minlength {params.min_length} \
        --minquality {params.min_quality} \
        --qualitymax {params.quality_max} \
        --minadapteroverlap {params.min_overlap} \
        --seed {params.seed} \
        --collapse \
        --gzip 2> {log}
    """


rule adapter_removal_se:
    """
    Perform Adapter Trimming for Illumina Paired-End sequuencing data. Contrary to some other
    workflows, we don't output a combined fq file, to allow specific mapping using bwa aln.
    """
    input:
        fastq       = "original-data/samples/{sample}/{run}/{sample}_R1.fastq.gz",
        metadata    = "results/meta/pipeline-metadata.yml",
    output:
        truncated   = temp("results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.truncated.gz"),
        discarded   = temp("results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.discarded.gz"),
    params:
        base_name   = "results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}",
        min_overlap = config["preprocess"]["trimming"]["min-overlap"],
        min_length  = config["preprocess"]["trimming"]["min-overlap"],
        min_quality = config["preprocess"]["trimming"]["min-quality"],
        quality_max = config["preprocess"]["trimming"]["qualitymax"],
        seed        = assign_adapter_removal_seed
    log: "logs/01-preprocess/01-adapter_removal/adapter_removal_pe/{sample}/{run}.log"
    benchmark: "benchmarks/01-preprocess/01-adapter_removal/adapter_removal_pe/{sample}/{run}-bench.tsv"
    conda: "../envs/adapterremoval-2.3.3.yml"
    priority: 4
    threads:  4
    shell: """
        AdapterRemoval \
        --threads {threads} \
        --file1 {input.fastq} \
        --basename {params.base_name} \
        --minlength {params.min_length} \
        --minquality {params.min_quality} \
        --qualitymax {params.quality_max} \
        --minadapteroverlap {params.min_overlap} \
        --seed {params.seed} \
        --gzip 2> {log}
    """


# ----------------------------------------------------------------------------------------------- #
# ---- 01-A. Align to genome using BWA aln

rule sam_to_tmp_bam:
    """
    Convert a samfile to a **temporary** bam file.
    """
    input:
        sam = "{directory}/{file}.sam"
    output:
        bam = temp("{directory}/{file}.tmp-bam")
    threads: 1
    conda:   "../envs/samtools-1.15.yml"
    shell: """
        samtools view -@ {threads} -OBAM {input.sam} > {output.bam}
    """


rule bwa_aln:
    """
    Use bwa aln to prepare alignment and generate suffix array. Gives off very robust results, at 
    the cost of an impoverished runtime performance compared to bwa mem.
    See:    Oliva, A., Tobler, R., Llamas, B. and Souilmi, Y. (2021), Additional evaluations show that specific
            BWA-aln settings still outperform BWA-mem for ancient DNA data alignment. Ecol Evol, 11: 18743-18748. 
            https://doi.org/10.1002/ece3.8297

    extension: '.collapsed' || '.collapsed.truncated' || '.pair1.truncated' || '.pair2.truncated' || '.singleton.truncated'

    # Benchmarks max RSS:
    | depth | .collapsed | .collapsed.truncated | .pair{1,2} | .singleton.truncated |
    | ----- | ---------- | -------------------- | ---------- | -------------------- |
    | 0.01X | 6911       | 2964                 | 3211       | 2964                 |
    | 0.05X | 6997       | 2964                 | 3292       | 2964                 |
    """
    input:
        trimmed       = "results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.{extension}.gz",
        reference     = config["reference"],
        bwt           = rules.index_reference_genome.output.bwt
    output:
        sai           = temp(pipe("results/01-preprocess/02-align/{sample}/{run}/{sample}.bwaaln.{extension}.sai"))
    params:
        seed_length   = config['preprocess']['bwa']['bwa-aln']['seed-length'],
        max_open_gap  = config['preprocess']['bwa']['bwa-aln']['max-open-gap'],
        missing_prob  = config['preprocess']['bwa']['bwa-aln']['max-miss-prob'],
        max_seed_diff = config['preprocess']['bwa']['bwa-aln']['max-seed-diff']
    resources:
        mem_mb        = 7000
    log:       "logs/01-preprocess/02-align/bwa_aln/{extension}/{sample}/{run}.log"
    benchmark: "benchmarks/01-preprocess/02-align/bwa_aln/{extension}/{sample}/{run}-bench.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    #group:     "bwaaln"
    threads:   16
    priority:  50
    shell: """
        bwa aln \
        {input.reference} \
        {input.trimmed} \
        -t {threads} \
        -l {params.seed_length} \
        -n {params.missing_prob} \
        -k {params.max_seed_diff} \
        -o {params.max_open_gap}  > {output.sai} 2> {log}
    """

rule bwa_samse:
    """
    Perform single-end read mapping on a 'bwa aln' suffix array.
    extension: '.collapsed' || '.collapsed.truncated' || '.singleton.truncated'

    # Benchmarks max RSS:
    | depth | .collapsed | .collapsed.truncated | .singleton.truncated |
    | ----- | ---------- | -------------------- | -------------------- |
    | 0.01X | 4588       | 7.43                 | 7.35                 |
    | 0.05X | 4989       | 7.44                 | 18.23                |
    """
    input:
        trimmed   = "results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.{extension}.gz",
        sai       = rules.bwa_aln.output.sai,
        reference = config["reference"]
    output:
        sam       = temp(pipe("results/01-preprocess/02-align/{sample}/{run}/{sample}.bwaaln.{extension}.sam"))
    params:
        RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina' ## ADD IT LATER (-r argument)
    resources:
        mem_mb    = 5000
    log:       "logs/01-preprocess/02-align/bwa_samse/{sample}/{run}.{extension}.log"
    benchmark: "benchmarks/01-preprocess/02-align/bwa_samse/{extension}/{sample}/{run}-bench.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    #group:     "bwaaln"
    priority:  50
    threads:   1
    shell: """
        bwa samse -r '{params.RG}' {input.reference} {input.sai} {input.trimmed} > {output.sam} 2> {log}
    """

rule bwa_sampe:
    """
    Perform paired-end read mapping on a 'bwa aln' suffix array.
    extension: '.pair1.truncated' || '.pair2.truncated'

    # Benchmark max RSS:
    | depth | .paired |
    | ----- | ------- | 
    | 0.01X | 4478    |
    | 0.05X | 4618    |
    """
    input:
        pair1     = "results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.pair1.truncated.gz",
        pair2     = "results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.pair2.truncated.gz",
        sai1      = "results/01-preprocess/02-align/{sample}/{run}/{sample}.bwaaln.pair1.truncated.sai",
        sai2      = "results/01-preprocess/02-align/{sample}/{run}/{sample}.bwaaln.pair2.truncated.sai",
        reference = config["reference"]
    output:
        sam       = temp(pipe("results/01-preprocess/02-align/{sample}/{run}/{sample}.bwaaln.paired.sam"))
    params:
        RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina' 
    resources:
        mem_mb    = 5000
    log:       "logs/01-preprocess/02-align/bwa_sampe/{sample}/{run}.log"
    benchmark: "benchmarks/01-preprocess/02-align/bwa_sampe/{sample}/{run}.paired-bench.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    #group:     "bwaaln"
    priority: 50
    threads: 1
    shell: """
        bwa sampe -r '{params.RG}' {input.reference} {input.sai1} {input.sai2} {input.pair1} {input.pair2} > {output.sam} 2> {log}
    """


rule samtools_merge_aln:
    """
    Merge the outputs of bwa_samse & bwa_sampe for a given sample, and output a single merged bam file.

    # Benchmarks max RSS:
    | depth | .paired |
    | ----- | ------- | 
    | 0.01X | 122.77  | (VMS)
    | 0.05X | 12.81   |
    """
    input:
        paired_end = "results/01-preprocess/02-align/{sample}/{run}/{sample}.bwaaln.paired.tmp-bam",
        single_end = expand(
            "results/01-preprocess/02-align/{{sample}}/{{run}}/{{sample}}.bwaaln.{extension}.tmp-bam",
            extension=["collapsed", "collapsed.truncated", "singleton.truncated"]
        )
    output:
        merged    = "results/01-preprocess/02-align/{sample}/{run}/{sample}.bwaaln.merged.bam"
    params:
        RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina'
    log:       "logs/01-preprocess/02-align/samtools_merge/{sample}/{run}.log"
    benchmark: "benchmarks/01-preprocess/02-align/samtools_merge/{sample}/{run}-bench.tsv"
    conda:     "../envs/samtools-1.15.yml"
    #group:     "bwaaln"
    threads:   1
    priority:  60
    shell: """
        samtools merge -@ {threads} -o - {input.paired_end} {input.single_end} \
        | samtools addreplacerg -r '{params.RG}' -w -OBAM -o {output.merged} - 2> {log}
    """


# ----------------------------------------------------------------------------------------------- #
# ---- 01-B. Align to genome using BWA mem

"""
    Use bwa mem to perform read alignment and generate suffix array. This is much much faster than
    the traditional bwa aln + offers multithreading capability. but gives off less reliable results.
    Useful for prototyping, debugging, or when exactness isn't paramount.
    
    See:  Xu, W, Lin, Y, Zhao, K, et al. An efficient pipeline for ancient DNA mapping and recovery
          of endogenous ancient DNA from whole-genome sequencing data.
          https://doi.org/10.1002/ece3.7056 

    See:  Oliva, A., Tobler, R., Llamas, B. and Souilmi, Y. (2021), Additional evaluations show that specific
          BWA-aln settings still outperform BWA-mem for ancient DNA data alignment. 
          https://doi.org/10.1002/ece3.8297
"""

rule bwa_mem_se:
    """
    Perform a one pass single-end read mapping using bwa mem
    extension: '.collapsed' || '.collapsed.truncated' || '.singleton.truncated'
    """
    input:
        trimmed   = "results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.{extension}.gz",
        reference = config["reference"],
        bwt       = rules.index_reference_genome.output.bwt
    output:
        sam       = temp(pipe("results/01-preprocess/02-align/{sample}/{run}/{sample}.bwamem.{extension}.sam"))
    params:
        RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina'
    log:       "logs/01-preprocess/02-align/bwa_mem_se/{sample}/{run}.bwamem.{extension}.log"
    benchmark: "benchmarks/01-preprocess/02-align/bwa_mem_se/{extension}/{sample}/{run}.bwamem-bench.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    #group:     "bwamem"
    priority:  50
    threads:   4
    shell: """
        bwa mem -M -t {threads} -R \'{params.RG}\' -p {input.reference} {input.trimmed} > {output.sam} 2> {log}
    """


rule bwa_mem_pe:
    """
    Perform a one pass paired-end read mapping using bwa mem
    extension: '.pair1.truncated' || '.pair2.truncated'
    """
    input:
        pair1     = "results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.pair1.truncated.gz",
        pair2     = "results/01-preprocess/01-adapter_removal/{sample}/{run}/{sample}.pair2.truncated.gz",
        reference = config["reference"],
        bwt       = rules.index_reference_genome.output.bwt
    output:
        sam       = temp(pipe("results/01-preprocess/02-align/{sample}/{run}/{sample}.bwamem.paired.sam"))
    params:
        RG = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina'
    log:       "logs/01-preprocess/02-align/bwa_mem_pe/{sample}/{run}.bwamem.paired.log"
    benchmark: "benchmarks/01-preprocess/02-align/bwa_mem_pe/{sample}/{run}.paired-bench.tsv"
    conda:     "../envs/bwa-0.7.17.yml"
    #group:     "bwamem"
    priority:  50
    threads:   4
    shell: """
        bwa mem -M -t {threads} -R \'{params.RG}\' {input.reference} {input.pair1} {input.pair2} > {output.sam} 2> {log}
    """


rule samtools_merge_mem:
    """
    Merge the outputs of bwa_mem_se & bwa_mem_pe for a given sample, and output a single merged bam file.
    """
    input:
        paired_end = "results/01-preprocess/02-align/{sample}/{run}/{sample}.bwamem.paired.tmp-bam",
        single_end = expand(
            "results/01-preprocess/02-align/{{sample}}/{{run}}/{{sample}}.bwamem.{extension}.tmp-bam",
            extension=["collapsed", "collapsed.truncated", "singleton.truncated"]
        )
    output:
        merged     = "results/01-preprocess/02-align/{sample}/{run}/{sample}.bwamem.merged.bam"
    params:
        RG        = '@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina'
    log:       "logs/01-preprocess/02-align/samtools_merge/{sample}/{run}.log"
    benchmark: "benchmarks/01-preprocess/02-align/samtools_merge/{sample}/{run}.bwamem-bench.tsv"
    conda:     "../envs/samtools-1.15.yml"
    #group:     "bwamem"
    priority:  60
    threads:   4
    shell: """
        samtools merge -@ {threads} -o - {input.paired_end} {input.single_end} \
        | samtools addreplacerg -r '{params.RG}' -w -OBAM -o {output.merged} - 2> {log}
    """
