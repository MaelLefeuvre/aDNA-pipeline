# ---- Set config variables
configfile: "./config/config.yml"

# ------------------------------------------------------------------------------------------------------------------- #
# ---- Generic / Utility rules 

rule samtools_faidx:
    """
    Generate a `.fai` fasta index using samtools faidx
    """
    input:
        fasta = "{directory}/{fasta}.fa"
    output:
        fai   = "{directory}/{fasta}.fa.fai"
    conda: "../envs/samtools-1.15.yml"
    shell: """
        samtools faidx {input.fasta} --fai-idx {output.fai}
    """

# ------------------------------------------------------------------------------------------------------------------- #
# ---- 01. Decompress the reference genome.
#          @TODO: This could be a generic rule.

rule decompress_reference_genome:
    """
    Decrompress a reference file from .fa.gz -> .fa 
    """
    input:
        # refgen = rules.download_reference_genome.output.refgen
        refgen = "data/refgen/GRCh37/{reference}.fa.gz"
    output:
        refgen = protected("data/refgen/GRCh37/{reference}.fa")
    log: "logs/00-preprocess-reference/decompress_reference_genome/{reference}.log"
    shell: """
        gunzip -v {input.refgen} > {log} 2>&1
    """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 02. Index the reference genome.

rule index_reference_genome:
    """
    Generate the Burrows-Wheeler transform on the reference genome.
    """
    input:
        refgen = config["reference"]
    output:
        amb = config["reference"] + ".amb",
        ann = config["reference"] + ".ann",
        bwt = config["reference"] + ".bwt",
        pac = config["reference"] + ".pac",
        sa  = config["reference"] + ".sa",
    conda: "../envs/bwa-0.7.17.yml"
    log: "logs/00-preprocess-reference/index_reference_genome.log"
    shell: """
        bwa index {input.refgen} >2 {log}
    """

rule picard_create_sequence_dictionary:
    input:
        reference = config['reference']
    output:
        dictionary = os.path.splitext(config['reference'])[0] + ".dict"
    params:
        tmpdir = config["tempdir"]
    log:   "logs/00-preprocess-reference/picard_create_sequence_dictionary.log"
    conda: "../envs/picard-2.27.4.yml"
    shell: """
        picard CreateSequenceDictionary --TMP_DIR {params.tmpdir} --REFERENCE {input.reference} --OUTPUT {output.dictionary} > {log} 2>&1
    """


# ------------------------------------------------------------------------------------------------------------------- #
# ---- 03. Split the reference genome on a per-chromosome basis

rule split_reference_genome:
    """
    Split the reference genome according to chromosome.
    """
    input:
        refgen   = config["reference"]
    output:
        splitted = expand("data/refgen/GRCh37/splitted/{chr}.fasta", chr=range(1,23))
    log: "logs/00-preprocess-reference/split_reference_genome.log"
    shell: """
        curr_wd=`pwd`
        cd $(dirname {output.splitted[0]}) 2> {log}
        csplit -s -z $curr_wd/{input.refgen} '/>/' '{{*}}'
        for i in xx* ; do                                  \
            n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ;         \
            mv "$i" "$n.fasta" ;                           \
        done 2>> {log}
    """