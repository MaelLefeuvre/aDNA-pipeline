# ---- Set config variables
configfile: "./config/config.yml"

# ------------------------------------------------------------------------------------------------------------------- #
# ---- Generic / Utility rules 

rule samtools_faidx:
    """
    Generate a `.fai` fasta index using samtools faidx
    """
    input:
        fasta = ReferenceGenome.get_path()
    output:
        fai   = ReferenceGenome.get_path() + ".fai"
    log: f"logs/00-preprocess-reference/samtools_faidx/{ReferenceGenome.get_path()}.log"
    conda: "../envs/samtools-1.15.yml"
    shell: """
        samtools faidx {input.fasta} --fai-idx {output.fai} > {log} 2>&1
    """

rule decompress_reference_genome:
    """
    Decompress a reference file from .fa.gz -> .fa
    """
    input:
        refgen = ReferenceGenome.get_path() + ".gz"
    output:
        refgen = ReferenceGenome.get_path()
    log: "logs/00-preprocess-reference/decompress_reference_genome.log"
    conda: "../envs/coreutils-9.1.yml"
    shell: """
        gunzip -v {input.refgen} > {log} 2>&1 || test $? = 2
    """

rule index_reference_genome:
    """
    Index a referecnce genome using Burrows-Wheeler transform.
    """
    input:
        refgen = ReferenceGenome.get_path()
    output:
        amb = ReferenceGenome.get_path() + ".amb",
        ann = ReferenceGenome.get_path() + ".ann",
        bwt = ReferenceGenome.get_path() + ".bwt",
        pac = ReferenceGenome.get_path() + ".pac",
        sa  = ReferenceGenome.get_path() + ".sa",
    conda: "../envs/bwa-0.7.17.yml"
    log: "logs/00-preprocess-reference/index_reference_genome.log"
    shell: """
        bwa index {input.refgen} > {log} 2>&1
    """

rule picard_create_sequence_dictionary:
    input:
        reference = ReferenceGenome.get_path()  
    output:
        dictionary = os.path.splitext(ReferenceGenome.get_path())[0] + ".dict"
    log:   "logs/00-preprocess-reference/picard_create_sequence_dictionary.log"
    conda: "../envs/picard-2.27.4.yml"
    shell: """
        picard CreateSequenceDictionary --TMP_DIR {resources.tmpdir} --REFERENCE {input.reference} --OUTPUT {output.dictionary} > {log} 2>&1
    """

rule split_reference_genome:
    """
    Split a reference genome on a per-chromosome basis
    """
    input:
        refgen = ReferenceGenome.get_path()
    output:
        splitted  = expand(f"{os.path.dirname(ReferenceGenome.get_path())}/splitted/{{chr}}.fasta", chr=range(1, 23))
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