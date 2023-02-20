from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

from os.path import dirname

FTP  = FTPRemoteProvider(retry=config['netrules']['retries']) # Anonymous 
HTTP = HTTPRemoteProvider() # Anonymous 

configfile: "config/config.yml"

# ------------------------------------------------------------------------------------------------ #
# ---- 1000g-phase3 dataset

rule download_1000_genomes:
    """
    Download 1000genomes phase 3 SNPs
    http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
    """
    input:
        vcf = FTP.remote("ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
    output:
        vcf = "data/vcf/1000g-phase3/00-original/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    log: "logs/00-netrules/download_1000_genomes/download_1000_genomes-chr{chr}.log"
    shell: """
        mv {input.vcf} {output.vcf} 2> {log}
    """

rule fetch_samples_panel:
    """
    Download samples metadata from the 1000g FTP website
    """
    input:
        panel = FTP.remote("ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")
    output:
        panel = "data/vcf/1000g-phase3/samples-list/integrated_call_samples_v3.20130502.ALL.panel"
    log: "logs/00-netrules/fetch_samples_panel.log"
    shell: """
        mv {input.panel} {output.panel} > {log} 2>&1
    """

# ------------------------------------------------------------------------------------------------ #
# ---- Reich Lab's 1240K-Chip variant callset and panel

rule download_reich_1240K:
    """
    Download the 1240K dataset from Reich Lab's website.
    """
    input:
        tarball    = HTTP.remote("https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V{major}/V{major}.{minor}/SHARE/public.dir/v{major}.{minor}_1240K_public.tar")
    output:
        eigenstrat = multiext("data/Reich-dataset/1240K/v{major}.{minor}/v{major}.{minor}_1240K_public", ".snp", ".ind", ".geno")
    params:
        output_dir = lambda wildcards, output: dirname(output.eigenstrat[0])
    log: "logs/00-netrules/download_reich_1240K-v{major}.{minor}.log"
    shell:"""
        tar -xvf {input.tarball} -C {params.output_dir} > {log} 2>&1
    """

# ------------------------------------------------------------------------------------------------ #
# ---- Reference genomes

rule download_reference_genome:
    """
    Download a reference genome from a predefined ftp URL
    human_g1k_v37 is malformed... see: https://github.com/hammerlab/biokepi/issues/117
    """
    input:
        refgen = FTP.remote("ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/{reference}.fa.gz")
    output:
        refgen = "data/refgen/GRCh37/{reference}.fa.gz"
    log: "logs/00-netrules/download_reference_genome/{reference}.log"
    shell: """
        mv {input.refgen} {output.refgen}
    """


# ------------------------------------------------------------------------------------------------ #
# ---- Genetic maps

rule download_HapMapII_recombination_map:
    """
    Download the 2011 HapMapII recombination map from ncbi.
    """
    input:
        tarball = HTTP.remote("http://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz")
    output:
        map     = expand("data/recombination-maps/HapMapII_GRCh37/genetic_map_GRCh37_chr{chr}.txt", chr=range(1, 23)),
        exclude = temp(expand("data/recombination-maps/HapMapII_GRCh37/genetic_map_GRCh37_chr{chr}.txt", chr=["X", "X_par1", "X_par2"])),
        readme  = "data/recombination-maps/HapMapII_GRCh37/README.txt"
    params:
        output_dir = lambda wildcards, output: dirname(output.map[0])
    log: "logs/00-netrules/download_HapMapII_recombination_map.log"
    shell: """
        tar -xvzf {input.tarball} -C {params.output_dir} 2>  {log}
        rm {input.tarball}                               2>> {log}
    """

# ------------------------------------------------------------------------------------------------ #
# ---- Miscellaneous

rule download_TKGWV2_support_files:
    """
    Download Daniel Fernandes' 22M SNP panel from his public google drive.
    """
    output:
        support_files = expand("{directory}/{dataset}", 
            directory = config['netrules']['TKGWV2']['support-files-dir'],
            dataset = [
                "1240K/1000GP3_EUR_1240K.frq",
                "genomeWideVariants_hg19/1000GP3_22M_noFixed_noChr.bed",
                "genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed.bed",
                "genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed.bim",
                "genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed.fam"
            ]
        )
    params:
        url = config['netrules']['TKGWV2']['support-files-url'],
        output_dir = config['netrules']['TKGWV2']['support-files-dir']
    conda: "../envs/gdown-4.6.0.yml"
    log: "logs/04-kinship/TKGWV2/TKGWV2_download_support_files.log"
    shell: """
        gdown "{params.url}" -O {params.output_dir} --folder > {log} 2>&1
    """