from requests.packages import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

configfile: "config/config.yml"
configfile: "config/netrules.yml"

class ReferenceGenome:

    _name = config["reference-genome"].lower()
    _dict = config["netrules"]["reference-genomes"]

    @staticmethod
    def get_name():
        return config["reference-genome"].lower()
    @classmethod
    def get_url(cls):
        return cls._dict[ReferenceGenome.get_name()]["url"]

    @classmethod
    def get_path(cls):
        return cls._dict[ReferenceGenome.get_name()]["path"]

    @classmethod
    def list_available_references(cls):
        for key in cls._dict.keys():
            print(f" - {key}")


def get_sample_names():
    # Order guaranteed to reflect the order of insertion, as of Python3.7+
    return [key for key in config["samples"]["raw-fastq"].keys()]

def get_illumina_protocol(wildcards):
    return config["samples"]["raw-fastq"][wildcards.sample][wildcards.run]["protocol"]

def get_requested_sample_runs(wildcards):
    return [key for key in config["samples"]["raw-fastq"][wildcards.sample].keys()]


_DEFAULT_FASTQ_R1 = "original-data/samples/{sample}/{run}/{sample}_R1.fastq.gz"
_DEFAULT_FASTQ_R2 = "original-data/samples/{sample}/{run}/{sample}_R2.fastq.gz"

def get_fastq_sample_path(wildcards):
    sample = config["samples"]["raw-fastq"][wildcards.sample][wildcards.run]
    if "path" not in sample or sample["path"] is None:
        sample["path"] = {
            "r1" : _DEFAULT_FASTQ_R1,
            "r2" : None if sample["protocol"] == "single" else _DEFAULT_FASTQ_R2
        } 
    return sample["path"]


def resolve_aadr_version():
    """
    Split a floating point version identifier into two major / minor wildcards
    """
    version = config["variant-calling"]["targets"]["aadr-version"]
    (major, minor) = str(version).split(".", 1)
    return {'major': major, 'minor': minor}

def parse_snp_targets_file(wildcards):
    requested_target_file = config["variant-calling"]["targets"]["custom-file"]
    if requested_target_file is not None:
        return requested_target_file

    return "data/Reich-dataset/1240K/v{major}.{minor}/v{major}.{minor}_1240K_public.snp".format(
        **resolve_aadr_version()
    )
