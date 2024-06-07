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
    return [key for key in config["samples"].keys()]

def get_illumina_protocol(wildcards):
    return config["samples"][wildcards.sample][wildcards.run]["protocol"]

def get_requested_sample_runs(wildcards):
    return [key for key in config["samples"][wildcards.sample].keys()]

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
