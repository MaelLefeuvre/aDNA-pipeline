def get_sample_names(wildcards):
    return [key for key in config['samples'].keys()]
    #return [list(sample.keys())[0] for sample in config['samples']]


def get_illumina_protocol(wildcards):
    return config['samples'][wildcards.sample][wildcards.run]['protocol']

def get_requested_sample_runs(wildcards):
    return [key for key in config['samples'][wildcards.sample].keys()]