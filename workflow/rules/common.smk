SAMP_NAMES = list(config['samples'].keys())

CTL_NAMES = list(config['control_samples'])

nctl = len(CTL_NAMES)

def get_index_name():
    genome = config["genome_fasta"]
    index = str(genome) + ".fai"
    return index

FORMAT = config['FORMAT']

ALIGNER = config['use_hisat3n']

SAMP_NAMES = list(config['samples'].keys())

CPU_num = config['cpus']

PAIRS = [1, 2]

CTL_NAMES = list(config['control_samples'])

nctl = len(CTL_NAMES)

def get_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]