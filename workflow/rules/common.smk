import glob
import os

SAMP_NAMES = list(config['samples'].keys())

CTL_NAMES = list(config['control_samples'])

nctl = len(CTL_NAMES)

def get_index_name():
    genome = config["genome_fasta"]
    index = str(genome) + ".fai"
    return index

FORMAT = config['FORMAT']

ALIGNER = config['use_hisat3n']

NORMALIZE = config['normalize']

PAIRS = [1, 2]

nctl = len(CTL_NAMES)

def get_input_fastqs(wildcards):
    fastq_path = config["samples"][wildcards.sample]
    fastq_files = glob.glob(f"{fastq_path}/*.fastq*")
    return fastq_files

