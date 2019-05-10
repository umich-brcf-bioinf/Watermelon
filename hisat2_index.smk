'''Build a HISAT2 index for all genomes in config

e.g.:
First adjust config/genome_references.yaml to set
  references:hisat_2, fasta, gtf
Then,
$ snakemake --use-conda --snakefile hisat2_index.smk --cores=8 -p

Finally: Jazzhands near the Whole Foods hot bar!
'''

from os.path import join
from argparse import Namespace
configfile: 'config/genome_references.yaml'

def get_genome_config(wildcards):
    genome = config[wildcards.genome]
    genome_config = Namespace()
    genome_config.output_dir = join(genome['references']['hisat2_index'], '')
    genome_config.gtf = genome['references']['gtf']
    genome_config.fasta = genome['references']['fasta']
    genome_config.spliced = genome['hisat2_index_spliced']
    genome_config.exon_file = genome_config.output_dir + 'exon.tmp'
    genome_config.ss_file = genome_config.output_dir + 'ss.tmp'
    genome_config.spliced_commands = _build_spliced_commands(genome_config)
    genome_config.spliced_flags = _build_spliced_flags(genome_config)
    return genome_config

def _build_spliced_commands(genome_config):
    commands = ''
    if genome_config.spliced:
        commands = '''
            hisat2_extract_exons.py {gtf} > {exon_file} && hisat2_extract_splice_sites.py {gtf} > {ss_file}''' \
                .format(gtf=genome_config.gtf,
                        exon_file=genome_config.exon_file,
                        ss_file=genome_config.ss_file)
    return commands

def _build_spliced_flags(genome_config):
    flag = ''
    if genome_config.spliced:
        flag = '--exon {} --ss {}'\
            .format(genome_config.exon_file, genome_config.ss_file)
    return flag

get_output = lambda genome: join(config[genome]['references']['hisat2_index'], genome+'.summary')
OUTPUTS = [get_output(genome) for genome in config]

rule all:
    input:
        OUTPUTS,

rule align_create_hisat2_index:
    '''
We run a hisat2-inspect following the index creation to
a) validate the index and
b) gracefully allow the index builder to choose between a small index
   (suffix .ht2) and large index (suffix .ht21) without confusing Snakemake.
'''
    input:
        gtf = lambda x: get_genome_config(x).gtf,
        fasta = lambda x: get_genome_config(x).fasta,
    output:
        '{output_dir}/{genome}.summary',
    log:
        '{output_dir}/.log/align_create_hisat2_index.{genome}.log',
    benchmark:
        '{output_dir}/benchmarks/align_create_hisat2_index.{genome}.benchmark.txt'
    conda:
        'rules/envs/align_hisat2_stringtie.yaml'
    params:
        output_dir = lambda x: get_genome_config(x).output_dir,
        spliced_commands = lambda x: get_genome_config(x).spliced_commands,
        spliced_flags = lambda x: get_genome_config(x).spliced_flags,
    threads: 8
    shell:
        '''(
        {params.spliced_commands}
        hisat2-build -p {threads} \
            {params.spliced_flags} {input.fasta} \
            {params.output_dir}/genome
        hisat2-inspect --summary {params.output_dir}/genome \
            | tee {output}.tmp
        mv {output}.tmp {output}
        ) 2>&1 | tee {log}'''
