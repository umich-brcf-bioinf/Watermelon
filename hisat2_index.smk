'''Build a HISAT2 index for a specified genome

e.g.:
First adjust config/genome_references.yaml to set
  references:hisat_2, fasta, gtf
Then,
$ snakemake --use-conda --snakefile hisat2_index.smk --cores=8 -p --config genome=ce10

Finally: Jazzhands!
'''

from os.path import join
configfile: 'config/genome_references.yaml'

GENOME = config[config['genome']]

OUTPUT_DIR = join(GENOME['references']['hisat2_index'], '')
TMP_EXON_FILE = OUTPUT_DIR + 'exon.tmp'
TMP_SS_FILE = OUTPUT_DIR + 'ss.tmp'

def build_spliced_commands(gtf):
    commands = ''
    if GENOME['hisat2_index_spliced']:
        commands = '''
            hisat2_extract_exons.py {gtf} > {exon_file}
            hisat2_extract_splice_sites.py {gtf} > {ss_file}''' \
                .format(gtf=gtf,
                        exon_file=TMP_EXON_FILE,
                        ss_file=TMP_SS_FILE)
    return commands

def build_spliced_flags(wildcards):
    flag = ''
    if GENOME['hisat2_index_spliced']:
        flag = '--exon {} --ss {}'.format(TMP_EXON_FILE, TMP_SS_FILE)
    return flag

rule all:
    input:
        OUTPUT_DIR + 'hisat2-inspect.summary',

rule align_create_hisat2_index:
    '''
We run a hisat2-inspect following the index creation to
a) validate the index and
b) gracefully allow the index builder to choose between a small index
   (suffix .ht2) and large index (suffix .ht21) without confusing snakemake.
'''
    input:
        gtf = GENOME['references']['gtf'],
        fasta = GENOME['references']['fasta'],
    output:
        OUTPUT_DIR + 'hisat2-inspect.summary',
    log:
        OUTPUT_DIR + '.log/align_create_hisat2_index.log',
    benchmark:
        OUTPUT_DIR + 'benchmarks/align_create_hisat2_index.benchmark.txt'
    conda:
        'rules/envs/align_hisat2_stringtie.yaml'
    params:
        spliced_commands = lambda wildcards, input: build_spliced_commands(input.gtf),
        spliced_flags = lambda wildcards: build_spliced_flags(wildcards),
    threads: 8
    shell:
        '''(
        {params.spliced_commands}
        hisat2-build -p {threads} \
            {params.spliced_flags} {input.fasta} \
            {OUTPUT_DIR}/genome
        hisat2-inspect --summary {OUTPUT_DIR}/genome \
            | tee {OUTPUT_DIR}/.hisat2-inspect.summary.tmp
        mv {OUTPUT_DIR}/.hisat2-inspect.summary.tmp {OUTPUT_DIR}/hisat2-inspect.summary
        ) 2>&1 | tee {log}'''
