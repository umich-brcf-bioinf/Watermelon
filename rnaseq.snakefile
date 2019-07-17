from __future__ import print_function, absolute_import, division

from collections import defaultdict, OrderedDict
from itertools import combinations, repeat
from shutil import copyfile
import csv
import os
import pandas as pd
import subprocess
import yaml

from scripts import rnaseq_snakefile_helper, config_validator, __version__

WAT_VER = __version__
WATERMELON_CONFIG_DIR = os.path.join(os.environ.get('WATERMELON_CONFIG_DIR', srcdir('config')), '')
WATERMELON_SCRIPTS_DIR = os.path.join(os.environ.get('WATERMELON_SCRIPTS_DIR', srcdir('scripts')), '')

_DIRS = config.get("dirs", {})
INPUT_DIR = os.path.join(_DIRS.get("input", "inputs"), "")
ALIGNMENT_DIR = os.path.join(_DIRS.get("alignment_output", "alignment_results"), "")
DIFFEX_DIR = os.path.join(_DIRS.get("diffex_output", "diffex_results"), "")
DELIVERABLES_DIR = os.path.join(_DIRS.get("deliverables_output", "deliverables"), "")
DESEQ2_DIR = os.path.join(DIFFEX_DIR, "deseq2", "")
BALLGOWN_DIR = os.path.join(DIFFEX_DIR, "ballgown", "")

CONFIGFILE_PATH = workflow.overwrite_configfile

SAMPLES_KEY = 'samples'

#Load in samplesheet
samplesheet = pd.read_csv(config["sample_description_file"]).set_index("sample", drop=True)

PHENOTYPES = (list(samplesheet.columns))

#Add sample info to config
#sample_phenotype_value_dict = samplesheet.to_dict(orient='index')
#TWS - for some reason, adding the above dict to the config breaks the Rscript functionality
#TWS - instead, adding only the sample names
config[SAMPLES_KEY] = list(samplesheet.index)

rnaseq_snakefile_helper.init_references(config["references"])

SAMPLE_READS = rnaseq_snakefile_helper.flattened_sample_reads(INPUT_DIR, config[SAMPLES_KEY])

#Perform config validation
schema_filepath = os.path.join(workflow.basedir, 'config', 'config_schema.yaml')
config_validation_exit_code = config_validator.main(CONFIGFILE_PATH, schema_filepath)
if config_validation_exit_code == 1:
    exit(1)

#Set up emailing functionality
def email(subject_prefix):
    msg = 'config file:\n{}\nlog file:\n{}'.format(workflow.overwrite_configfile, logger.get_logfile())
    email_config = config.get('email')
    if not email_config:
        print(subject_prefix, msg)
    else:
        command = "echo '{msg}' | mutt -s '{subject_prefix}{subject}' {to}".format(
                to=email_config['to'],
                subject_prefix=subject_prefix,
                subject=email_config['subject'],
                msg=msg,
                )
        shell(command)

onstart:
    email('Watermelon started: ')
onsuccess:
    email('Watermelon completed ok: ')
onerror:
    email('Watermelon ERROR: ')

#Set-up output targets
if 'fastq_screen' in config:
    FASTQ_SCREEN_CONFIG = config['fastq_screen']
    FASTQ_SCREEN_ALIGNMENT = [
            rnaseq_snakefile_helper.expand_sample_read_endedness(\
                ALIGNMENT_DIR + "03-fastq_screen/multi_species/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
            rnaseq_snakefile_helper.expand_sample_read_endedness(\
                ALIGNMENT_DIR + "03-fastq_screen/biotype/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS)
        ]
    FASTQ_SCREEN_DELIVERABLES = [
            rnaseq_snakefile_helper.expand_sample_read_endedness(
                DELIVERABLES_DIR + "alignment/fastq_screen/multi_species/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
            rnaseq_snakefile_helper.expand_sample_read_endedness(
                DELIVERABLES_DIR + "alignment/fastq_screen/biotype/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
        ]
else:
    FASTQ_SCREEN_CONFIG = defaultdict(str)
    FASTQ_SCREEN_ALIGNMENT = []
    FASTQ_SCREEN_DELIVERABLES = []

DESEQ2_CONTRAST_DICT = rnaseq_snakefile_helper.DESeq2_contrasts(config['diffex'])
DESeq2_ALL = [
    #deseq2_counts
    DIFFEX_DIR + 'deseq2/counts/txi_rsem_genes.rda',
    DIFFEX_DIR + 'deseq2/counts/count_data.rda',
    DIFFEX_DIR + 'deseq2/counts/raw_counts.txt',
    DIFFEX_DIR + 'deseq2/counts/depth_normalized_counts.txt',
    DIFFEX_DIR + 'deseq2/counts/rlog_normalized_counts.txt',
    #deseq2_contrasts
    rnaseq_snakefile_helper.expand_model_contrast_filenames(\
        DIFFEX_DIR + 'deseq2/gene_lists/{model_name}/{contrast}.txt',
        DESEQ2_CONTRAST_DICT),
    #deseq2_plots_by_phenotype
    expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/PCAplot_{dim}_top{ngenes}.pdf',
        phenotype = PHENOTYPES,
        dim = ['12','23'],
        ngenes = ['100','500']),
    expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/MDSplot_{dim}_top{ngenes}.pdf',
        phenotype = PHENOTYPES,
        dim = ['12','23'],
        ngenes = ['100','500']),
    expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/ScreePlot_top{ngenes}.pdf',
        phenotype = PHENOTYPES,
        ngenes = ['100','500']),
    expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/BoxPlot.pdf', phenotype = PHENOTYPES),
    expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/SampleHeatmap.pdf', phenotype = PHENOTYPES),
    expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/Heatmap_TopVar.pdf', phenotype = PHENOTYPES),
    expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/Heatmap_TopExp.pdf', phenotype = PHENOTYPES),
    #deseq2_comparison_plots
    rnaseq_snakefile_helper.expand_model_contrast_filenames(\
        DIFFEX_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.pdf',
        DESEQ2_CONTRAST_DICT),
    #deseq2_annotation
    rnaseq_snakefile_helper.expand_model_contrast_filenames(\
        DIFFEX_DIR + 'deseq2/annotated/{model_name}/{contrast}.annot.txt',
        DESEQ2_CONTRAST_DICT),
    #deseq2_excel
    rnaseq_snakefile_helper.expand_model_contrast_filenames(\
        DIFFEX_DIR + 'deseq2/excel/{model_name}/{contrast}.xlsx',
        DESEQ2_CONTRAST_DICT),
    #deseq2_summary
    DIFFEX_DIR + "deseq2/summary/deseq2_summary.txt",
    DIFFEX_DIR + "deseq2/summary/deseq2_summary.xlsx"
]

RSEM_ALL = [
    expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genes.results', sample=config[SAMPLES_KEY]),
    expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.isoforms.results', sample=config[SAMPLES_KEY]),
    expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genome.bam', sample=config[SAMPLES_KEY]),
    expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.transcript.bam', sample=config[SAMPLES_KEY])
]

DELIVERABLES = [
    #align_deliverables
    rnaseq_snakefile_helper.expand_sample_read_endedness(
        DELIVERABLES_DIR + "alignment/sequence_reads_fastqc/{sample}_trimmed_{read_endedness}_fastqc.html",
        SAMPLE_READS),
    expand(DELIVERABLES_DIR + "alignment/aligned_reads_fastqc/{sample}.genome_fastqc.html",
        sample=config["samples"]),
    DELIVERABLES_DIR + "alignment/alignment_qc.html",
    #deseq2_deliverables
    DELIVERABLES_DIR + "deseq2/gene_lists/deseq2_summary.txt" #TWS - why does this rule only have one provided output?
]



ALL = RSEM_ALL + DESeq2_ALL + FASTQ_SCREEN_ALIGNMENT + FASTQ_SCREEN_DELIVERABLES + DELIVERABLES


include: 'rules/align_concat_reads.smk'

if rnaseq_snakefile_helper.cutadapt_options(config["trimming_options"]):
    include: 'rules/align_cutadapt_SE.smk'
    include: 'rules/align_cutadapt_PE.smk'
else:
    include: 'rules/align_pseudotrim_SE.smk'
    include: 'rules/align_pseudotrim_PE.smk'

include: 'rules/align_fastq_screen_biotype.smk'
include: 'rules/align_fastq_screen_multi_species.smk'
include: 'rules/align_fastqc_trimmed_reads.smk'
include: 'rules/align_fastqc_align.smk'
include: 'rules/align_qc.smk'
include: 'rules/align_deliverables_alignment.smk'
include: 'rules/align_deliverables_fastq_screen.smk'

include: 'rules/deseq2_counts.smk'
include: 'rules/deseq2_init.smk'
include: 'rules/deseq2_contrasts.smk'
include: 'rules/deseq2_plots_by_phenotype.smk'
include: 'rules/deseq2_comparison_plots.smk'
include: 'rules/deseq2_annotation.smk'
include: 'rules/deseq2_excel.smk'
include: 'rules/deseq2_summary.smk'

include: 'rules/align_rsem_star.smk'
include: 'rules/rsem_star_genome_generate.smk'

include: 'rules/deliverables_deseq2.smk'
#include: 'rules/deliverables_run_info.smk'


rule all:
    input:
        *ALL
