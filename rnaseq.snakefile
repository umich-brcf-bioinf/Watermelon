from __future__ import print_function, absolute_import, division

from collections import defaultdict, OrderedDict
from itertools import combinations, repeat
from shutil import copyfile
import csv
import logging
import os
import pandas as pd
import subprocess
import yaml

from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler
from scripts import rnaseq_snakefile_helper, config_validator, __version__

WAT_VER = __version__
WATERMELON_CONFIG_DIR = os.path.join(os.environ.get('WATERMELON_CONFIG_DIR', srcdir('config')), '')
WATERMELON_SCRIPTS_DIR = os.path.join(os.environ.get('WATERMELON_SCRIPTS_DIR', srcdir('scripts')), '')


_DIRS = config.get("dirs", {})
INPUT_DIR = os.path.join(_DIRS.get("input", "inputs"), "")
ALIGNMENT_DIR = os.path.join(_DIRS.get("alignment_output", "alignment_results"), "")
DIFFEX_DIR = os.path.join(_DIRS.get("diffex_output", "diffex_results"), "")
DELIVERABLES_DIR = os.path.join(_DIRS.get("deliverables_output", "deliverables"), "")

CONFIGFILE_PATH = workflow.overwrite_configfile
WORKFLOW_BASEDIR = workflow.basedir
CONFIG_SCHEMA_PATH = os.path.join(workflow.basedir, 'config', 'config_schema.yaml')

SAMPLES_KEY = 'samples'

#Load in samplesheet
samplesheet = pd.read_csv(config["sample_description_file"], comment='#').set_index("sample", drop=True)

PHENOTYPES = (list(samplesheet.columns))

#Add sample info to config
#sample_phenotype_value_dict = samplesheet.to_dict(orient='index')
#TWS - for some reason, adding the above dict to the config breaks the Rscript functionality
#TWS - instead, adding only the sample names
config[SAMPLES_KEY] = list(samplesheet.index)

rnaseq_snakefile_helper.init_references(config["references"])

SAMPLE_READS = rnaseq_snakefile_helper.flattened_sample_reads(INPUT_DIR, config[SAMPLES_KEY])


#Set up emailing functionality
def email(subject_prefix, msg="", attachment=""):

    email_config = config.get('email')
    if not email_config:
        print(subject_prefix, msg)
    else:
        command = "echo '{msg}' | mutt -s '{subject_prefix}{subject}' {attachment} {to}".format(
                to=email_config['to'],
                subject_prefix=subject_prefix,
                subject=email_config['subject'],
                attachment=attachment,
                msg=msg,
                )
        shell(command)

#function to call config_validator
def validate_config(config_fp, schema_fp):
    config_validation_exit_code = config_validator.main(config_fp, schema_fp)
    if config_validation_exit_code:
        exit(config_validation_exit_code)



#Perform config validation if dryrun
if set(['-n', '--dryrun']).intersection(set(sys.argv)):
    validate_config(CONFIGFILE_PATH, CONFIG_SCHEMA_PATH)

onstart:
    #Perform config validation before starting
    validate_config(CONFIGFILE_PATH, CONFIG_SCHEMA_PATH)

    #Set up error-only logfile with same basename as snakemake logfile
    #Doing this within onstart because logger.get_logfile evaluates to None for some reason during workflow execution
    ERR_ONLY_LOGFILE = os.path.splitext(logger.get_logfile())[0] + "-errors.log"
    #Create filehandler to create error-only logfile for error-level logging events
    fh = logging.FileHandler(ERR_ONLY_LOGFILE, delay=True, mode='w')
    fh.setLevel(logging.ERROR)
    #Add this filehandler to existing snakemake logger
    logger.logger.addHandler(fh)

    #Setup handler to trigger email upon creation of error-only logfile
    class MyWatchdogHandler(PatternMatchingEventHandler):
        def on_created(self, event):
            super(MyWatchdogHandler, self).on_created(event)
            message = "Watermelon has detected at least one error. Attaching nascent error-only log-file"
            attach_str = "-a " + ERR_ONLY_LOGFILE + " --"
            email('Watermelon error notification: ', msg=message, attachment=attach_str)
    event_handler = MyWatchdogHandler(patterns=[ERR_ONLY_LOGFILE])
    observer = Observer()
    observer.schedule(event_handler, path=os.path.dirname(ERR_ONLY_LOGFILE), recursive=False)
    observer.start()

    email('Watermelon started: ')
onsuccess:
    message = 'config file:\n{}\nlog file:\n{}'.format(workflow.overwrite_configfile, logger.get_logfile())
    email('Watermelon completed ok: ', msg=message)
onerror:
    message = "Watermelon completed with errors. Full log file attached"
    attach_str = "-a " + logger.get_logfile() + " --"
    email('Watermelon completed with errors: ', msg=message, attachment = attach_str)



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

DESEQ2_CONTRAST_DICT = rnaseq_snakefile_helper.diffex_contrasts(config['diffex'])
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
    expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.transcript.bam', sample=config[SAMPLES_KEY]),
    ALIGNMENT_DIR + "07-qc/alignment_qc.html",
    expand(ALIGNMENT_DIR + '05-combine_counts/{feature}_{metric}.txt',
        feature=['gene','isoform'],
        metric=['expected_count', 'FPKM', 'TPM'])
]

DELIVERABLES = [
    #align deliverables
    rnaseq_snakefile_helper.expand_sample_read_endedness(
        DELIVERABLES_DIR + "alignment/sequence_reads_fastqc/{sample}_trimmed_{read_endedness}_fastqc.html",
        SAMPLE_READS),
    expand(DELIVERABLES_DIR + "alignment/aligned_reads_fastqc/{sample}.genome_fastqc.html",
        sample=config["samples"]),
    expand(DELIVERABLES_DIR + "counts/gene_{type}.txt",
        type=['expected_count', 'FPKM', 'TPM']),
    DELIVERABLES_DIR + "alignment/alignment_qc.html",
    #deseq2 deliverables
    expand(DELIVERABLES_DIR + 'deseq2/counts/{name}.txt',
        name=['raw_counts', 'depth_normalized_counts', 'rlog_normalized_counts']),
    rnaseq_snakefile_helper.expand_model_contrast_filenames(\
        DELIVERABLES_DIR + 'deseq2/gene_lists/{model_name}/{contrast}.txt',
        DESEQ2_CONTRAST_DICT),
    rnaseq_snakefile_helper.expand_model_contrast_filenames(\
        DELIVERABLES_DIR + "deseq2/excel/{model_name}/{contrast}.xlsx",
        DESEQ2_CONTRAST_DICT),
    expand(DELIVERABLES_DIR + 'deseq2/plots/by_phenotype/{phenotype}/PCAplot_{dim}_top{ngenes}.pdf',
            phenotype = PHENOTYPES,
            dim = ['12','23'],
            ngenes = ['100','500']),
    expand(DELIVERABLES_DIR + 'deseq2/plots/by_phenotype/{phenotype}/ScreePlot_top{ngenes}.pdf',
            phenotype = PHENOTYPES,
            ngenes = ['100','500']),
    expand(DELIVERABLES_DIR + 'deseq2/plots/by_phenotype/{phenotype}/{plotType}.pdf',
        phenotype = PHENOTYPES, plotType = ['BoxPlot', 'SampleHeatmap', 'Heatmap_TopVar', 'Heatmap_TopExp']),
    rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DELIVERABLES_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.pdf',
            DESEQ2_CONTRAST_DICT),
    DELIVERABLES_DIR + "deseq2/summary/deseq2_summary.txt",
    DELIVERABLES_DIR + "deseq2/summary/deseq2_summary.xlsx",
    #run info deliverables
    DELIVERABLES_DIR + "run_info/env_software_versions.yaml",
    DELIVERABLES_DIR + "run_info/" + os.path.basename(CONFIGFILE_PATH),
    DELIVERABLES_DIR + "run_info/" + os.path.basename(config['sample_description_file'])
]



ALL = RSEM_ALL + DESeq2_ALL + FASTQ_SCREEN_ALIGNMENT + FASTQ_SCREEN_DELIVERABLES + DELIVERABLES


include: 'rules/align_concat_reads.smk'
include: 'rules/align_standardize_gz.smk'

if "trimming_options" in config:
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
include: 'rules/align_rsem_star.smk'
include: 'rules/align_rsem_star_genome_generate.smk'
include: 'rules/align_rsem_star_combined_count_matrices.smk'

include: 'rules/deseq2_counts.smk'
include: 'rules/deseq2_init.smk'
include: 'rules/deseq2_contrasts.smk'
include: 'rules/deseq2_plots_by_phenotype.smk'
include: 'rules/deseq2_comparison_plots.smk'
include: 'rules/deseq2_annotation.smk'
include: 'rules/deseq2_excel.smk'
include: 'rules/deseq2_summary.smk'

include: 'rules/deliverables_deseq2.smk'
include: 'rules/deliverables_run_info.smk'


rule all:
    input:
        *ALL
