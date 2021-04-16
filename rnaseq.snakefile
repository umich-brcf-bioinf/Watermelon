from __future__ import print_function, absolute_import, division

from collections import defaultdict, OrderedDict
from itertools import chain, combinations, repeat
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
WORKFLOW_BASEDIR = workflow.basedir
WATERMELON_CONFIG_DIR = os.path.join(WORKFLOW_BASEDIR, 'config', '')
WATERMELON_SCRIPTS_DIR = os.path.join(WORKFLOW_BASEDIR, 'scripts', '')
# WATERMELON_CONFIG_DIR = os.path.join(os.environ.get('WATERMELON_CONFIG_DIR', srcdir('config')), '')
# WATERMELON_SCRIPTS_DIR = os.path.join(os.environ.get('WATERMELON_SCRIPTS_DIR', srcdir('scripts')), '')

#TWS TODO: Consider refactoring - Should we really return "" if key doesn't exist? Added requirements to schema for now
_DIRS = config.get("dirs", {})
INPUT_DIR = os.path.join(_DIRS.get("input", "inputs"), "")
ALIGNMENT_DIR = os.path.join(_DIRS.get("alignment_output", "alignment_results"), "")
DIFFEX_DIR = os.path.join(_DIRS.get("diffex_output", "diffex_results"), "")
DELIVERABLES_DIR = os.path.join(_DIRS.get("deliverables_output", "deliverables"), "")
REPORT_DIR = os.path.join(_DIRS.get("report_output", "report"), "")

JOB_LOG_DIR = os.path.join(os.getcwd(), "job_logs", "")

CLUSTER_LOG_DIR = os.path.join(os.getcwd(), "cluster_logs")

#CONFIGFILE_PATH = workflow.overwrite_configfile

CONFIG_SCHEMA_PATH = os.path.join(WATERMELON_CONFIG_DIR, 'config_schema.yaml')

SAMPLES_KEY = 'samples'

#Load in samplesheet
samplesheet = pd.read_csv(config["samplesheet"], comment='#', dtype='object').set_index("sample", drop=True)

PHENOTYPES = (list(samplesheet.columns))

#Add sample info to config
#sample_phenotype_value_dict = samplesheet.to_dict(orient='index')
#TWS - for some reason, adding the above dict to the config breaks the Rscript functionality
#TWS - instead, adding only the sample names
config[SAMPLES_KEY] = list(samplesheet.index)

rnaseq_snakefile_helper.init_references(config["references"])

if config.get('capture_regex'):
    capture_regex = config['capture_regex']
else:
    capture_regex = r'.*_R(\d+)[_0-9]*\.fastq.*'


if not config.get('count_matrix'):
    INPUT_MANAGER = rnaseq_snakefile_helper.InputFastqManager(input_dir=INPUT_DIR, capture_regex=capture_regex)


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


#Get config file path
if '--configfile' in list(sys.argv):
    cfg_idx = sys.argv.index('--configfile') + 1
elif '--configfiles' in list(sys.argv): # Cluster environment seems to use configfiles plural form under the hood
    cfg_idx = sys.argv.index('--configfiles') + 1
else:
    msg = '--configfile not specified in {}'.format(sys.argv)
    raise ValueError(msg)

CONFIGFILE_PATH = sys.argv[cfg_idx]

#Perform config validation if dryrun
if set(['-n', '--dryrun']).intersection(set(sys.argv)) and not 'skip_validation' in config:
    validate_config(CONFIGFILE_PATH, CONFIG_SCHEMA_PATH)

#Check for cluster_logs folder if run in cluster environment
if set(['-c', '--cluster']).intersection(set(sys.argv)):
    if not os.path.exists(CLUSTER_LOG_DIR):
        msg = "Cluster log folder {} not found. Creating it before running in a cluster environment."
        logger.logger.info(msg.format(CLUSTER_LOG_DIR))
        os.mkdir(CLUSTER_LOG_DIR)
elif '--profile' in sys.argv:
    profile_idx = sys.argv.index('--profile') + 1
    profile = sys.argv[profile_idx]
    profile_config = snakemake.get_profile_file(profile, "config.yaml")
    with open(profile_config) as fh:
        profile_dict = yaml.load(fh, Loader=yaml.SafeLoader)
        if 'cluster' in profile_dict and not os.path.exists(CLUSTER_LOG_DIR):
            msg = "Cluster log folder {} not found. Creating it before running in a cluster environment."
            logger.logger.info(msg.format(CLUSTER_LOG_DIR))
            os.mkdir(CLUSTER_LOG_DIR)


onstart:
    #Perform config validation before starting
    if not 'skip_validation' in config:
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
    message = 'config file:\n{}\nlog file:\n{}'.format(CONFIGFILE_PATH, logger.get_logfile())
    email('Watermelon completed ok: ', msg=message)
onerror:
    message = "Watermelon completed with errors. Full log file attached"
    attach_str = "-a " + logger.get_logfile() + " --"
    email('Watermelon completed with errors: ', msg=message, attachment = attach_str)


if config.get('fastq_screen'):
    FASTQ_SCREEN_CONFIG = config['fastq_screen']
    FASTQ_SCREEN_ALIGNMENT = expand(ALIGNMENT_DIR + "03-fastq_screen/{screen_type}/{basename}_trimmed_screen.html",
                                screen_type=['multi_species', 'biotype'],
                                basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
                            )
    FASTQ_SCREEN_DELIVERABLES = expand(DELIVERABLES_DIR + "alignment/fastq_screen/{screen_type}/{basename}_trimmed_screen.html",
                                    screen_type=['multi_species', 'biotype'],
                                    basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
                                )
else:
    FASTQ_SCREEN_CONFIG = defaultdict(str)
    FASTQ_SCREEN_ALIGNMENT = []
    FASTQ_SCREEN_DELIVERABLES = []

if config.get('diffex'):
    PHENOTYPE_MANAGER = rnaseq_snakefile_helper.PhenotypeManager(config)
    DIFFEX_MODEL_INFO, DESEQ2_CONTRAST_DICT = rnaseq_snakefile_helper.diffex_model_info(config['diffex'])
    DESeq2_ALL = [
        #deseq2_counts
        DIFFEX_DIR + 'deseq2/counts/count_data.rda',
        expand(DIFFEX_DIR + 'deseq2/counts/deseq2_{name}.txt',
            name=['raw_counts', 'depth_normalized_counts', 'rlog_normalized_counts']),
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
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/BoxPlot_{transformation}.pdf', phenotype = PHENOTYPES, transformation=['raw', 'rlog']),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/SampleHeatmap.pdf', phenotype = PHENOTYPES),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/Heatmap_TopVar.pdf', phenotype = PHENOTYPES),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/Heatmap_TopExp.pdf', phenotype = PHENOTYPES),
        #deseq2_comparison_plots
        rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DIFFEX_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.pdf',
            DESEQ2_CONTRAST_DICT),
        rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DIFFEX_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.png',
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

    DESeq2_DELIVERABLES = [
        #deseq2 deliverables
        expand(DELIVERABLES_DIR + 'counts/deseq2_{name}.txt',
            name=['raw_counts', 'depth_normalized_counts', 'rlog_normalized_counts']),
        rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DELIVERABLES_DIR + 'deseq2/gene_lists/{model_name}/{contrast}.annot.txt',
            DESEQ2_CONTRAST_DICT),
        rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DELIVERABLES_DIR + "deseq2/gene_lists/{model_name}/{contrast}.xlsx",
            DESEQ2_CONTRAST_DICT),
        expand(DELIVERABLES_DIR + 'deseq2/plots/by_phenotype/{phenotype}/PCAplot_{dim}_top{ngenes}.pdf',
                phenotype = PHENOTYPES,
                dim = ['12','23'],
                ngenes = ['100','500']),
        expand(DELIVERABLES_DIR + 'deseq2/plots/by_phenotype/{phenotype}/ScreePlot_top{ngenes}.pdf',
                phenotype = PHENOTYPES,
                ngenes = ['100','500']),
        expand(DELIVERABLES_DIR + 'deseq2/plots/by_phenotype/{phenotype}/{plotType}.pdf',
            phenotype = PHENOTYPES,
            plotType = ['BoxPlot_raw', 'BoxPlot_rlog', 'SampleHeatmap', 'Heatmap_TopVar', 'Heatmap_TopExp']),
        rnaseq_snakefile_helper.expand_model_contrast_filenames(\
                DELIVERABLES_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.pdf',
                DESEQ2_CONTRAST_DICT),
        DELIVERABLES_DIR + "deseq2/summary/deseq2_summary.txt",
        DELIVERABLES_DIR + "deseq2/summary/deseq2_summary.xlsx"
    ]

else:
    DESeq2_ALL = []
    DESeq2_DELIVERABLES = []


# If  running from counts onward, don't populate these
if config.get('count_matrix'):
    ALIGN_DELIVERABLES = []
    RSEM_ALL = []
else:
    ALIGN_DELIVERABLES = [
        #align deliverables
        expand(DELIVERABLES_DIR + "alignment/trimmed_reads/{basename}_trimmed.fastq.gz",
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
        ),
        expand(DELIVERABLES_DIR + "alignment/aligned_bams/{sample}.genome.bam",
            sample=config["samples"]),
        expand(DELIVERABLES_DIR + "alignment/trimmed_reads_fastqc/{basename}_trimmed_fastqc.html",
             basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])),
        expand(DELIVERABLES_DIR + "counts/gene_{type}.annot.txt",
            type=['FPKM', 'TPM', 'expected_count']),
        DELIVERABLES_DIR + "alignment/alignment_qc.html"
    ]
    RSEM_ALL = [
        expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.{result_type}',
            sample=config[SAMPLES_KEY],
            result_type=['genes.results', 'isoforms.results', 'genome.bam', 'transcript.bam']
        ),
        ALIGNMENT_DIR + "07-qc/alignment_qc.html",
        expand(ALIGNMENT_DIR + '06-annotate_combined_counts/{feature}_{metric}.annot.txt',
            feature=['gene','isoform'],
            metric=['expected_count', 'FPKM', 'TPM']
        )
    ]

RUN_INFO_DELIVERABLES = [
    #run info deliverables
    DELIVERABLES_DIR + "run_info/env_software_versions.yaml",
    DELIVERABLES_DIR + "run_info/" + os.path.basename(CONFIGFILE_PATH),
    DELIVERABLES_DIR + "run_info/" + os.path.basename(config['samplesheet'])
]

REPORT_ALL = [
    REPORT_DIR + 'report_draft.md',
    REPORT_DIR + 'report_draft.html'
]


# Various target lists may or may not be populated, see if/else defs above
ALL = RSEM_ALL + ALIGN_DELIVERABLES + RUN_INFO_DELIVERABLES + FASTQ_SCREEN_ALIGNMENT + FASTQ_SCREEN_DELIVERABLES + DESeq2_ALL + DESeq2_DELIVERABLES + REPORT_ALL


include: 'rules/align_concat_reads.smk'
include: 'rules/align_standardize_gz.smk'

if 'trimming_options' in config:
    if config['trimming_options'].get('paired_end_mode'):
        include: 'rules/align_cutadapt_PE.smk'
    else:
        include: 'rules/align_cutadapt_SE.smk'
else:
    include: 'rules/align_pseudotrim.smk'

if config.get('count_matrix'):
    include: 'rules/deseq2_counts_from_matrix.smk'
    include: 'rules/report_from_counts.smk'
else:
    include: 'rules/align_fastq_screen_biotype.smk'
    include: 'rules/align_fastq_screen_multi_species.smk'
    include: 'rules/align_fastqc_trimmed_reads.smk'
    include: 'rules/align_fastqc_align.smk'
    include: 'rules/align_multiqc.smk'
    include: 'rules/align_deliverables_alignment.smk'
    include: 'rules/align_deliverables_fastq_screen.smk'
    include: 'rules/align_rsem_star.smk'
    include: 'rules/align_rsem_star_genome_generate.smk'
    include: 'rules/align_combine_counts_to_matrices.smk'
    include: 'rules/align_annotate_combined_counts.smk'

if config.get('diffex'):
    if not config.get('count_matrix'):
        include: 'rules/deseq2_counts_from_tximport_rsem.smk'
        # TODO: Find a way to manage optional inputs within the reporting rule. For now using different reporting rules
        include: 'rules/report_align_diffex.smk'
    include: 'rules/deseq2_init.smk'
    include: 'rules/deseq2_contrasts.smk'
    include: 'rules/deseq2_plots_by_phenotype.smk'
    include: 'rules/deseq2_volcano_plots.smk'
    include: 'rules/deseq2_annotation.smk'
    include: 'rules/deseq2_excel.smk'
    include: 'rules/deseq2_summary.smk'
    include: 'rules/deliverables_deseq2.smk'

else:
    include: 'rules/report_align_only.smk'

include: 'rules/deliverables_run_info.smk'
include: 'rules/report_finalize.smk'


rule all:
    input:
        *ALL
