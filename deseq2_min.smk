
import logging
import os
import pandas as pd
import yaml

from scripts import rnaseq_snakefile_helper as helper
#from scripts import config_validator

# Set up all of the directories
# Dirs relative to pipeline
WORKFLOW_BASEDIR = workflow.basedir
WATERMELON_CONFIG_DIR = os.path.join(WORKFLOW_BASEDIR, 'config', '')
WATERMELON_SCRIPTS_DIR = os.path.join(WORKFLOW_BASEDIR, 'scripts', '')
# Output directories from config
_DIRS = config.get("dirs", {})
DIFFEX_DIR = os.path.join(_DIRS.get("diffex_results", "diffex_results"), "")
DELIVERABLES_DIR = os.path.join(_DIRS.get("deliverables", "deliverables"), "")
REPORT_DIR = os.path.join(_DIRS.get("report", "report"), "")
# Logging directories
JOB_LOG_DIR = os.path.join(os.getcwd(), "job_logs", "")
CLUSTER_LOG_DIR = os.path.join(os.getcwd(), "cluster_logs")

# Load in samplesheet
# Drop either input_dir or input_glob if they exist - transition period dir -> glob
SAMPLESHEET = pd.read_csv(config["samplesheet"], comment='#', dtype='object') \
    .drop("input_dir", axis=1, errors="ignore") \
    .drop("input_glob", axis=1, errors="ignore") \
    .set_index("sample", drop=True)

if 'phenotypes' in config:
    PHENOTYPES = config['phenotypes']
else:
    PHENOTYPES = (list(SAMPLESHEET.columns))

PHENOTYPE_SAMPLE_LIST = helper.phenotype_sample_list(SAMPLESHEET)

#Get config file (first need to grab index from argv)
if '--configfile' in list(sys.argv):
    cfg_idx = sys.argv.index('--configfile') + 1
elif '--configfiles' in list(sys.argv): # Cluster environment uses configfiles (plural) under the hood
    cfg_idx = sys.argv.index('--configfiles') + 1
else:
    msg = '--configfile not specified in {}'.format(sys.argv)
    raise ValueError(msg)
# Now with the index, get the value
CONFIGFILE_PATH = sys.argv[cfg_idx]

#function to call config_validator
def validate_config(config_fp, schema_fp):
    config_validation_exit_code = config_validator.main(config_fp, schema_fp)
    if config_validation_exit_code:
        exit(config_validation_exit_code)

# TODO: Perform config validation if dryrun?
if set(['-n', '--dryrun']).intersection(set(sys.argv)) and not 'skip_validation' in config:
    foo = "bar" #TWS DEBUG
    #validate_config(CONFIGFILE_PATH, CONFIG_SCHEMA_PATH)


onstart:
    helper.email(email_config=config.get('email', None), subject_prefix='Watermelon started: ')
onsuccess:
    message = 'config file:\n{}\nlog file:\n{}'.format(CONFIGFILE_PATH, logger.get_logfile())
    helper.email(email_config=config.get('email', None), subject_prefix='Watermelon completed ok: ', msg=message)
onerror:
    message = "Watermelon completed with errors. Full log file attached"
    attach_str = "-a " + logger.get_logfile() + " --"
    helper.email(email_config=config.get('email', None), subject_prefix='Watermelon completed with errors: ', msg=message, attachment = attach_str)


DESeq2_ALL = [
    #deseq2_counts
    DIFFEX_DIR + 'counts/count_data.rda',
    expand(DIFFEX_DIR + 'counts/deseq2_{name}.txt',
        name=['raw_counts', 'depth_normalized_counts', 'rlog_normalized_counts']),
    #deseq2_plots_by_phenotype
    expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{phenotype}/PCAplot_{dim}_top{ngenes}.pdf',
        phenotype = PHENOTYPES,
        dim = ['12','23'],
        ngenes = ['100','500']),
    expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{phenotype}/ScreePlot_top{ngenes}.pdf',
        phenotype = PHENOTYPES,
        ngenes = ['100','500']),
    expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{phenotype}/BoxPlot_{transformation}.pdf', phenotype = PHENOTYPES, transformation=['raw', 'rlog']),
    expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{plot_type}.pdf',
        plot_type=["SampleHeatmap", "Heatmap_TopVar", "Heatmap_TopExp"]),
    #deseq2_init rda data
    expand(DIFFEX_DIR + 'diffex_{model_name}/.deseq2_init_{model_name}.rda',
        model_name=helper.diffex_models(config['diffex']))
]

RUN_INFO_DELIVERABLES = [
    DELIVERABLES_DIR + "run_info/env_software_versions.yaml",
    DELIVERABLES_DIR + "run_info/" + os.path.basename(config['samplesheet'])
]

rule all:
    input:
        DESeq2_ALL + RUN_INFO_DELIVERABLES


# Includes put last - variables used in rules must be defined above
# version_info.smk gives access to VER_INFO with software version info
# and also to ENV_INFO - a singularity info dict
include: "version_info.smk"

if config.get('count_matrix'):
    include: 'rules/deseq2_counts_from_matrix.smk'

include: 'rules/deseq2_init.smk'
include: 'rules/deseq2_plots_by_phenotype.smk'
include: 'rules/deliverables_run_info.smk'
include: 'rules/report_finalize.smk'