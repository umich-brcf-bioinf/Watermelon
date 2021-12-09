
import logging
import os
import pandas as pd
import yaml

from argparse import Namespace

from scripts import rnaseq_snakefile_helper as helper
from scripts import get_run_stats

# Set up all of the directories
# Dirs relative to pipeline
WORKFLOW_BASEDIR = workflow.basedir
WATERMELON_CONFIG_DIR = os.path.join(WORKFLOW_BASEDIR, 'config', '')
WATERMELON_SCRIPTS_DIR = os.path.join(WORKFLOW_BASEDIR, 'scripts', '')
# Output directories from config
_DIRS = config.get("dirs", {})
ALIGNMENT_DIR = os.path.join(_DIRS.get("alignment_results", "alignment_results"), "")
DELIVERABLES_DIR = os.path.join(_DIRS.get("deliverables", "deliverables"), "")
REPORT_DIR = os.path.join(_DIRS.get("report", "report"), "")
# Logging directories
JOB_LOG_DIR = os.path.join(os.getcwd(), "job_logs", "")
CLUSTER_LOG_DIR = os.path.join(os.getcwd(), "cluster_logs")

#Load in samplesheet
SAMPLESHEET = pd.read_csv(config["samplesheet"], comment='#', dtype='object') \
    .set_index("sample", drop=True)

# Allow custom capture regex if specified in config
if config.get('capture_regex'):
    capture_regex = config['capture_regex']
else:
    capture_regex = r'.*_R(\d+)[_L0-9]*\.fastq.*'

FASTQS_TO_CONCAT = helper.fastqs_to_concat(SAMPLESHEET, capture_regex)
SAMPLE_BNAMES = helper.sample_bnames_from_filenames(SAMPLESHEET, capture_regex, '{}_R{}')

# Determine if run in cluster environment
# If so, make cluster_logs folder if it's missing
if set(['-c', '--cluster']).intersection(set(sys.argv)):
    if not os.path.exists(CLUSTER_LOG_DIR):
        msg = "Cluster log folder {} not found. Creating it before running in a cluster environment."
        logger.logger.info(msg.format(CLUSTER_LOG_DIR))
        os.mkdir(CLUSTER_LOG_DIR)
# May also need to parse the profile to see cluster argument
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

onstart:
    helper.email(email_config=config.get('email', None), subject_prefix='Watermelon started: ')
onsuccess:
    message = 'config file:\n{}\nlog file:\n{}'.format(CONFIGFILE_PATH, logger.get_logfile())
    helper.email(email_config=config.get('email', None), subject_prefix='Watermelon completed ok: ', msg=message)
    # Gather run stats
    runstats_args = Namespace(infile=logger.get_logfile(), outfile=None)
    get_run_stats.main(runstats_args) # Do this only when run on cluster?
onerror:
    message = "Watermelon completed with errors. Full log file attached"
    attach_str = "-a " + logger.get_logfile() + " --"
    helper.email(email_config=config.get('email', None), subject_prefix='Watermelon completed with errors: ', msg=message, attachment = attach_str)
    # Gather run stats
    runstats_args = Namespace(infile=logger.get_logfile(), outfile=None)
    get_run_stats.main(runstats_args) # Do this only when run on cluster?


# Defining targets
ALIGN_DELIVERABLES = [
    expand(DELIVERABLES_DIR + "trimmed/trimmed_reads/{basename}_trimmed.fastq.gz",
        basename=helper.gather_basenames(SAMPLE_BNAMES, SAMPLESHEET.index)),
    expand(DELIVERABLES_DIR + "alignment/aligned_bams/{sample}.genome.{ext}",
        sample=SAMPLESHEET.index, ext=["bam", "bai"]),
    expand(DELIVERABLES_DIR + "trimmed/trimmed_reads_fastqc/{basename}_trimmed_fastqc.html",
         basename=helper.gather_basenames(SAMPLE_BNAMES, SAMPLESHEET.index)),
    expand(DELIVERABLES_DIR + "counts/gene_{type}.annot.txt",
        type=['FPKM', 'TPM', 'expected_count']),
    DELIVERABLES_DIR + "alignment/alignment_qc.html"
]
RSEM_ALL = [
    expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.{result_type}',
        sample=SAMPLESHEET.index,
        result_type=['genes.results', 'isoforms.results', 'genome.bam', 'transcript.bam']
    ),
    ALIGNMENT_DIR + "07-qc/alignment_qc.html",
    expand(ALIGNMENT_DIR + '06-annotate_combined_counts/{feature}_{metric}.annot.txt',
        feature=['gene','isoform'],
        metric=['expected_count', 'FPKM', 'TPM']
    )
]

RUN_INFO_DELIVERABLES = [
    DELIVERABLES_DIR + "run_info/env_software_versions.yaml",
    DELIVERABLES_DIR + "run_info/" + os.path.basename(config['samplesheet'])
]

REPORT_ALL = [
    REPORT_DIR + 'report_draft.md',
    REPORT_DIR + 'report_draft.html'
]

if config.get('fastq_screen'):
    FASTQ_SCREEN_ALIGNMENT = expand(ALIGNMENT_DIR + "03-fastq_screen/{screen_type}/{basename}_trimmed_screen.html",
        screen_type=['multi_species', 'biotype'],
        basename=helper.gather_basenames(SAMPLE_BNAMES, SAMPLESHEET.index)),
        # ^ This is useful for the inclusion/exclusion in multiqc rule as necessary
    FASTQ_SCREEN_DELIVERABLES = expand(DELIVERABLES_DIR + "fastq_screen/{screen_type}/{basename}_trimmed_screen.html",
        screen_type=['multi_species', 'biotype'],
        basename=helper.gather_basenames(SAMPLE_BNAMES, SAMPLESHEET.index))
else:
    FASTQ_SCREEN_ALIGNMENT = []
    # ^ This is useful for the inclusion/exclusion in multiqc rule as necessary
    FASTQ_SCREEN_DELIVERABLES = []

if "rseqc" in config:
    RSEQC_RIBO_EST = expand(ALIGNMENT_DIR + "05-rseqc/{sample}_ribo_count.txt", sample=SAMPLESHEET.index)
else:
    RSEQC_RIBO_EST = []
    # ^ This is useful for the inclusion/exclusion in multiqc rule as necessary

# Target list - may or may not be populated, see if/else defs above
ALL = RSEM_ALL + ALIGN_DELIVERABLES + RUN_INFO_DELIVERABLES + FASTQ_SCREEN_DELIVERABLES + REPORT_ALL


rule all:
    input:
        *ALL

# Includes put last - variables used in rules must be defined above
# version_info.smk gives access to VER_INFO with software version info
# and also to ENV_INFO - a singularity info dict
include: "version_info.smk"

# Include rule snakefiles
include: 'rules/align_concat_reads.smk'
include: 'rules/align_standardize_gz.smk'
#
if 'trimming' in config:
    if config['trimming'].get('paired_end_mode'):
        include: 'rules/align_cutadapt_PE.smk'
    else:
        include: 'rules/align_cutadapt_SE.smk'
else:
    include: 'rules/align_pseudotrim.smk'
if 'fastq_screen' in config:
    include: 'rules/align_fastq_screen_biotype.smk'
    include: 'rules/align_fastq_screen_multi_species.smk'
if "rseqc" in config:
    include: 'rules/align_rseqc_estimate_ribo.smk'
include: 'rules/align_fastqc_trimmed_reads.smk'
include: 'rules/align_fastqc_align.smk'
include: 'rules/align_multiqc.smk'
include: 'rules/align_deliverables_alignment.smk'
include: 'rules/align_deliverables_fastq_screen.smk'
include: 'rules/align_rsem_star.smk'
include: 'rules/align_rsem_star_genome_generate.smk'
include: 'rules/align_combine_counts_to_matrices.smk'
include: 'rules/align_annotate_combined_counts.smk'
#
include: 'rules/report_align_qc.smk'
#
include: 'rules/deliverables_run_info.smk'
include: 'rules/report_finalize.smk'
