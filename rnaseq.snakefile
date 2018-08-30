from __future__ import print_function, absolute_import, division

from collections import defaultdict, OrderedDict
from itertools import combinations, repeat
from shutil import copyfile
import csv
import os
import subprocess
import yaml

import scripts
import scripts.rnaseq_snakefile_helper as rnaseq_snakefile_helper
import scripts.deseq2_helper as deseq2_helper

WAT_VER = scripts.__version__
WATERMELON_CONFIG_DIR = os.path.join(os.environ.get('WATERMELON_CONFIG_DIR', srcdir('config')), '')
WATERMELON_SCRIPTS_DIR = os.path.join(os.environ.get('WATERMELON_SCRIPTS_DIR', srcdir('scripts')), '')

rnaseq_snakefile_helper.transform_config(config)
_DIRS = config.get("dirs", {})
INPUT_DIR = os.path.join(_DIRS.get("input", "inputs"), "")
ALIGNMENT_DIR = os.path.join(_DIRS.get("alignment_output", "alignment_results"), "")
DIFFEX_DIR = os.path.join(_DIRS.get("diffex_output", "diffex_results"), "")
DELIVERABLES_DIR = os.path.join(_DIRS.get("deliverables_output", "deliverables"), "")
DESEQ2_DIR = os.path.join(DIFFEX_DIR, "deseq2", "")
TUXEDO_DIR = os.path.join(DIFFEX_DIR, "tuxedo", "")
CONFIG_CHECKSUMS_DIR = os.path.join(".config_checksums", "")

COMPARISON_INFIX = '_v_'
DELIMITER = '^'
SAMPLES_KEY = 'samples'
PHENOTYPES_KEY = 'phenotypes'
COMPARISONS_KEY ='comparisons'

phenotypeManager = rnaseq_snakefile_helper.PhenotypeManager(config,
                                                            DELIMITER,
                                                            COMPARISON_INFIX)

ALL_PHENOTYPE_NAMES = phenotypeManager.phenotypes_comparisons_all_tuple.phenotypes
ALL_COMPARISON_GROUPS = phenotypeManager.phenotypes_comparisons_all_tuple.comparisons

REPLICATE_PHENOTYPE_NAMES = phenotypeManager.phenotypes_comparisons_replicates_tuple.phenotypes
REPLICATE_COMPARISON_GROUPS = phenotypeManager.phenotypes_comparisons_replicates_tuple.comparisons


rnaseq_snakefile_helper.init_references(config["references"])
rnaseq_snakefile_helper.checksum_reset_all(CONFIG_CHECKSUMS_DIR,
                                           config=config,
                                           phenotype_comparisons=config['comparisons'],
                                           phenotype_samples=phenotypeManager.phenotype_sample_list)

SAMPLE_READS = rnaseq_snakefile_helper.flattened_sample_reads(INPUT_DIR, config[SAMPLES_KEY])

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

DESEQ2_ALL = []
if REPLICATE_PHENOTYPE_NAMES:
    DESEQ2_ALL = [
        DESEQ2_DIR + "01-metadata_contrasts/sample_metadata.txt",
        DESEQ2_DIR + "01-metadata_contrasts/contrasts.txt",
        expand(DESEQ2_DIR + "02-deseq2_diffex/gene_lists/{phenotype_name}/{comparison}.txt",
               zip,
               phenotype_name=REPLICATE_PHENOTYPE_NAMES,
               comparison=REPLICATE_COMPARISON_GROUPS),
        expand(DESEQ2_DIR + "03-annotation/{phenotype_name}/{comparison}.annot.txt",
               zip,
               phenotype_name=REPLICATE_PHENOTYPE_NAMES,
               comparison=REPLICATE_COMPARISON_GROUPS),
        expand(DESEQ2_DIR + "05-excel/{phenotype_name}/{comparison}.xlsx",
                zip,
                phenotype_name=REPLICATE_PHENOTYPE_NAMES,
                comparison=REPLICATE_COMPARISON_GROUPS),
        DESEQ2_DIR + "06-summary/deseq2_summary.txt",
        DESEQ2_DIR + "06-summary/deseq2_summary.xlsx",
        DELIVERABLES_DIR + "deseq2/gene_lists/deseq2_summary.txt",
        ]

OPTIONAL_ALL = DESEQ2_ALL + FASTQ_SCREEN_ALIGNMENT + FASTQ_SCREEN_DELIVERABLES

ALL = [OPTIONAL_ALL]

include: 'rules/align_concat_reads.smk'
include: 'rules/align_cutadapt_SE.smk'
include: 'rules/align_cutadapt_PE.smk'
include: 'rules/align_fastq_screen_biotype.smk'
include: 'rules/align_fastq_screen_multi_species.smk'
include: 'rules/align_fastqc_trimmed_reads.smk'
include: 'rules/align_insert_size_PE.smk'
include: 'rules/align_create_transcriptome_index.smk'
include: 'rules/align_build_tophat_sample_options.smk'
include: 'rules/align_tophat.smk'
include: 'rules/align_fastqc_tophat_align.smk'
include: 'rules/align_stringtie.smk'
include: 'rules/align_stringtie_prepDE.smk'
include: 'rules/align_qc.smk'
include: 'rules/align_deliverables_alignment.smk'
include: 'rules/align_deliverables_fastq_screen.smk'

include: 'rules/tuxedo_cuffdiff.smk'
include: 'rules/tuxedo_flip.smk'
include: 'rules/tuxedo_flag.smk'
include: 'rules/tuxedo_annotate.smk'
include: 'rules/tuxedo_group_replicates.smk'
include: 'rules/tuxedo_cummerbund.smk'
include: 'rules/tuxedo_split.smk'
include: 'rules/tuxedo_last_split.smk'
include: 'rules/tuxedo_run_info.smk'
include: 'rules/tuxedo_excel.smk'
include: 'rules/tuxedo_summary.smk'

include: 'rules/deseq2_metadata_contrasts.smk'
include: 'rules/deseq2_diffex.smk'
include: 'rules/deseq2_annotation.smk'
include: 'rules/deseq2_run_info.smk'
include: 'rules/deseq2_excel.smk'
include: 'rules/deseq2_summary.smk'

include: 'rules/deliverables_tuxedo.smk'
include: 'rules/deliverables_deseq2.smk'
include: 'rules/deliverables_combined_summary.smk'

rule all:
    input:
        *ALL
