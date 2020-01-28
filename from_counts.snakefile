
import pandas as pd

from scripts import rnaseq_snakefile_helper, __version__


WAT_VER = __version__

_DIRS = config.get("dirs", {})
#INPUT_DIR = os.path.join(_DIRS.get("input", "inputs"), "")
#ALIGNMENT_DIR = os.path.join(_DIRS.get("alignment_output", "alignment_results"), "")
DIFFEX_DIR = os.path.join(_DIRS.get("diffex_output", "diffex_results"), "")
DELIVERABLES_DIR = os.path.join(_DIRS.get("deliverables_output", "deliverables"), "")
WATERMELON_SCRIPTS_DIR = os.path.join(os.environ.get('WATERMELON_SCRIPTS_DIR', srcdir('scripts')), '')
WORKFLOW_BASEDIR = workflow.basedir

CONFIGFILE_PATH = workflow.overwrite_configfile

samplesheet = pd.read_csv(config["samplesheet"], comment='#', dtype='object').set_index("sample", drop=True)

PHENOTYPES = (list(samplesheet.columns))


DESEQ2_CONTRAST_DICT = rnaseq_snakefile_helper.diffex_contrasts(config['diffex'])
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

DELIVERABLES = [
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
        phenotype = PHENOTYPES, plotType = ['BoxPlot_raw', 'BoxPlot_rlog', 'SampleHeatmap', 'Heatmap_TopVar', 'Heatmap_TopExp']),
    rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DELIVERABLES_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.pdf',
            DESEQ2_CONTRAST_DICT),
    DELIVERABLES_DIR + "deseq2/summary/deseq2_summary.txt",
    DELIVERABLES_DIR + "deseq2/summary/deseq2_summary.xlsx",
    #run info deliverables
    DELIVERABLES_DIR + "run_info/env_software_versions.yaml",
    DELIVERABLES_DIR + "run_info/" + os.path.basename(CONFIGFILE_PATH),
    DELIVERABLES_DIR + "run_info/" + os.path.basename(config['samplesheet'])
]


ALL = DESeq2_ALL + DELIVERABLES

include: 'rules/deseq2_counts_from_matrix.smk'
include: 'rules/deseq2_init.smk'
#include: 'rules/deseq2_init_from_countsTable.smk'
include: 'rules/deseq2_contrasts.smk'
include: 'rules/deseq2_plots_by_phenotype.smk'
include: 'rules/deseq2_volcano_plots.smk'
include: 'rules/deseq2_annotation.smk'
include: 'rules/deseq2_excel.smk'
include: 'rules/deseq2_summary.smk'

include: 'rules/deliverables_deseq2.smk'
include: 'rules/deliverables_run_info.smk'

rule all:
    input: DESeq2_ALL + DELIVERABLES
