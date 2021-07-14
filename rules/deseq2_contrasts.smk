rule deseq2_contrasts:
    input:
        rda = DIFFEX_DIR + 'deseq2_init_{model_name}.rda',
        gene_info = config['references']['annotation_tsv'],
        glossary = WATERMELON_SCRIPTS_DIR + 'deseq2_glossary.txt',
    output:
        gene_list = DIFFEX_DIR + 'diffex_{model_name}/{contrast}.txt',
        annot_results = DIFFEX_DIR + "diffex_{model_name}/{contrast}.annot.txt",
        annot_results_xlsx = DIFFEX_DIR + "diffex_{model_name}/{contrast}.annot.xlsx",
        volcano_plot_pdf = DIFFEX_DIR + 'volcano_plots_{model_name}/VolcanoPlot_{contrast}.pdf',
        volcano_plot_png = DIFFEX_DIR + 'volcano_plots_{model_name}/VolcanoPlot_{contrast}.png'
    log:
        JOB_LOG_DIR + 'deseq2_contrast_{model_name}_{contrast}.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    singularity: 'docker://umichbfxcore/wat_diffex:0.2.0'
    resources: cpus=8, mem_mb=8000
    params:
        project_name = config['report_info']['project_name'],
        snakemake_rdata = DIFFEX_DIR + 'diffex_{model_name}/.{contrast}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_contrasts.R'
