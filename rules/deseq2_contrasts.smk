rule deseq2_contrasts:
    input:
        rda = DIFFEX_DIR + 'deseq2_init_{model_name}.rda',
        gene_info = config['references']['annotation_tsv'],
        glossary = WATERMELON_SCRIPTS_DIR + 'deseq2_glossary.txt',
    output:
        gene_list = DIFFEX_DIR + '{model_name}/gene_lists/{contrast}.txt',
        annot_results = DIFFEX_DIR + "{model_name}/annotated/{contrast}.annot.txt",
        annot_results_xlsx = DIFFEX_DIR + "{model_name}/excel/{contrast}.xlsx",
        volcano_plot_pdf = DIFFEX_DIR + '{model_name}/volcano_plots/VolcanoPlot_{contrast}.pdf',
        volcano_plot_png = DIFFEX_DIR + '{model_name}/volcano_plots/VolcanoPlot_{contrast}.png'
        #rda = DIFFEX_DIR + '{model_name}/DESeq2/{contrast}_data.rda'
    log:
        JOB_LOG_DIR + 'deseq2_contrast_{model_name}_{contrast}.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    singularity: 'docker://umichbfxcore/wat_diffex:0.2.0'
    resources: cpus=8, mem_mb=8000
    params:
        project_name = config['report_info']['project_name'],
        snakemake_rdata = DIFFEX_DIR + '{model_name}/gene_lists/.{contrast}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_contrasts.R'
