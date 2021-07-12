#This rule creates a volcano plot for each comparison
rule deseq2_volcano_plots:
    input:
        gene_list = DIFFEX_DIR + '{model_name}/annotated/{contrast}.annot.txt'
    output:
        volcano_plot_pdf = DIFFEX_DIR + '{model_name}/volcano_plots/VolcanoPlot_{contrast}.pdf',
        volcano_plot_png = DIFFEX_DIR + '{model_name}/volcano_plots/VolcanoPlot_{contrast}.png'
    log:
        JOB_LOG_DIR + 'deseq2_volcano_plots_{model_name}_{contrast}.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    singularity: 'docker://umichbfxcore/wat_diffex:0.1.1'
    params:
        project_name = config['report_info']['project_name'],
        method = "deseq2",
        snakemake_rdata = DIFFEX_DIR + '{model_name}/volcano_plots/.volcano_plot_{contrast}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/diffex_volcano_plot.R'
