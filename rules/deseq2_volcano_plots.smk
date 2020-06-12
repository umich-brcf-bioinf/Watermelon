#This rule creates a volcano plot for each comparison
rule deseq2_volcano_plots:
    input:
        gene_list = DIFFEX_DIR + 'deseq2/annotated/{model_name}/{contrast}.annot.txt'
    output:
        volcano_plot_pdf = DIFFEX_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.pdf',
        volcano_plot_png = DIFFEX_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.png'
    log:
        DIFFEX_DIR + 'deseq2/plots/comparison_plots/{model_name}/.log/{contrast}_deseq2_comparison_plots.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    singularity: 'docker://umichbfxcore/wat_diffex:0.1.1'
    params:
        method = "deseq2",
        snakemake_rdata = DIFFEX_DIR + 'deseq2/plots/comparison_plots/{model_name}/.comparison_plot_{contrast}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/diffex_volcano_plot.R'
