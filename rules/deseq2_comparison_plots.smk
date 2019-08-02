#This rule creates a volcano plot for each comparison
rule deseq2_comparison_plots:
    input:
        gene_list = DIFFEX_DIR + 'deseq2/gene_lists/{factor_name}/{contrast}.txt'
    output:
        volcano_plot = DIFFEX_DIR + 'deseq2/plots/comparison_plots/{factor_name}/VolcanoPlot_{contrast}.pdf'
    log:
        DIFFEX_DIR + 'deseq2/plots/comparison_plots/{factor_name}/.log/{contrast}_deseq2_comparison_plots.log'
    conda:
        'envs/diffex.yaml'
    params:
        method = "deseq2",
        snakemake_rdata = DIFFEX_DIR + 'deseq2/plots/comparison_plots/{factor_name}/.comparison_plot_{contrast}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/diffex_volcano_plot.R'
