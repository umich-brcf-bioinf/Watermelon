rule deseq2_plots_by_phenotype:
    input:
        rda = DIFFEX_DIR + 'deseq2/counts/count_data.rda'
    output:
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/PCAplot_{dim}_top{ngenes}.pdf',
            phenotype = PHENOTYPES,
            dim = ['12','23'],
            ngenes = ['100','500']),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/MDSplot_{dim}_top{ngenes}.pdf',
            phenotype = PHENOTYPES,
            dim = ['12','23'],
            ngenes = ['100','500']),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/ScreePlot_top{ngenes}.pdf',
            phenotype = PHENOTYPES,
            ngenes = ['100','500']),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/BoxPlot.pdf', phenotype = PHENOTYPES),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/SampleHeatmap.pdf', phenotype = PHENOTYPES),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/Heatmap_TopVar.pdf', phenotype = PHENOTYPES),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/Heatmap_TopExp.pdf', phenotype = PHENOTYPES)
    log:
        DIFFEX_DIR + 'deseq2/plots/by_phenotype/.log/deseq2_plots_by_phenotype.log'
    conda:
        'envs/diffex.yaml'
    params:
        phenotypes = PHENOTYPES,
        diffex_dir = DIFFEX_DIR,
        snakemake_rdata = DIFFEX_DIR + 'deseq2/plots/by_phenotype/.deseq2_plots_by_phenotype_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_plots_by_phenotype.R'
