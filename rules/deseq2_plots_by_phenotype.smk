rule deseq2_plots_by_phenotype:
    input:
        rda = DIFFEX_DIR + 'deseq2/counts/count_data.rda'
    output:
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/PCAplot_{dim}_top{ngenes}.{extension}',
            phenotype = PHENOTYPES,
            dim = ['12','23'],
            ngenes = ['100','500'],
            extension = ['pdf', 'png']),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/ScreePlot_top{ngenes}.{extension}',
            phenotype = PHENOTYPES,
            ngenes = ['100','500'],
            extension = ['pdf', 'png']),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/BoxPlot_{transformation}.{extension}', phenotype = PHENOTYPES, transformation=['raw', 'rlog'], extension = ['pdf', 'png']),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/SampleHeatmap.{extension}', phenotype = PHENOTYPES, extension = ['pdf', 'png']),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/Heatmap_TopVar.{extension}', phenotype = PHENOTYPES, extension = ['pdf', 'png']),
        expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/Heatmap_TopExp.{extension}', phenotype = PHENOTYPES, extension = ['pdf', 'png'])
    log:
        DIFFEX_DIR + 'deseq2/plots/by_phenotype/.log/deseq2_plots_by_phenotype.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    singularity: 'docker://umichbfxcore/wat_diffex'
    params:
        phenotypes = PHENOTYPES,
        diffex_dir = DIFFEX_DIR,
        snakemake_rdata = DIFFEX_DIR + 'deseq2/plots/by_phenotype/.deseq2_plots_by_phenotype_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_plots_by_phenotype.R'
