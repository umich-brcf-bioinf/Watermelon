rule deseq2_plots_by_phenotype:
    input:
        rda = DIFFEX_DIR + 'counts/count_data.rda'
    output:
        expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{plot_type}.{extension}',
            plot_type = ['SampleHeatmap', 'Heatmap_TopVar', 'Heatmap_TopExp'],
            extension = ['pdf', 'png']),
        expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{phenotype}/PCAplot_{dim}_top{ngenes}.{extension}',
            phenotype = PHENOTYPES,
            dim = ['12','23'],
            ngenes = ['100','500'],
            extension = ['pdf', 'png']),
        expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{phenotype}/ScreePlot_top{ngenes}.{extension}',
            phenotype = PHENOTYPES,
            ngenes = ['100','500'],
            extension = ['pdf', 'png']),
        expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{phenotype}/BoxPlot_{transformation}.{extension}',
            phenotype = PHENOTYPES,
            transformation=['raw', 'rlog'],
            extension = ['pdf', 'png']),
    log:
        JOB_LOG_DIR + 'deseq2_plots_by_phenotype.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    container: 'docker://umichbfxcore/wat_diffex:0.2.0'
    params:
        project_name = config['report_info']['project_name'],
        sample_phenotypes = PHENOTYPE_SAMPLE_LIST,
        diffex_dir = DIFFEX_DIR,
        snakemake_rdata = DIFFEX_DIR + 'plots_labeled_by_pheno/.deseq2_plots_by_phenotype_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_plots_by_phenotype.R'
