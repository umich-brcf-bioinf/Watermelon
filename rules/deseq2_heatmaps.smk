rule deseq2_heatmaps:
    input:
        rda = DIFFEX_DIR + 'counts/count_data.rda'
    output:
        expand(DIFFEX_DIR + 'plots_labeled_by_pheno/SampleHeatmap.{extension}', extension = ['pdf', 'png']),
        expand(DIFFEX_DIR + 'plots_labeled_by_pheno/Heatmap_TopVar.{extension}', extension = ['pdf', 'png']),
        expand(DIFFEX_DIR + 'plots_labeled_by_pheno/Heatmap_TopExp.{extension}', extension = ['pdf', 'png'])
    log:
        JOB_LOG_DIR + 'deseq2_heatmaps.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    container: 'docker://umichbfxcore/wat_diffex:0.2.0'
    params:
        project_name = config['report_info']['project_name'],
        sample_phenotypes = PHENOTYPE_MANAGER.phenotype_sample_list,
        diffex_dir = DIFFEX_DIR,
        snakemake_rdata = DIFFEX_DIR + 'plots_labeled_by_pheno/.deseq2_heatmaps_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_heatmaps.R'
