rule deseq2_contrasts:
    input:
        rda = DIFFEX_DIR + 'deseq2_init_{model_name}.rda'
    output:
        gene_list = DIFFEX_DIR + '{model_name}/gene_lists/{contrast}.txt',
        #rda = DIFFEX_DIR + '{model_name}/DESeq2/{contrast}_data.rda'
    log:
        JOB_LOG_DIR + 'deseq2_contrast_{model_name}_{contrast}.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    singularity: 'docker://umichbfxcore/wat_diffex:0.1.1'
    resources: cpus=8
    params:
        project_name = config['report_info']['project_name'],
        snakemake_rdata = DIFFEX_DIR + '{model_name}/gene_lists/.{contrast}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_contrasts.R'
