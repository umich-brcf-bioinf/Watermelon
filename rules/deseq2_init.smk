#This rule initializes the DESeq2 object for each given model.
rule deseq2_init:
    input:
        data_import = DIFFEX_DIR + 'counts/count_data.rda',
    output:
        rda = DIFFEX_DIR + 'diffex_{model_name}/.deseq2_init_{model_name}.rda'
    log:
        JOB_LOG_DIR + 'deseq2_init_{model_name}.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    container: 'docker://umichbfxcore/wat_diffex:0.3.1'
    resources: cpus=8, mem_mb=8000
    params:
        project_name = config['report_info']['project_name'],
        snakemake_rdata = DIFFEX_DIR + 'diffex_{model_name}/.deseq2_init_{model_name}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_init.R'
