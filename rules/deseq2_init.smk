#This rule initializes the DESeq2 object for each given model.
rule deseq2_init:
    input:
        data_import = DIFFEX_DIR + 'deseq2/counts/count_data.rda',
    output:
        rda = DIFFEX_DIR + 'deseq2/deseq2_init_{model_name}.rda'
    log:
        JOB_LOG_DIR + 'deseq2_init_{model_name}.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    singularity: 'docker://umichbfxcore/wat_diffex:0.1.1'
    resources: cpus=8
    params:
        snakemake_rdata = DIFFEX_DIR + 'deseq2/.deseq2_init_{model_name}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_init.R'
