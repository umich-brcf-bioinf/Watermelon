#This rule initializes the DESeq2 object for each given model.
rule deseq2_init:
    input:
        data_import = DIFFEX_DIR + 'deseq2/counts/count_data.rda',
    output:
        rda = DIFFEX_DIR + 'deseq2/deseq2_init_{model_name}.rda'
    threads: 8
    log:
        DIFFEX_DIR + 'deseq2/.log/deseq2_init_{model_name}.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    singularity: 'docker://umichbfxcore/wat_diffex:0.1.1'
    params:
        snakemake_rdata = DIFFEX_DIR + 'deseq2/.deseq2_init_{model_name}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_init.R'
