rule deseq2_contrasts:
    input:
        rda = DIFFEX_DIR + '{model_name}/DESeq2/deseq2_init.rda'
    output:
        gene_list = DIFFEX_DIR + '{model_name}/DESeq2/{contrast}_gene.results',
        #rda = DIFFEX_DIR + '{model_name}/DESeq2/{contrast}_data.rda'
    threads: 8
    log:
        DIFFEX_DIR + '{model_name}/DESeq2/.log/{contrast}_deseq2_contrast.log'
    conda:
        'envs/diffex.yaml'
    params:
        snakemake_rdata = DIFFEX_DIR + '{model_name}/DESeq2/{contrast}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_contrasts.R'
