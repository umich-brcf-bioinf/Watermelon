rule deseq2_contrasts:
    input:
        rda = DIFFEX_DIR + 'deseq2/deseq2_init_{factor_name}.rda'
    output:
        gene_list = DIFFEX_DIR + 'deseq2/gene_lists/{factor_name}/{contrast}.txt',
        #rda = DIFFEX_DIR + '{factor_name}/DESeq2/{contrast}_data.rda'
    threads: 8
    log:
        DIFFEX_DIR + 'deseq2/gene_lists/.log/{factor_name}_{contrast}_deseq2_contrast.log'
    conda:
        'envs/diffex.yaml'
    params:
        snakemake_rdata = DIFFEX_DIR + 'deseq2/gene_lists/{factor_name}/{contrast}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_contrasts.R'
