#This rule initializes the DESeq2 object for each given model.
rule deseq2_init:
    input:
        data_import = DIFFEX_DIR + 'deseq2/counts/txi_rsem_genes.rda',
    output:
        #expand(DIFFEX_DIR + '{{model_name}}/DESeq2/{contrast}_gene.results', contrast=rnaseq_snakefile_helper.get_DESeq2_model_contrasts(config['diffex'], wildcards.model_name))
        rda = DIFFEX_DIR + 'deseq2/deseq2_init_{model_name}.rda'
    threads: 8
    log:
        DIFFEX_DIR + 'deseq2/.log/deseq2_init_{model_name}.log'
    conda:
        'envs/diffex.yaml'
    params:
        snakemake_rdata = DIFFEX_DIR + 'deseq2/deseq2_init_{model_name}_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_init.R'