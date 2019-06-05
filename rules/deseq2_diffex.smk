rule deseq2_diffex:
    input:
        expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genes.results', sample=config[SAMPLES_KEY]),
        expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.isoforms.results', sample=config[SAMPLES_KEY])
    output:
        iso_results = DIFFEX_DIR + '{model_name}/DESeq2/{contrast}_isoform.results',
        gene_results = DIFFEX_DIR + '{model_name}/DESeq2/{contrast}_gene.results',
        rda = DIFFEX_DIR + '{model_name}/DESeq2/{contrast}.rda'
    threads: 8
    log:
        DIFFEX_DIR + '{model_name}/DESeq2/.log/{contrast}_deseq2_DESeq2Diffex.log'
    conda:
        'envs/diffex.yaml'
    params:
        rsem_dir = ALIGNMENT_DIR + '04-rsem_star_align',
        snakemake_rdata = DIFFEX_DIR + '{model_name}/DESeq2/{contrast}_snakemake.rda' #TWS DEBUG
    shell:
        '../scripts/deseq2_diffex.R'
