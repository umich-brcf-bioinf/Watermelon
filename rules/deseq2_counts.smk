'''This rule imports DESeq2 data from RSEM,
    outputs the raw counts and normalized counts as calculated by DESeq2
'''
rule deseq2_counts:
    input:
        expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genes.results', sample=config[SAMPLES_KEY]),
        expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.isoforms.results', sample=config[SAMPLES_KEY])
    output:
        txi = DIFFEX_DIR + 'counts/txi_rsem_genes.rda',
        count_tables_rda = DIFFEX_DIR + 'counts/count_data.rda',
        raw = DIFFEX_DIR + 'counts/raw_counts.txt',
        norm = DIFFEX_DIR + 'counts/depth_normalized_counts.txt',
        rlog = DIFFEX_DIR + 'counts/rlog_normalized_counts.txt'
    log:
        DIFFEX_DIR + 'counts/.log/DESeq2_counts.log'
    conda:
        'envs/diffex.yaml'
    params:
        rsem_dir = ALIGNMENT_DIR + '04-rsem_star_align',
        snakemake_rdata = DIFFEX_DIR + 'counts/deseq2_counts_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_counts.R'
