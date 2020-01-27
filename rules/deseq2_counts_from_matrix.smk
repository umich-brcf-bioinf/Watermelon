'''This rule imports count data from a count matrix,
    outputs the raw counts and normalized counts as calculated by DESeq2
'''
rule deseq2_counts_from_matrix:
    input:
        config['count_matrix']
    output:
        count_data_rda = DIFFEX_DIR + 'deseq2/counts/count_data.rda',
        raw = DIFFEX_DIR + 'deseq2/counts/deseq2_raw_counts.txt',
        norm = DIFFEX_DIR + 'deseq2/counts/deseq2_depth_normalized_counts.txt',
        rlog = DIFFEX_DIR + 'deseq2/counts/deseq2_rlog_normalized_counts.txt'
    log:
        DIFFEX_DIR + 'deseq2/counts/.log/DESeq2_counts.log'
    conda:
        'envs/diffex.yaml'
    params:
        snakemake_rdata = DIFFEX_DIR + 'deseq2/counts/.deseq2_counts_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_counts.R'
