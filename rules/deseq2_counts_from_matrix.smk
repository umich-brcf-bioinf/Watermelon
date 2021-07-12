'''This rule imports count data from a count matrix,
    outputs the raw counts and normalized counts as calculated by DESeq2
'''
rule deseq2_counts_from_matrix:
    input:
        count_matrix = config['count_matrix']
    output:
        count_data_rda = DIFFEX_DIR + 'counts/count_data.rda',
        raw = DIFFEX_DIR + 'counts/deseq2_raw_counts.txt',
        norm = DIFFEX_DIR + 'counts/deseq2_depth_normalized_counts.txt',
        rlog = DIFFEX_DIR + 'counts/deseq2_rlog_normalized_counts.txt'
    log:
        JOB_LOG_DIR + 'deseq2_counts_from_matrix.log'
    conda: 'envs/WAT_diffex/WAT_diffex.yaml'
    singularity: 'docker://umichbfxcore/wat_diffex:0.1.1'
    params:
        project_name = config['report_info']['project_name'],
        snakemake_rdata = DIFFEX_DIR + 'counts/.deseq2_counts_snakemake.rda' #TWS DEBUG
    script:
        '../scripts/deseq2_counts.R'
