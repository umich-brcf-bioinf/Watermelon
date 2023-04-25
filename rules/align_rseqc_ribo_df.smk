rule align_rseqc_ribo_df:
    input:
        expand(ALIGNMENT_DIR + "05-rseqc/{sample}_ribo_pct.txt", sample=SAMPLESHEET.index)
    output:
        ALIGNMENT_DIR + "05-rseqc/ribo_pct_summary.txt"
    params:
        project_name = config['report_info']['project_name']
    log:
        JOB_LOG_DIR + "align_rseqc_ribo_df.log"
    run:
        sample_rows = []
        for fn in input:
            sample = pd.read_table(fn, index_col='Samplename')
            sample_rows.append(sample)
        pct_df = pd.concat(sample_rows)
        pct_df.to_csv(output[0], sep='\t')
