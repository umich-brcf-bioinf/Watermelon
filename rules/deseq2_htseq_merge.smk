rule deseq2_htseq_merge:
    input:
        sample_checksum = CONFIG_CHECKSUMS_DIR + "config-samples.watermelon.md5",
        sample_count_files = expand(DESEQ2_DIR + "01-htseq/{sample}_counts.txt",
                                    sample=config[SAMPLES_KEY])
    output:
        counts_filename = DESEQ2_DIR + "01-htseq/htseq_merged.txt",
        stats_filename = DESEQ2_DIR + "01-htseq/htseq_stats.txt",
    params:
        input_dir = DESEQ2_DIR + "01-htseq/",
        counts_filename = DESEQ2_DIR + "01-htseq/htseq_merged.txt",
        stats_filename = DESEQ2_DIR + "01-htseq/htseq_stats.txt",
    log:
        DESEQ2_DIR + "01-htseq/.log/htseq_merge.log"
    shell:
       """(python {WATERMELON_SCRIPTS_DIR}/deseq2_htseq_merge.py \
       --htseq_dir={params.input_dir} \
       --suffix=_counts.txt \
       --counts_filename={output.counts_filename} \
       --stats_filename={output.stats_filename}
       ) 2>&1 | tee {log}"""
