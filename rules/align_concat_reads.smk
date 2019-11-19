rule align_concat_reads:
    input:
        lambda wildcards: INPUT_MANAGER.get_files_to_concat_by_filename(wildcards.sample, wildcards.read),
    output:
        ALIGNMENT_DIR + "01-raw_reads/{sample}_R{read}.fastq.gz", # TWS - Despite the suffix, these may or may not be gzipped
    shell:
        "cat {input} > {output}"
