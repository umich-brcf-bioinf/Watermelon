rule align_concat_reads:
    input:
        lambda wildcards: INPUT_MANAGER.fastqs_to_concat_dict[wildcards.sample][wildcards.read],
    output:
        ALIGNMENT_DIR + "01-raw_reads/{sample}_R{read}.fastq.gz", # TWS - Despite the suffix, these may or may not be gzipped
    shell:
        "cat {input} > {output}"
