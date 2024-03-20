rule align_concat_reads:
    input:
        lambda wildcards: FASTQS_TO_CONCAT[wildcards.sample][wildcards.read]
    output:
        ALIGNMENT_DIR + "01-raw_reads/{sample}_R{read}.fastq.gz", # TWS - Despite the suffix, these may or may not be gzipped
    resources: runtime=180
    params:
        project_name = config['report_info']['project_name']
    shell:
        "cat {input} > {output}"
