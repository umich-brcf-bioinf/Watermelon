rule align_concat_reads:
    input:
        INPUT_DIR + "{sample}/",
    output:
        ALIGNMENT_DIR + "01-raw_reads/{sample}_R{read}.fastq.gz", # TWS - Despite the suffix, these may or may not be gzipped
    params:
        source_glob = lambda wildcards: "{}{}/*_R{}*.fastq*".format(INPUT_DIR, wildcards.sample, wildcards.read[:1]),
    shell:
        "cat {params.source_glob} > {output}"
