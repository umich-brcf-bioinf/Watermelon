SAMPLES = ["Sample_01", "Sample_02"]

rule all:
    input:
        expand("01-raw_reads/{sample}_R1.fastq.gz", sample=SAMPLES)

rule concat_reads:
    input:
        "00-multiplexed_reads/{sample}/"
    output:
        "01-raw_reads/{sample}_R1.fastq.gz"
    shell:
        "cat {input}/* > {output}"
