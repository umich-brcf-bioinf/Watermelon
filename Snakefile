#SAMPLES = ["Sample_01", "Sample_02"]

configfile: "config.yaml"

rule all:
    input:
        expand("01-raw_reads/{sample}_R1.fastq.gz", sample=config["samples"]),
        expand("02-fastqc_reads/{sample}_R1_fastqc.html", sample=config["samples"]),
        expand("03-tophat/{sample}/accepted_hits.bam", sample=config["samples"])
        
rule concat_reads:
    input:
        "00-multiplexed_reads/{sample}/"
    output:
        "01-raw_reads/{sample}_R1.fastq.gz"
    shell:
        "cat {input}/* > {output}"

rule fastqc:
    input:
        "01-raw_reads/{sample}_R1.fastq.gz"
    output:
        "02-fastqc_reads/{sample}_R1_fastqc.html"
    shell:
        "fastqc {input} -o 02-fastqc_reads"

rule align:
    output:
        "03-tophat/{sample}/accepted_hits.bam",
        "03-tophat/{sample}/align_summary.txt"
    input:
        "01-raw_reads/{sample}_R1.fastq.gz"
    params:
        gtf_file = config["gtf"], 
        transcriptome_index= config["transcriptome_index"],
        bowtie2_index = config["bowtie2_index"],
        sample = lambda wildcards: wildcards.sample
    shell: "module load rnaseq && "
            "tophat -p 8 "
            " --b2-very-sensitive "
            " --no-coverage-search "
            " --library-type fr-unstranded "
            " -I 500000 "
            " -G {params.gtf_file} "
            " --transcriptome-index={params.transcriptome_index} "
            " -T "
            " --no-novel-juncs "
            " -o {params.sample} "
            " {params.bowtie2_index} "
            " {input} "
 
