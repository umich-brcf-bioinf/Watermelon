rule align_fastqc_trimmed_reads:
    input:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_{read_endedness}.fastq.gz"
    output:
        touch(ALIGNMENT_DIR + "03-fastqc_reads/reads_fastq.done"),
        ALIGNMENT_DIR + "03-fastqc_reads/{sample}_trimmed_{read_endedness}_fastqc.html"
    log:
        ALIGNMENT_DIR + "03-fastqc_reads/.log/{sample}_trimmed_{read_endedness}_fastqc.log"
    params:
        fastqc_dir = ALIGNMENT_DIR + "03-fastqc_reads"
    shell:
        "module purge && module load watermelon_dependencies && "
        "fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log}"
