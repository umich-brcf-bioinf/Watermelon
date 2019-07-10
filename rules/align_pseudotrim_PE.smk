rule align_pseudotrim_PE:
    input:
        raw_fastq_R1 = ALIGNMENT_DIR + "01-raw_reads/{sample}_R1_PE.fastq.gz",
        raw_fastq_R2 = ALIGNMENT_DIR + "01-raw_reads/{sample}_R2_PE.fastq.gz",
    output:
        R1 = ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_R1_PE.fastq.gz",
        R2 = ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_R2_PE.fastq.gz",
    log:
        ALIGNMENT_DIR + "02-cutadapt/.log/{sample}_cutadapt.log"
    shell:
        '''(
        ln -sf ../../{input.raw_fastq_R1} {output.R1}
        ln -sf ../../{input.raw_fastq_R2} {output.R2}
        echo No trimming done
        ) 2>&1 | tee {log}
        '''