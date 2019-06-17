rule align_pseudotrim_SE:
    input:
        raw_fastq = ALIGNMENT_DIR + "01-raw_reads/{sample}_{read}_SE.fastq.gz",
    output:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_{read}_SE.fastq.gz",
    log:
        ALIGNMENT_DIR + "02-cutadapt/.log/{sample}_{read}.align_cutadapt_SE.log"
    shell:
        '''(ln -sf ../../{input.raw_fastq} {output}
            echo No trimming done
            ) 2>&1 |tee {log} '''
