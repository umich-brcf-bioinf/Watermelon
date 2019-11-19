rule align_pseudotrim:
    input:
        raw_fastq = ALIGNMENT_DIR + "02-gz_reads/{sample}_R{read}.fastq.gz",
    output:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_R{read}_trimmed.fastq.gz",
    log:
        ALIGNMENT_DIR + "02-cutadapt/.log/{sample}_R{read}_pseudotrim_SE.log"
    shell:
        '''(ln -sTr {input.raw_fastq} {output}
            echo No trimming done
            ) 2>&1 |tee {log} '''
