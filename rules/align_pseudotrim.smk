rule align_pseudotrim:
    input:
        raw_fastq = ALIGNMENT_DIR + "02-gz_reads/{sample}_R{read}.fastq.gz",
    output:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_R{read}_trimmed.fastq.gz",
    shell:
        "ln -sTr {input.raw_fastq} {output}"
