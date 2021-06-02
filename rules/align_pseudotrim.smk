rule align_pseudotrim:
    input:
        raw_fastq = ALIGNMENT_DIR + "02-gz_reads/{sample}_R{read}.fastq.gz",
    output:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_R{read}_trimmed.fastq.gz",
    params:
        project_name = config['report_info']['project_name']
    shell:
        "ln -sTr {input.raw_fastq} {output}"
