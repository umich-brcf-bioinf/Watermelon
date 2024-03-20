rule align_fastqc_trimmed_reads:
    input:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_R{read}_trimmed.fastq.gz"
    output:
        ALIGNMENT_DIR + "03-fastqc_reads/{sample}_R{read}_trimmed_fastqc.html"
    log:
        JOB_LOG_DIR + "align_fastqc_trimmed_reads_{sample}_R{read}.log"
    container: ENV_INFO['fastqc']['image_str']
    params:
        project_name = config['report_info']['project_name'],
        fastqc_dir = ALIGNMENT_DIR + "03-fastqc_reads"
    resources: runtime=240, cpus_per_task=2
    shell:
        'fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log}'
