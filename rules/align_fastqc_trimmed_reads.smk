rule align_fastqc_trimmed_reads:
    input:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_R{read}_trimmed.fastq.gz"
    output:
        ALIGNMENT_DIR + "03-fastqc_reads/{sample}_R{read}_trimmed_fastqc.html"
    log:
        JOB_LOG_DIR + "align_fastqc_trimmed_reads_{sample}_R{read}.log"
    conda: 'envs/fastqc/fastqc.yaml'
    singularity: 'docker://umichbfxcore/fastqc'
    params:
        fastqc_dir = ALIGNMENT_DIR + "03-fastqc_reads"
    resources: time_str='04:00:00', cpus=2
    shell:
        'fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log}'
