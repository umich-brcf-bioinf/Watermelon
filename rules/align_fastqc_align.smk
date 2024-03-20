rule align_fastqc_align:
    input:
        ALIGNMENT_DIR + "04-rsem_star_align/{sample}.genome.bam"
    output:
        ALIGNMENT_DIR + "05-fastqc_align/{sample}.genome_fastqc.html"
    params:
        project_name = config['report_info']['project_name'],
        fastqc_dir =  ALIGNMENT_DIR + "05-fastqc_align"
    log:
        JOB_LOG_DIR + "align_fastqc_align_{sample}.log"
    container: ENV_INFO['fastqc']['image_str']
    resources: runtime=240, cpus_per_task=2 # fastqc is not multithreaded, but Java spawns way too many processes, so this keeps Snakemake from overruning the process limit.
    shell:
        'fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log} '
