rule align_bamCoverage:
    input:
        ALIGNMENT_DIR + "04-rsem_star_align/{sample}.genome.bam"
    output:
        ALIGNMENT_DIR + "05-bamcoverage/{sample}.genome.bw"
    params:
        project_name = config['report_info']['project_name'],
    log:
        JOB_LOG_DIR + "align_bamCoverage_{sample}.log"
    container: ENV_INFO['deeptools']['image_str']
    resources: time_min=120, cpus=4
    shell:
        'bamCoverage -p {resources.cpus} -b {input} -o {output} 2>&1 | tee {log}'