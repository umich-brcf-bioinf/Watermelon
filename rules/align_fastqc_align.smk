#TWS - Maybe we can remove this rule entirely - doesn't seem like it has a real purpose
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
    conda: 'envs/fastqc/fastqc.yaml'
    singularity: ENV_INFO['fastqc']['image_str']
    resources: time_min=240, cpus=2 # fastqc is not multithreaded, but Java spawns way too many processes, so this keeps Snakemake from overruning the process limit.
    shell:
        'fastqc {input} -o {params.fastqc_dir} 2>&1 | tee {log} '
