rule align_fastq_screen_biotype:
    input:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_R{read}_trimmed.fastq.gz"
    output:
        ALIGNMENT_DIR + "03-fastq_screen/biotype/{sample}_R{read}_trimmed_screen.html",
        ALIGNMENT_DIR + "03-fastq_screen/biotype/{sample}_R{read}_trimmed_screen.txt",
    log:
        JOB_LOG_DIR + "align_fastq_screen_biotype_{sample}_R{read}_trimmed_screen.log"
    conda: 'envs/fastq_screen/fastq_screen.yaml'
    singularity: ENV_INFO['fastq_screen']['image_str']
    resources: cpus=8, mem_mb=8000
    params:
        project_name = config['report_info']['project_name'],
        aligner = config['fastq_screen']['aligner'],
        subset = config['fastq_screen']['subset'],
        biotype_output_dir = ALIGNMENT_DIR + "03-fastq_screen/biotype",
        biotype_config_file = config['fastq_screen']['reference_basedir'] +'/' + config['fastq_screen']['species'] + '.conf'
    shell:
        '''(fastq_screen --version
        fastq_screen \
            --threads {resources.cpus} \
            --subset {params.subset} \
            --conf {params.biotype_config_file} \
            --outdir {params.biotype_output_dir} \
            --aligner {params.aligner} \
            {input}
        ) 2>&1 | tee {log}'''
