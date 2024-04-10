rule align_fastq_screen_multi_species:
    input:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_R{read}_trimmed.fastq.gz"
    output:
        ALIGNMENT_DIR + "03-fastq_screen/multi_species/{sample}_R{read}_trimmed_screen.html",
        ALIGNMENT_DIR + "03-fastq_screen/multi_species/{sample}_R{read}_trimmed_screen.txt",
    log:
        JOB_LOG_DIR + "align_fastq_screen_multi_species_{sample}_R{read}.log"
    container: ENV_INFO['fastq_screen']['image_str']
    resources: cpus=8, mem_mb=8000
    params:
        project_name = config['report_info']['project_name'],
        aligner = config['fastq_screen']['aligner'],
        subset = config['fastq_screen']['subset'],
        multi_species_output_dir = ALIGNMENT_DIR + "03-fastq_screen/multi_species",
        multi_species_config_file = config['fastq_screen']['reference_basedir'] +"/multi_species.conf"
    shell:
        '''(fastq_screen --version
        fastq_screen \
            --threads {resources.cpus} \
            --subset {params.subset} \
            --conf {params.multi_species_config_file} \
            --outdir {params.multi_species_output_dir} \
            --aligner {params.aligner} \
            {input}
        ) 2>&1 | tee {log}'''
