rule align_multiqc:
    input:
        align_summary_files = expand(ALIGNMENT_DIR + "04-rsem_star_align/{sample}.log",
                                     sample=config['samples']),
        align_fastq_files = expand(ALIGNMENT_DIR + "03-fastqc_reads/{basename}_trimmed_fastqc.html",
                                   basename=INPUT_MANAGER.gather_basenames(config['samples'])),
        fastq_screen_alignment = FASTQ_SCREEN_ALIGNMENT,
    output:
        ALIGNMENT_DIR + "07-qc/alignment_qc.html",
        ALIGNMENT_DIR + "07-qc/alignment_qc_data/multiqc_general_stats.txt",
        ALIGNMENT_DIR + "07-qc/alignment_qc_data/multiqc_star.txt",
    params:
        project_name = config['report_info']['project_name'],
        output_dir = ALIGNMENT_DIR + "07-qc/",
        output_filename = "alignment_qc.html",
        multiqc_config_filename = WATERMELON_CONFIG_DIR + "multiqc_config.yaml",
    conda: 'envs/multiqc/multiqc.yaml'
    resources: mem_mb=4000
    singularity: ENV_INFO['multiqc']['image_str']
    log:
        JOB_LOG_DIR + "align_multiqc.log"
    shell:
        '''(multiqc --version
        multiqc --force \
            --interactive \
            --exclude general_stats \
            --exclude cutadapt \
            --exclude bowtie2 \
            --exclude rsem \
            --export \
            --config {params.multiqc_config_filename} \
            --outdir {params.output_dir} \
            --filename {params.output_filename} \
            {ALIGNMENT_DIR}
        ) 2>&1 | tee {log} '''
