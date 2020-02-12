rule align_multiqc:
    input:
        raw_read_fastqc_files = expand(ALIGNMENT_DIR + "03-fastqc_reads/{basename}_trimmed_fastqc.html",
                                    basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
                                ),
        align_summary_files = expand(ALIGNMENT_DIR + "04-rsem_star_align/{sample}.stat/{sample}.cnt",
                                     sample=config[SAMPLES_KEY]),
        align_fastq_files = expand(ALIGNMENT_DIR + "05-fastqc_align/{sample}.genome_fastqc.html",
                                   sample=config[SAMPLES_KEY]),
        fastq_screen_alignment = FASTQ_SCREEN_ALIGNMENT,
    output:
        ALIGNMENT_DIR + "07-qc/alignment_qc.html"
    params:
        output_dir = ALIGNMENT_DIR + "07-qc/",
        output_filename = "alignment_qc.html",
        multiqc_config_filename = WATERMELON_CONFIG_DIR + "multiqc_config.yaml",
    conda: 'envs/multiqc/multiqc.yaml'
    singularity: 'docker://umichbfxcore/multiqc'
    log:
        ALIGNMENT_DIR + "07-qc/.log/align_qc.log"
    shell:
        '''(multiqc --version
        multiqc --force \
            --exclude cutadapt \
            --exclude bowtie2 \
            --exclude rsem \
            --config {params.multiqc_config_filename} \
            --outdir {params.output_dir} \
            --filename {params.output_filename} \
            {ALIGNMENT_DIR}
        ) 2>&1 | tee {log} '''
