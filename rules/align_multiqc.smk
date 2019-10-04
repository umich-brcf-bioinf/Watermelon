rule align_qc:
    input:
        raw_read_fastq_files = rnaseq_snakefile_helper.expand_sample_read_endedness(\
            ALIGNMENT_DIR + "03-fastqc_reads/{sample}_trimmed_{read_endedness}_fastqc.html",
            SAMPLE_READS),
        align_summary_files = expand(ALIGNMENT_DIR + "04-rsem_star_align/{sample}.stat/{sample}.cnt",
                                     sample=config["samples"]),
        align_fastq_files = expand(ALIGNMENT_DIR + "05-fastqc_align/{sample}.genome_fastqc.html",
                                   sample=config["samples"]),
        fastq_screen_alignment = FASTQ_SCREEN_ALIGNMENT,
    output:
        ALIGNMENT_DIR + "07-qc/alignment_qc.html"
    params:
        output_dir = ALIGNMENT_DIR + "07-qc/",
        output_filename = "alignment_qc.html",
        multiqc_config_filename = WATERMELON_CONFIG_DIR + "multiqc_config.yaml",
    conda:
        'envs/align_qc.yaml'
    log:
        ALIGNMENT_DIR + "07-qc/.log/align_qc.log"
    shell:
        '''(multiqc --version
        multiqc --force \
            --exclude cutadapt \
            --exclude bowtie2 \
            --config {params.multiqc_config_filename} \
            --outdir {params.output_dir} \
            --filename {params.output_filename} \
            {ALIGNMENT_DIR}
        ) 2>&1 | tee {log} '''
