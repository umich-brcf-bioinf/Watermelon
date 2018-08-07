rule align_qc:
    input:
        sample_checksum = CONFIG_CHECKSUMS_DIR + "config-samples.watermelon.md5",
        raw_read_fastq_files = rnaseq_snakefile_helper.expand_sample_read_endedness(\
            ALIGNMENT_DIR + "03-fastqc_reads/{sample}_trimmed_{read_endedness}_fastqc.html",
            SAMPLE_READS),
        align_summary_files = expand(ALIGNMENT_DIR + "04-tophat/{sample}/{sample}_align_summary.txt",
                                     sample=config["samples"]),
        align_fastq_files = expand(ALIGNMENT_DIR + "05-fastqc_align/{sample}_accepted_hits_fastqc.html",
                                     sample=config["samples"]),
        fastq_screen_alignment = FASTQ_SCREEN_ALIGNMENT,
    output:
        ALIGNMENT_DIR + "06-qc/alignment_qc.html"
    params:
        output_dir = ALIGNMENT_DIR + "06-qc/",
        output_filename = "alignment_qc.html",
        multiqc_config_filename = WATERMELON_CONFIG_DIR + "multiqc_config.yaml",
    log:
        ALIGNMENT_DIR + "06-qc/.log/align_qc.log"
    shell:
        '''(module purge
        module load watermelon &&
        echo 'watermelon|version|multiqc|'`multiqc --version | cut -d' ' -f2-`
        multiqc --force \
            --exclude cutadapt \
            --exclude bowtie2 \
            --config {params.multiqc_config_filename} \
            --outdir {params.output_dir} \
            --filename {params.output_filename} \
            {ALIGNMENT_DIR}
        ) 2>&1 | tee {log} '''
