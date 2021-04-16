rule align_deliverables_alignment:
    input:
        #Trimmed fastqs
        expand(ALIGNMENT_DIR + "02-cutadapt/{basename}_trimmed.fastq.gz",
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
        ),
        #Aligned BAMs
        expand(ALIGNMENT_DIR + "04-rsem_star_align/{sample}.genome.bam",
            sample=config["samples"]
        ),
        #fastqc reads
        expand(ALIGNMENT_DIR + "03-fastqc_reads/{basename}_trimmed_fastqc.html",
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])),
        #combined count matrices (gene-level only for now)
        combined_counts = expand(ALIGNMENT_DIR + "06-annotate_combined_counts/gene_{type}.annot.txt",
            type=['FPKM', 'TPM', 'expected_count']),
        #multiQC
        alignment_html = ALIGNMENT_DIR + "07-qc/alignment_qc.html"

    output:
        expand(DELIVERABLES_DIR + "alignment/trimmed_reads/{basename}_trimmed.fastq.gz",
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
        ),
        expand(DELIVERABLES_DIR + "alignment/aligned_bams/{sample}.genome.bam",
            sample=config["samples"]
        ),
        expand(DELIVERABLES_DIR + "alignment/trimmed_reads_fastqc/{basename}_trimmed_fastqc.html",
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])),
        expand(DELIVERABLES_DIR + "counts/gene_{type}.annot.txt",
            type=['FPKM', 'TPM', 'expected_count']
        ),
        alignment_html = DELIVERABLES_DIR + "alignment/alignment_qc.html"
    params:
        trimmed_reads_input_dir = ALIGNMENT_DIR + "02-cutadapt",
        trimmed_reads_output_dir = DELIVERABLES_DIR + "alignment/trimmed_reads",
        aligned_bams_input_dir = ALIGNMENT_DIR + "04-rsem_star_align",
        aligned_bams_output_dir = DELIVERABLES_DIR + "alignment/aligned_bams",
        align_fastqc_input_dir  =  ALIGNMENT_DIR + "03-fastqc_reads",
        align_fastqc_output_dir =  DELIVERABLES_DIR + "alignment/trimmed_reads_fastqc",
        combined_counts_output_dir = DELIVERABLES_DIR + "counts"
    shell:
        #Copy link the fastqs and bams will take a while
        "cp {params.trimmed_reads_input_dir}/* {params.trimmed_reads_output_dir} ; "
        "cp {params.aligned_bams_input_dir}/*.genome.bam {params.aligned_bams_output_dir} ; "
        #For fastqc, need to copy dirs and html files
        "cp -r {params.align_fastqc_input_dir}/* {params.align_fastqc_output_dir} ; "
        "for i in {input.combined_counts} ; do cp $i {params.combined_counts_output_dir} ; done ; "
        "cp {input.alignment_html} {output.alignment_html} "
