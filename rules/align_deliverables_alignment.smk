rule align_deliverables_alignment:
    input:
        #"Hidden" deliverable - just needs to be stated as input even if it's not used
        ALIGNMENT_DIR + "05-combine_counts/gene_expected_count.txt",
        #Trimmed fastqs
        expand(ALIGNMENT_DIR + "02-cutadapt/{basename}_trimmed.fastq.gz",
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
        ),
        #fastqc reads
        expand(ALIGNMENT_DIR + "03-fastqc_reads/{basename}_trimmed_fastqc.html",
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
        ),
        expand(ALIGNMENT_DIR + "05-fastqc_align/{sample}.genome_fastqc.html",
                sample=config["samples"]),
        #combined count matrices (gene-level only for now)
        combined_counts = expand(ALIGNMENT_DIR + "06-annotate_combined_counts/gene_{type}.annot.txt",
            type=['FPKM', 'TPM']),
        #multiQC
        alignment_html = ALIGNMENT_DIR + "07-qc/alignment_qc.html"

    output:
        expand(DELIVERABLES_DIR + "alignment/trimmed_reads/{basename}_trimmed.fastq.gz",
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
        ),
        expand(DELIVERABLES_DIR + "alignment/sequence_reads_fastqc/{basename}_trimmed_fastqc.html",
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
        ),
        expand(DELIVERABLES_DIR + "alignment/aligned_reads_fastqc/{sample}.genome_fastqc.html",
                sample=config["samples"]),
        expand(DELIVERABLES_DIR + "counts/gene_{type}.annot.txt",
            type=['FPKM', 'TPM']),
        alignment_html = DELIVERABLES_DIR + "alignment/alignment_qc.html"
    params:
        trimmed_reads_input_dir = ALIGNMENT_DIR + "02-cutadapt",
        trimmed_reads_output_dir = DELIVERABLES_DIR + "alignment/trimmed_reads",
        raw_fastqc_input_dir    =  ALIGNMENT_DIR + "03-fastqc_reads",
        raw_fastqc_output_dir   =  DELIVERABLES_DIR + "alignment/sequence_reads_fastqc",
        align_fastqc_input_dir  =  ALIGNMENT_DIR + "05-fastqc_align",
        align_fastqc_output_dir =  DELIVERABLES_DIR + "alignment/aligned_reads_fastqc",
        combined_counts_output_dir = DELIVERABLES_DIR + "counts"
    shell:
        #For fastqc, need to copy dirs and html files
        "cp {params.trimmed_reads_input_dir}/* {params.trimmed_reads_output_dir} ; "
        "cp -r {params.raw_fastqc_input_dir}/* {params.raw_fastqc_output_dir} ; "
        "cp -r {params.align_fastqc_input_dir}/* {params.align_fastqc_output_dir} ; "
        "for i in {input.combined_counts} ; do cp $i {params.combined_counts_output_dir} ; done ; "
        "cp {input.alignment_html} {output.alignment_html} "
