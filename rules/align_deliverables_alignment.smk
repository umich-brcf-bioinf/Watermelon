rule align_deliverables_alignment:
    input:
        #fastqc reads
        rnaseq_snakefile_helper.expand_sample_read_endedness(
                ALIGNMENT_DIR + "03-fastqc_reads/{sample}_trimmed_{read_endedness}_fastqc.html",
                SAMPLE_READS),
        #fastqc aligned
        expand(ALIGNMENT_DIR + "05-fastqc_align/{sample}.genome_fastqc.html",
                sample=config["samples"]),
        #combined count matrices (gene-level only for now)
        combined_counts = expand(ALIGNMENT_DIR + "06-annotate_combined_counts/gene_{type}.annot.txt",
            type=['FPKM', 'TPM']),
        #multiQC
        alignment_html = ALIGNMENT_DIR + "07-qc/alignment_qc.html",
        alignment_rsem_stats = ALIGNMENT_DIR + "07-qc/alignment_qc_data/multiqc_rsem.txt"

    output:
        rnaseq_snakefile_helper.expand_sample_read_endedness(
                DELIVERABLES_DIR + "alignment/sequence_reads_fastqc/{sample}_trimmed_{read_endedness}_fastqc.html",
                SAMPLE_READS),
        expand(DELIVERABLES_DIR + "alignment/aligned_reads_fastqc/{sample}.genome_fastqc.html",
                sample=config["samples"]),
        expand(DELIVERABLES_DIR + "counts/gene_{type}.annot.txt",
            type=['FPKM', 'TPM']),
        alignment_html = DELIVERABLES_DIR + "alignment/alignment_qc.html",
        alignment_rsem_stats = DELIVERABLES_DIR + "alignment/multiqc_rsem.txt"
    params:
        raw_fastqc_input_dir    =  ALIGNMENT_DIR + "03-fastqc_reads",
        raw_fastqc_output_dir   =  DELIVERABLES_DIR + "alignment/sequence_reads_fastqc",
        align_fastqc_input_dir  =  ALIGNMENT_DIR + "05-fastqc_align",
        align_fastqc_output_dir =  DELIVERABLES_DIR + "alignment/aligned_reads_fastqc",
        combined_counts_output_dir = DELIVERABLES_DIR + "counts"
    shell:
        #For fastqc, need to copy dirs and html files
        "cp -r {params.raw_fastqc_input_dir}/* {params.raw_fastqc_output_dir} ; "
        "cp -r {params.align_fastqc_input_dir}/* {params.align_fastqc_output_dir} ; "
        "for i in {input.combined_counts} ; do cp $i {params.combined_counts_output_dir} ; done ; "
        "cp {input.alignment_html} {output.alignment_html} ; "
        "cp {input.alignment_rsem_stats} {output.alignment_rsem_stats} "
