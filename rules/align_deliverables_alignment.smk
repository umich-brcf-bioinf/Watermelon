rule align_deliverables_alignment:
    input:
        #Trimmed fastqs
        expand(ALIGNMENT_DIR + "02-cutadapt/{basename}_trimmed.fastq.gz",
            basename=helper.gather_basenames(SAMPLE_BNAMES, SAMPLESHEET.index)),
        #Aligned BAMs
        expand(ALIGNMENT_DIR + "04-rsem_star_align/{sample}.genome.{ext}",
            sample=SAMPLESHEET.index, ext=["bam", "bai"]),
        #fastqc reads
        expand(ALIGNMENT_DIR + "03-fastqc_reads/{basename}_trimmed_fastqc.html",
            basename=helper.gather_basenames(SAMPLE_BNAMES, SAMPLESHEET.index)),
        #combined count matrices (gene-level only for now)
        combined_counts = expand(ALIGNMENT_DIR + "06-annotate_combined_counts/gene_{type}.annot.txt",
            type=['FPKM', 'TPM', 'expected_count']),
        #RSeQC ribosomal content table
        ribo_df = RSEQC_RIBO_DF,
        #multiQC
        alignment_html = ALIGNMENT_DIR + "07-qc/alignment_qc.html"

    output:
        expand(DELIVERABLES_DIR + "trimmed/trimmed_reads/{basename}_trimmed.fastq.gz",
            basename=helper.gather_basenames(SAMPLE_BNAMES, SAMPLESHEET.index)),
        expand(DELIVERABLES_DIR + "alignment/aligned_bams/{sample}.genome.{ext}",
            sample=SAMPLESHEET.index, ext=["bam", "bai"]),
        expand(DELIVERABLES_DIR + "trimmed/trimmed_reads_fastqc/{basename}_trimmed_fastqc.html",
            basename=helper.gather_basenames(SAMPLE_BNAMES, SAMPLESHEET.index)),
        expand(DELIVERABLES_DIR + "counts/gene_{type}.annot.txt",
            type=['FPKM', 'TPM', 'expected_count']),
        ribo_df = RSEQC_RIBO_DELIVERABLE,
        alignment_html = DELIVERABLES_DIR + "alignment/alignment_qc.html"
    resources: mem_mb=4000, runtime=480
    params:
        project_name = config['report_info']['project_name'],
        trimmed_reads_input_dir = ALIGNMENT_DIR + "02-cutadapt",
        trimmed_reads_output_dir = DELIVERABLES_DIR + "trimmed/trimmed_reads",
        aligned_bams_input_dir = ALIGNMENT_DIR + "04-rsem_star_align",
        aligned_bams_output_dir = DELIVERABLES_DIR + "alignment/aligned_bams",
        fastqc_input_dir  =  ALIGNMENT_DIR + "03-fastqc_reads",
        fastqc_output_dir =  DELIVERABLES_DIR + "trimmed/trimmed_reads_fastqc",
        combined_counts_output_dir = DELIVERABLES_DIR + "counts"
    shell:
        #Copy the fastqs and bams will take a while
        "cp {params.trimmed_reads_input_dir}/* {params.trimmed_reads_output_dir} ; "
        "cp {params.aligned_bams_input_dir}/*.genome.ba[mi] {params.aligned_bams_output_dir} ; "
        #For fastqc, need to copy dirs and html files
        "cp -r {params.fastqc_input_dir}/* {params.fastqc_output_dir} ; "
        "for i in {input.combined_counts} ; do cp $i {params.combined_counts_output_dir} ; done ; "
        "cp {input.alignment_html} {output.alignment_html} ; "
        "if [ {input.ribo_df} ] ; then cp {input.ribo_df} {output.ribo_df} ; fi "
