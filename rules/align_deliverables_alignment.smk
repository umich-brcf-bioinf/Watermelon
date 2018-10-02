ALL.append([rnaseq_snakefile_helper.expand_sample_read_endedness(
        DELIVERABLES_DIR + "alignment/sequence_reads_fastqc/{sample}_trimmed_{read_endedness}_fastqc.html",
        SAMPLE_READS),
    expand(DELIVERABLES_DIR + "alignment/aligned_reads_fastqc/{sample}_fastqc.html",
        sample=config["samples"]),
    DELIVERABLES_DIR + "alignment/alignment_qc.html"])

rule align_deliverables_alignment:
    input:
        rnaseq_snakefile_helper.expand_sample_read_endedness(
                ALIGNMENT_DIR + "03-fastqc_reads/{sample}_trimmed_{read_endedness}_fastqc.html",
                SAMPLE_READS),
        expand(ALIGNMENT_DIR + "05-fastqc_align/{sample}_fastqc.html",
                sample=config["samples"]),
        alignment_stats = ALIGNMENT_DIR + "07-qc/alignment_qc.html",
    output:
        rnaseq_snakefile_helper.expand_sample_read_endedness(
                DELIVERABLES_DIR + "alignment/sequence_reads_fastqc/{sample}_trimmed_{read_endedness}_fastqc.html",
                SAMPLE_READS),
        expand(DELIVERABLES_DIR + "alignment/aligned_reads_fastqc/{sample}_fastqc.html",
                sample=config["samples"]),
        alignment_stats = DELIVERABLES_DIR + "alignment/alignment_qc.html",
    params:
        raw_fastqc_input_dir    =  ALIGNMENT_DIR + "03-fastqc_reads",
        raw_fastqc_output_dir   =  DELIVERABLES_DIR + "alignment/sequence_reads_fastqc",
        align_fastqc_input_dir  =  ALIGNMENT_DIR + "05-fastqc_align",
        align_fastqc_output_dir =  DELIVERABLES_DIR + "alignment/aligned_reads_fastqc",
    shell:
        "cp -r {params.raw_fastqc_input_dir}/* {params.raw_fastqc_output_dir} && "
        "cp -r {params.align_fastqc_input_dir}/* {params.align_fastqc_output_dir} && "
        "cp -r {input.alignment_stats} {output.alignment_stats} "
