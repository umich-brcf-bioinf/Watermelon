rule align_deliverables_fastq_screen:
    input:
        rnaseq_snakefile_helper.expand_sample_read_endedness(\
                ALIGNMENT_DIR + "03-fastq_screen/multi_species/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
        rnaseq_snakefile_helper.expand_sample_read_endedness(\
                ALIGNMENT_DIR + "03-fastq_screen/biotype/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
    output:
        rnaseq_snakefile_helper.expand_sample_read_endedness(
                DELIVERABLES_DIR + "alignment/fastq_screen/multi_species/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
        rnaseq_snakefile_helper.expand_sample_read_endedness(
                DELIVERABLES_DIR + "alignment/fastq_screen/biotype/{sample}_trimmed_{read_endedness}_screen.html",
                SAMPLE_READS),
    params:
        source_multi_species    =  ALIGNMENT_DIR + "03-fastq_screen/multi_species/",
        dest_multi_species    =  DELIVERABLES_DIR + "alignment/fastq_screen/multi_species/",
        source_biotype    =  ALIGNMENT_DIR + "03-fastq_screen/biotype/",
        dest_biotype    =  DELIVERABLES_DIR + "alignment/fastq_screen/biotype/",
    shell:
        "cp -r {params.source_multi_species}/*_screen.html {params.dest_multi_species} && "
        "cp -r {params.source_biotype}/*_screen.html {params.dest_biotype}"
