rule align_deliverables_fastq_screen:
    input:
        rules.align_fastq_screen_biotype.output,
        rules.align_fastq_screen_multi_species.output
    output:
        DELIVERABLES_DIR + "alignment/fastq_screen/multi_species/{sample}_R{read}_trimmed_screen.html",
        DELIVERABLES_DIR + "alignment/fastq_screen/biotype/{sample}_R{read}_trimmed_screen.html"
    params:
        source_multi_species    =  ALIGNMENT_DIR + "03-fastq_screen/multi_species/",
        dest_multi_species    =  DELIVERABLES_DIR + "alignment/fastq_screen/multi_species/",
        source_biotype    =  ALIGNMENT_DIR + "03-fastq_screen/biotype/",
        dest_biotype    =  DELIVERABLES_DIR + "alignment/fastq_screen/biotype/",
    shell:
        "cp -r {params.source_multi_species}/*_screen.html {params.dest_multi_species} && "
        "cp -r {params.source_biotype}/*_screen.html {params.dest_biotype}"
