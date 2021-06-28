rule align_deliverables_fastq_screen:
    input:
        msp = expand(ALIGNMENT_DIR + "03-fastq_screen/multi_species/{basename}_trimmed_screen.html",
            basename=helper.gather_basenames(SAMPLE_BNAMES, SAMPLESHEET.index)),
        biot = expand(ALIGNMENT_DIR + "03-fastq_screen/biotype/{basename}_trimmed_screen.html",
            basename=helper.gather_basenames(SAMPLE_BNAMES, SAMPLESHEET.index))
    output:
        expand(DELIVERABLES_DIR + "fastq_screen/{screen_type}/{basename}_trimmed_screen.html",
            screen_type=['multi_species', 'biotype'],
            basename=helper.gather_basenames(SAMPLE_BNAMES, SAMPLESHEET.index))
    params:
        project_name = config['report_info']['project_name'],
        dest_msp = DELIVERABLES_DIR + "fastq_screen/multi_species/",
        dest_biot = DELIVERABLES_DIR + "fastq_screen/biotype/"
    shell:
        "cp {input.msp} {params.dest_msp} && "
        "cp {input.biot} {params.dest_biot}"
