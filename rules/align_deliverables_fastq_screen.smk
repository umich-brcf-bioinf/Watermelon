rule align_deliverables_fastq_screen:
    input:
        msp = expand(ALIGNMENT_DIR + "03-fastq_screen/multi_species/{basename}_trimmed_screen.html",
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
        ),
        biot = expand(ALIGNMENT_DIR + "03-fastq_screen/multi_species/{basename}_trimmed_screen.html",
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
        )
    output:
        expand(DELIVERABLES_DIR + "alignment/fastq_screen/{screen_type}/{basename}_trimmed_screen.html",
            screen_type=['multi_species', 'biotype'],
            basename=INPUT_MANAGER.gather_basenames(config[SAMPLES_KEY])
        )
    params:
        dest_msp = DELIVERABLES_DIR + "alignment/fastq_screen/multi_species/",
        dest_biot = DELIVERABLES_DIR + "alignment/fastq_screen/biotype/",

    shell:
        "cp {input.msp} {params.dest_msp} && "
        "cp {input.biot} {params.dest_biot}"
