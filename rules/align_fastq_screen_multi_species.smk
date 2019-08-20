rule align_fastq_screen_multi_species:
    input:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_{read_endedness}.fastq.gz"
    output:
        ALIGNMENT_DIR + "03-fastq_screen/multi_species/{sample}_trimmed_{read_endedness}_screen.html",
        ALIGNMENT_DIR + "03-fastq_screen/multi_species/{sample}_trimmed_{read_endedness}_screen.txt",
    log:
        ALIGNMENT_DIR + "03-fastq_screen/multi_species/.log/{sample}_trimmed_{read_endedness}_screen.log"
    conda:
        'envs/fastq_screen.yaml'
    threads:
        8
    params:
        aligner = FASTQ_SCREEN_CONFIG['aligner'],
        subset = FASTQ_SCREEN_CONFIG['subset'],
        multi_species_output_dir = ALIGNMENT_DIR + "03-fastq_screen/multi_species",
        multi_species_config_file = FASTQ_SCREEN_CONFIG['reference_basedir'] +"/multi_species.conf"
    shell:
        '''(fastq_screen --version
        fastq_screen \
            --threads {threads} \
            --subset {params.subset} \
            --conf {params.multi_species_config_file} \
            --outdir {params.multi_species_output_dir} \
            --aligner {params.aligner} \
            {input}
        ) 2>&1 | tee {log}'''
