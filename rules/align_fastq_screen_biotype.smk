rule align_fastq_screen_biotype:
    input:
        ALIGNMENT_DIR + "02-cutadapt/{sample}_trimmed_{read_endedness}.fastq.gz"
    output:
        ALIGNMENT_DIR + "03-fastq_screen/biotype/{sample}_trimmed_{read_endedness}_screen.html",
        ALIGNMENT_DIR + "03-fastq_screen/biotype/{sample}_trimmed_{read_endedness}_screen.txt",
    log:
        ALIGNMENT_DIR + "03-fastq_screen/biotype/.log/{sample}_trimmed_{read_endedness}_screen.log"
    conda:
        '../envs/fastq_screen.yaml'
    threads:
        8
    params:
        aligner = FASTQ_SCREEN_CONFIG['aligner'],
        subset = FASTQ_SCREEN_CONFIG['subset'],
        biotype_output_dir = ALIGNMENT_DIR + "03-fastq_screen/biotype",
        biotype_config_file = FASTQ_SCREEN_CONFIG['reference_basedir'] +'/' + FASTQ_SCREEN_CONFIG['species'] + '.conf'
    shell:
        '''(fastq_screen --version
        fastq_screen \
            --threads {threads} \
            --subset {params.subset} \
            --conf {params.biotype_config_file} \
            --outdir {params.biotype_output_dir} \
            --aligner {params.aligner} \
            {input}
        ) 2>&1 | tee {log}'''
