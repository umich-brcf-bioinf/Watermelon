rule ballgown_diffex:
    input:
        expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genes.results', sample=config[SAMPLES_KEY]),
        expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.isoforms.results', sample=config[SAMPLES_KEY])
    output:
        gene_files = expand(BALLGOWN_DIR + '01-ballgown_diffex/gene_lists/{phenotype}/{comparison}_gene.txt',
                        zip,
                        phenotype=ALL_PHENOTYPE_NAMES,
                        comparison=ALL_COMPARISON_GROUPS),
        iso_files = expand(BALLGOWN_DIR + '01-ballgown_diffex/gene_lists/{phenotype}/{comparison}_isoform.txt',
                        zip,
                        phenotype=ALL_PHENOTYPE_NAMES,
                        comparison=ALL_COMPARISON_GROUPS),
        gene_counts = BALLGOWN_DIR + '01-ballgown_diffex/counts/gene_fpkms.txt',
        iso_counts = BALLGOWN_DIR + '01-ballgown_diffex/counts/iso_fpkms.txt',
        rda = BALLGOWN_DIR + '01-ballgown_diffex/ballgown_data.rda',
    log:
        BALLGOWN_DIR + '01-ballgown_diffex/.log/ballgown_diffex.log'
    conda:
        'envs/diffex.yaml'
    params:
        rsem_dir = ALIGNMENT_DIR + '04-rsem_star_align/',
        configfile_path = CONFIGFILE_PATH
    shell:
        '''
        (Rscript {WATERMELON_SCRIPTS_DIR}/ballgown_diffex.R \
            --rsem_dir {params.rsem_dir} \
            --config_file {params.configfile_path}
        ) 2>&1 | tee {log}
        '''
