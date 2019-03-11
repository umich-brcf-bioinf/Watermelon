rule ballgown_diffex:
    input:
        gtfs = expand(ALIGNMENT_DIR + '06-stringtie/{sample}.gtf', sample=config[SAMPLES_KEY]),
        ballgown = expand(ALIGNMENT_DIR + '06-stringtie/ballgown/{sample}/{bgown_prefixes}.ctab', sample=config[SAMPLES_KEY], bgown_prefixes=['e2t','e_data','i2t','i_data','t_data']),
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
        'envs/ballgown_diffex.yaml'
    params:
        stringtie_dir = ALIGNMENT_DIR + '06-stringtie/',
        configfile_path = CONFIGFILE_PATH
    shell:
        '''
        (Rscript {WATERMELON_SCRIPTS_DIR}/ballgown_diffex.R \
            --stringtie_dir {params.stringtie_dir} \
            --config_file {params.configfile_path}
        ) 2>&1 | tee {log}
        '''
