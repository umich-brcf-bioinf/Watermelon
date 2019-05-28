rule deseq2_diffex:
    input:
        expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.genes.results', sample=config[SAMPLES_KEY]),
        expand(ALIGNMENT_DIR + '04-rsem_star_align/{sample}.isoforms.results' sample=config[SAMPLES_KEY])
    output:
        gene_lists = expand(DESEQ2_DIR + '02-deseq2_diffex/gene_lists/{phenotype}/{comparison}.txt',
                       zip,
                       phenotype=REPLICATE_PHENOTYPE_NAMES,
                       comparison=REPLICATE_COMPARISON_GROUPS),
        counts = expand(DESEQ2_DIR + '02-deseq2_diffex/counts/{count_type}_counts.txt',
                       count_type = ['raw','depth_normalized','rlog_normalized']),
        rda = DESEQ2_DIR + '02-deseq2_diffex/deseq2_data.rda',
    threads: 8
    log:
        DESEQ2_DIR + '02-deseq2_diffex/.log/deseq2_DESeq2Diffex.log'
    params:
        rsem_dir = ALIGNMENT_DIR + '04-rsem_star_align'
        dir = DESEQ2_DIR + '02-deseq2_diffex',
        configfile_path = CONFIGFILE_PATH
    conda:
        'envs/diffex.yaml'
    shell:
        '''({WATERMELON_SCRIPTS_DIR}/deseq2_diffex.R \
            --rsem_dir {params.rsem_dir} \
            --config_file {params.configfile_path} \
            --threads {threads}
        touch {params.dir}
        ) 2>&1 | tee {log}'''
