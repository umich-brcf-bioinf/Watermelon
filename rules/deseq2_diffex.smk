rule deseq2_diffex:
    input:
        counts = ALIGNMENT_DIR + '06-stringtie/gene_count_matrix.csv',
    output:
        gene_lists = expand(DESEQ2_DIR + '01-deseq2_diffex/gene_lists/{phenotype}/{comparison}.txt',
                       zip,
                       phenotype=REPLICATE_PHENOTYPE_NAMES,
                       comparison=REPLICATE_COMPARISON_GROUPS),
        counts = expand(DESEQ2_DIR + '01-deseq2_diffex/counts/{count_type}_counts.txt'
                       count_type = ['raw','depth_normalized','rlog_normalized']),
        rda = DESEQ2_DIR + '01-deseq2_diffex/deseq2_data.rda',
    threads: 8
    log:
        DESEQ2_DIR + '01-deseq2_diffex/.log/deseq2_DESeq2Diffex.log'
    params:
        dir = DESEQ2_DIR + '01-deseq2_diffex',
        configfile_path = CONFIGFILE_PATH
    conda:
        'envs/diffex.yaml'
    shell:
        '''({WATERMELON_SCRIPTS_DIR}/deseq2_diffex.R \
            --count_file {input.counts} \
            --config_file {params.configfile_path} \
            --threads {threads}
        touch {params.dir}
        ) 2>&1 | tee {log}'''
