rule deseq2_diffex:
    input:
        counts = ALIGNMENT_DIR + '06-stringtie/gene_count_matrix.csv',
        sample_metadata = DESEQ2_DIR + '01-metadata_contrasts/sample_metadata.txt',
        contrasts = DESEQ2_DIR + '01-metadata_contrasts/contrasts.txt',
    output:
        files = expand(DESEQ2_DIR + '02-deseq2_diffex/gene_lists/{phenotype}/{comparison}.txt',
                       zip,
                       phenotype=REPLICATE_PHENOTYPE_NAMES,
                       comparison=REPLICATE_COMPARISON_GROUPS),
        rda = DESEQ2_DIR + '02-deseq2_diffex/deseq2_data.rda',
    threads: 8
    log:
        DESEQ2_DIR + '02-deseq2_diffex/.log/deseq2_DESeq2Diffex.log'
    params:
        dir = DESEQ2_DIR + '02-deseq2_diffex',
        fold_change = config['fold_change'],
        adjusted_pvalue = config['deseq2_adjustedPValue'],
    conda:
        'envs/diffex.yaml'
    shell:
        '''(rm -rf {params.dir}/counts
        rm -rf {params.dir}/gene_lists
        rm -rf {params.dir}/.tmp/*
        {WATERMELON_SCRIPTS_DIR}/deseq2_diffex.R \
            -c {input.counts} \
            -m {input.sample_metadata} \
            -f {input.contrasts} \
            -o {params.dir}/.tmp \
            --foldChange={params.fold_change} \
            --adjustedPValue={params.adjusted_pvalue} \
            --threads={threads}
        mv {params.dir}/.tmp/* {params.dir}
        touch {params.dir}
        ) 2>&1 | tee {log}'''
