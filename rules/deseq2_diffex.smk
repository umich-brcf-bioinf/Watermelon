rule deseq2_diffex:
    input:
        counts = ALIGNMENT_DIR + '06-stringtie/gene_count_matrix.csv',
#        sample_metadata_file = config['diffex']['sample_metadata_file'],
#        comparisons_file = config['diffex']['comparisons_file'],
    output:
        files = expand(DESEQ2_DIR + "02-deseq2_diffex/gene_lists/{phenotype}/{comparison}.txt",
                       zip,
                       phenotype=REPLICATE_PHENOTYPE_NAMES,
                       comparison=REPLICATE_COMPARISON_GROUPS),
    threads: 8
    log:
        DESEQ2_DIR + "02-deseq2_diffex/.log/deseq2_DESeq2Diffex.log"
    params:
        dir = DESEQ2_DIR + "02-deseq2_diffex",
        fold_change = config["fold_change"],
        adjusted_pvalue = config["deseq2_adjustedPValue"],
    resources:
        memoryInGb = 16
    shell:
        #TODO this script will need to be adjusted to consume the new files
        '''(module purge
        module load watermelon_dependencies/{WAT_VER}
        rm -rf {params.dir}/counts
        rm -rf {params.dir}/plots
        rm -rf {params.dir}/gene_lists
        rm -rf {params.dir}/.tmp/*
        {WATERMELON_SCRIPTS_DIR}/deseq2_diffex.R \
            -c {input.counts} \
            -m input.sample_metadata_file \
            -f input.comparisons_file \
            -o {params.dir}/.tmp \
            --foldChange={params.fold_change} \
            --adjustedPValue={params.adjusted_pvalue} \
            --threads={threads} \
            --javaMemoryInGb={resources.memoryInGb} \
            --pandocMemoryInGb={resources.memoryInGb}
        mv {params.dir}/.tmp/* {params.dir}
        touch {params.dir}
        rm -f Rplots.pdf #Some part of R generates this empty (nuisance) plot
        ) 2>&1 | tee {log}'''
