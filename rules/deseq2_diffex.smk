rule deseq2_diffex:
    input:
        counts = ALIGNMENT_DIR + '06-stringtie/gene_count_matrix.csv',
        sample_metadata = DESEQ2_DIR + "01-metadata_contrasts/sample_metadata.txt",
        contrasts = DESEQ2_DIR + "01-metadata_contrasts/contrasts.txt",
    output:
        dir = DESEQ2_DIR + "02-deseq2_diffex",
        files = expand(DESEQ2_DIR + "02-deseq2_diffex/gene_lists/{phenotype}/{comparison}.txt",
                       zip,
                       phenotype=REPLICATE_PHENOTYPE_NAMES,
                       comparison=REPLICATE_COMPARISON_GROUPS),
    threads: 8
    log:
        DESEQ2_DIR + "02-deseq2_diffex/.log/deseq2_DESeq2Diffex.log"
    params:
        fold_change = config["fold_change"],
        adjusted_pvalue = config["deseq2_adjustedPValue"],
    resources:
        memoryInGb = 16
    shell:
        '''(module purge
        module load watermelon_dependencies/{WAT_VER}
        rm -rf {output.dir}/normalized_data
        rm -rf {output.dir}/plots
        rm -rf {output.dir}/gene_lists
        rm -rf {output.dir}/.tmp/*
        {WATERMELON_SCRIPTS_DIR}/deseq2_diffex.R \
            -c {input.counts} \
            -m {input.sample_metadata} \
            -f {input.contrasts} \
            -o {output.dir}/.tmp \
            --foldChange={params.fold_change} \
            --adjustedPValue={params.adjusted_pvalue} \
            --threads={threads} \
            --javaMemoryInGb={resources.memoryInGb} \
            --pandocMemoryInGb={resources.memoryInGb}
        mv {output.dir}/.tmp/* {output.dir}
        touch {output.dir}
        rm -f Rplots.pdf #Some part of R generates this empty (nuisance) plot
        ) 2>&1 | tee {log}'''
