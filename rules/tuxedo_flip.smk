rule tuxedo_flip:
    input:
        gene_cuffdiff = TUXEDO_DIR + "01-cuffdiff/{pheno}/gene_exp.diff",
        isoform_cuffdiff = TUXEDO_DIR + "01-cuffdiff/{pheno}/isoform_exp.diff"
    output:
        gene_flip = TUXEDO_DIR + "02-flip/{pheno}/gene_exp.flip.diff",
        isoform_flip = TUXEDO_DIR + "02-flip/{pheno}/isoform_exp.flip.diff"
    params:
        comparisons = lambda wildcards: phenotypeManager.separated_comparisons(',')[wildcards.pheno]
    log:
        TUXEDO_DIR + "02-flip/.log/{pheno}_flip.log"
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_flip.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " {input.gene_cuffdiff} "
        " {output.gene_flip} "
        " {params.comparisons} "
        " 2>&1 | tee {log} && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_flip.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " {input.isoform_cuffdiff} "
        " {output.isoform_flip} "
        " {params.comparisons} "
        " 2>&1 | tee >>{log} "
