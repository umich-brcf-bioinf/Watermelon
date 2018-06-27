rule tuxedo_split:
    input:
        gene = TUXEDO_DIR + "04-annotate/{phenotype_name}/{phenotype_name}_gene.flagged.annot.txt",
        isoform = TUXEDO_DIR + "04-annotate/{phenotype_name}/{phenotype_name}_isoform.flagged.annot.txt",
    output:
        TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_gene.txt",
        TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_isoform.txt",
    params:
        output_dir = TUXEDO_DIR + "07-split/{phenotype_name}",
        user_specified_comparison_list = lambda wildcards: phenotypeManager.separated_comparisons(',')[wildcards.phenotype_name],
    log:
        TUXEDO_DIR + "07-split/.log/{phenotype_name}_tuxedo_split.log"
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_split.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " -o _gene.txt "
        " {input.gene} "
        " {params.output_dir} "
        " {params.user_specified_comparison_list} "
        " 2>&1 | tee {log} && "

        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/tuxedo_split.py "
        " --comparison_infix {COMPARISON_INFIX} "
        " -o _isoform.txt "
        " {input.isoform} "
        " {params.output_dir} "
        " {params.user_specified_comparison_list} "
        " 2>&1 | tee >>{log} "
