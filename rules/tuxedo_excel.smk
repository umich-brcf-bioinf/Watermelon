rule tuxedo_excel:
    input:
        gene = TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_gene.txt",
        isoform = TUXEDO_DIR + "07-split/{phenotype_name}/{comparison}_isoform.txt",
        glossary = TUXEDO_DIR + "08-run_info/glossary.txt",
        run_info = TUXEDO_DIR + "08-run_info/run_info.txt"
    output:
        TUXEDO_DIR + "09-excel/{phenotype_name}/{comparison}.xlsx"
    log:
        TUXEDO_DIR + "09-excel/.log/{phenotype_name}_diffex_excel.log"
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py "
        " -g {input.gene}"
        " -i {input.isoform}"
        " --glossary {input.glossary} "
        " --info_filepath {input.run_info} "
        " {output} "
        " 2>&1 | tee {log} "
