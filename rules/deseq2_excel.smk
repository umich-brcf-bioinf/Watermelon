rule deseq2_excel:
    input:
        gene = DESEQ2_DIR + "04-annotation/{phenotype}/{comparison}.annot.txt",
        glossary = DESEQ2_DIR + "05-run_info/glossary.txt",
        run_info = DESEQ2_DIR + "05-run_info/run_info.txt"
    output:
        annotated_file = DESEQ2_DIR + "06-excel/{phenotype}/{comparison}.xlsx",
    log:
        DESEQ2_DIR + "06-excel/.log/{phenotype}_diffex_excel.log"
    params:
        output_dir = DESEQ2_DIR + "06-excel/",
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py "
        " -g {input.gene}"
        " --glossary {input.glossary} "
        " --info_filepath {input.run_info} "
        " {output.annotated_file} "
        " 2>&1 | tee {log} "
