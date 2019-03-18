rule deseq2_excel:
    input:
        gene = DESEQ2_DIR + "03-annotation/{phenotype}/{comparison}.annot.txt",
        glossary = WATERMELON_SCRIPTS_DIR + 'deseq2_glossary.txt',
    output:
        annotated_file = DESEQ2_DIR + "05-excel/{phenotype}/{comparison}.xlsx",
    log:
        DESEQ2_DIR + "05-excel/.log/{phenotype}_{comparison}.diffex_excel.log"
    params:
        output_dir = DESEQ2_DIR + "05-excel/",
    shell:
        "module purge && module load python/3.6.1 && "
        "python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py "
        " -g {input.gene}"
        " --glossary {input.glossary} "
        " {output.annotated_file} "
        " 2>&1 | tee {log} "
