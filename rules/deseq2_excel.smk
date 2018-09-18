rule deseq2_excel:
    input:
        gene = DESEQ2_DIR + "03-annotation/{phenotype}/{comparison}.annot.txt",
        glossary = DESEQ2_DIR + "04-run_info/glossary.txt",
        run_info = DESEQ2_DIR + "04-run_info/run_info.txt"
    output:
        annotated_file = DESEQ2_DIR + "05-excel/{phenotype}/{comparison}.xlsx",
    log:
        DESEQ2_DIR + "05-excel/.log/{phenotype}_{comparison}.diffex_excel.log"
    params:
        output_dir = DESEQ2_DIR + "05-excel/",
    conda:
        'envs/deseq2_excel.yaml'
    shell:'''(
        python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py \
            -g {input.gene} \
            --glossary {input.glossary} \
            --info_filepath {input.run_info} \
            {output.annotated_file}
        ) 2>&1 | tee {log} '''
