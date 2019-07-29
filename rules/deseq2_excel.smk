rule deseq2_excel:
    input:
        gene = DIFFEX_DIR + "deseq2/annotated/{factor_name}/{contrast}.annot.txt",
        glossary = WATERMELON_SCRIPTS_DIR + 'deseq2_glossary.txt',
    output:
        annotated_file = DIFFEX_DIR + "deseq2/excel/{factor_name}/{contrast}.xlsx",
    log:
        DIFFEX_DIR + "deseq2/excel/.log/{factor_name}_{contrast}.diffex_excel.log",
    conda:
        'envs/deseq2_excel.yaml'
    shell:'''(
        python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py \
            -g {input.gene} \
            --glossary {input.glossary} \
            {output.annotated_file}
        ) 2>&1 | tee {log} '''
