rule deseq2_excel:
    input:
        gene = DIFFEX_DIR + "deseq2/annotated/{model_name}/{contrast}.annot.txt",
        glossary = WATERMELON_SCRIPTS_DIR + 'deseq2_glossary.txt',
    output:
        annotated_file = DIFFEX_DIR + "deseq2/excel/{model_name}/{contrast}.xlsx",
    log:
        DIFFEX_DIR + "deseq2/excel/.log/{model_name}_{contrast}.diffex_excel.log",
    conda: 'envs/python3_pandas_excel/python3_pandas_excel.yaml'
    singularity: 'docker://umichbfxcore/python3_pandas_excel'
    shell:'''(
        python {WATERMELON_SCRIPTS_DIR}/diffex_excel.py \
            -g {input.gene} \
            --glossary {input.glossary} \
            {output.annotated_file}
        ) 2>&1 | tee {log} '''
