rule deseq2_summary:
    input:
        input_files = expand(DESEQ2_DIR + "03-annotation/{phenotype}/{comparison}.annot.txt",
                             zip,
                             phenotype=REPLICATE_PHENOTYPE_NAMES,
                             comparison=REPLICATE_COMPARISON_GROUPS),
    output:
        summary_txt = DESEQ2_DIR + "06-summary/deseq2_summary.txt",
        summary_xlsx = DESEQ2_DIR + "06-summary/deseq2_summary.xlsx",
    log:
        DESEQ2_DIR + "06-summary/.log/deseq2_summary.log"
    conda:
        '../envs/python_3.6.1.yaml'
    params:
        output_dir = DESEQ2_DIR + "06-summary/",
    shell:
        '''(python {WATERMELON_SCRIPTS_DIR}/diffex_summary.py \
            --annotation_column gene_id \
            --annotation_null . \
            --diffex_call_column Call \
            --diffex_call_pass YES \
            --output_file {output.summary_txt} \
            --output_xlsx {output.summary_xlsx} \
            {input.input_files}
        touch {params.output_dir}) 2>&1 | tee {log} '''
