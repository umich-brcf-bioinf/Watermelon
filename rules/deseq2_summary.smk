rule deseq2_summary:
    input:
        input_files = helper.expand_model_contrast_filenames(\
            DIFFEX_DIR + "deseq2/annotated/{model_name}/{contrast}.annot.txt",
            DESEQ2_CONTRAST_DICT)
    output:
        summary_txt = DIFFEX_DIR + "deseq2/summary/deseq2_summary.txt",
        summary_xlsx = DIFFEX_DIR + "deseq2/summary/deseq2_summary.xlsx",
    log:
        JOB_LOG_DIR + "deseq2_summary.log"
    conda: 'envs/python3_pandas_excel/python3_pandas_excel.yaml'
    singularity: 'docker://umichbfxcore/python3_pandas_excel'
    params:
        project_name = config['report_info']['project_name'],
        output_dir = DIFFEX_DIR + "deseq2/summary/"
    shell:
        '''(python {WATERMELON_SCRIPTS_DIR}/diffex_summary.py \
            --annotation_column entrezgene_id \
            --annotation_null . \
            --diffex_call_column Call \
            --diffex_call_pass YES \
            --output_file {output.summary_txt} \
            --output_xlsx {output.summary_xlsx} \
            {input.input_files}
        touch {params.output_dir}) 2>&1 | tee {log} '''
