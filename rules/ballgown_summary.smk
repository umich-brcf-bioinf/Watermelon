rule ballgown_summary:
    input:
        input_files = expand(BALLGOWN_DIR + '02-annotate/{phenotype}/{comparison}_gene.annot.txt',
                             zip,
                             phenotype=ALL_PHENOTYPE_NAMES,
                             comparison=ALL_COMPARISON_GROUPS) +
                       expand(BALLGOWN_DIR + '02-annotate/{phenotype}/{comparison}_isoform.annot.txt',
                              zip,
                              phenotype=ALL_PHENOTYPE_NAMES,
                              comparison=ALL_COMPARISON_GROUPS)
    output:
        summary_txt = BALLGOWN_DIR + '05-summary/ballgown_summary.txt',
        summary_xlsx = BALLGOWN_DIR + '05-summary/ballgown_summary.xlsx',
    log:
        BALLGOWN_DIR + '05-summary/.log/ballgown_summary.log'
    conda:
        '../envs/python_3.6.1.yaml'
    params:
        output_dir = BALLGOWN_DIR + '05-summary/',
    shell:
        '''(
        python {WATERMELON_SCRIPTS_DIR}/diffex_summary.py \
            --annotation_column gene_id \
            --annotation_null . \
            --diffex_call_column diff_exp \
            --diffex_call_pass Yes \
            --trim_suffix .txt \
            --output_file {output.summary_txt} \
            --output_xlsx {output.summary_xlsx} \
            {input.input_files}
        touch {params.output_dir}
        ) 2>&1 | tee {log} '''
