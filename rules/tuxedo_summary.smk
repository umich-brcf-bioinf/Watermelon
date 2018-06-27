rule tuxedo_summary:
    input:
        input_files = expand(TUXEDO_DIR + "07-split/{phenotype}/{comparison}_gene.txt",
                             zip,
                             phenotype=ALL_PHENOTYPE_NAMES,
                             comparison=ALL_COMPARISON_GROUPS) +
                       expand(TUXEDO_DIR + "07-split/{phenotype}/{comparison}_isoform.txt",
                              zip,
                              phenotype=ALL_PHENOTYPE_NAMES,
                              comparison=ALL_COMPARISON_GROUPS)
    output:
        summary_txt = TUXEDO_DIR + "10-summary/tuxedo_summary.txt",
        summary_xlsx = TUXEDO_DIR + "10-summary/tuxedo_summary.xlsx",
    log:
        TUXEDO_DIR + "10-summary/.log/tuxedo_summary.log"
    params:
        output_dir = TUXEDO_DIR + "10-summary/",
    shell:
        '''(module purge && module load python/3.6.1 &&
        set -v &&
        python {WATERMELON_SCRIPTS_DIR}/diffex_summary.py \
            --annotation_column gene_id \
            --annotation_null . \
            --diffex_call_column diff_exp \
            --diffex_call_pass Yes \
            --trim_suffix .txt \
            --output_file {output.summary_txt} \
            --output_xlsx {output.summary_xlsx} \
            {input.input_files} &&
        touch {params.output_dir}
        ) 2>&1 | tee {log} '''
