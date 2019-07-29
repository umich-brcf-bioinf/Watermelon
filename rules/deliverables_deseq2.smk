rule deliverables_deseq2:
    input:
        gene_lists = rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DIFFEX_DIR + "deseq2/excel/{factor_name}/{contrast}.xlsx",
            DESEQ2_CONTRAST_DICT),
        summary_txt = DIFFEX_DIR + "deseq2/summary/deseq2_summary.txt",
        summary_xlsx = DIFFEX_DIR + "deseq2/summary/deseq2_summary.xlsx",
    output:
        summary_gene_lists = DELIVERABLES_DIR + "deseq2/gene_lists/deseq2_summary.txt",
    params:
        diffex_dir = DIFFEX_DIR + "deseq2",
        final_dir = DELIVERABLES_DIR + "deseq2",
        tmp_dir = DELIVERABLES_DIR + "deseq2.tmp",
        source_gene_list_dir = DIFFEX_DIR + "deseq2/excel"
    shell:
        """rm -rf {params.final_dir} {params.tmp_dir}
        mkdir -p {params.tmp_dir}
        cp -r {params.diffex_dir}/counts/*counts.txt {params.tmp_dir}
        cp -r {params.diffex_dir}/plots {params.tmp_dir}
        cp -r {params.source_gene_list_dir} {params.tmp_dir}/gene_lists
        cp {input.summary_txt} {params.tmp_dir}/gene_lists
        cp {input.summary_xlsx} {params.tmp_dir}/gene_lists
        mv {params.tmp_dir} {params.final_dir} """
