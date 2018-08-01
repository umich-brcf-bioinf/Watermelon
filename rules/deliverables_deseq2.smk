rule deliverables_deseq2:
    input:
        diffex_dir = DESEQ2_DIR + "03-deseq2_diffex",
        gene_lists = expand(DESEQ2_DIR + "06-excel/{phenotype_name}/{comparison}.xlsx",
                            zip,
                            phenotype_name=REPLICATE_PHENOTYPE_NAMES,
                            comparison=REPLICATE_COMPARISON_GROUPS),
        summary_txt = DESEQ2_DIR + "07-summary/deseq2_summary.txt",
        summary_xlsx = DESEQ2_DIR + "07-summary/deseq2_summary.xlsx",
    output:
        deliverables_dir = DELIVERABLES_DIR + "deseq2",
        summary_gene_lists = DELIVERABLES_DIR + "deseq2/gene_lists/deseq2_summary.txt",

    params:
        tmp_dir = DELIVERABLES_DIR + "deseq2.tmp",
        source_gene_list_dir = DESEQ2_DIR + "06-excel"
    shell:
        """rm -rf {output.deliverables_dir} {params.tmp_dir}
        mkdir -p {params.tmp_dir}
        cp -r {input.diffex_dir}/counts {params.tmp_dir}
        cp -r {input.diffex_dir}/plots {params.tmp_dir}
        cp -r {params.source_gene_list_dir} {params.tmp_dir}/gene_lists
        cp {input.summary_txt} {params.tmp_dir}/gene_lists
        cp {input.summary_xlsx} {params.tmp_dir}/gene_lists
        mv {params.tmp_dir} {output.deliverables_dir} """
