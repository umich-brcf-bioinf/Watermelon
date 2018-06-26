rule deliverables_tuxedo:
    input:
        excel = expand(TUXEDO_DIR + "09-excel/{phenotype_name}/{comparison}.xlsx",
                       zip,
                       phenotype_name=ALL_PHENOTYPE_NAMES,
                       comparison=ALL_COMPARISON_GROUPS),
        diffex_raw_counts = expand(TUXEDO_DIR + "06-cummerbund/{phenotype_name}/{phenotype_name}_repRawCounts.txt",
                                   phenotype_name=ALL_PHENOTYPE_NAMES),
        plots = expand(TUXEDO_DIR + "06-cummerbund/{phenotype_name}/Plots",
                       phenotype_name=ALL_PHENOTYPE_NAMES),
        summary_txt = TUXEDO_DIR + "10-summary/tuxedo_summary.txt",
        summary_xlsx = TUXEDO_DIR + "10-summary/tuxedo_summary.xlsx",
    output:
        deliverables_dir = DELIVERABLES_DIR + "tuxedo",
        summary_gene_lists = DELIVERABLES_DIR + "tuxedo/gene_lists/tuxedo_summary.txt",
    params:
        tmp_dir = DELIVERABLES_DIR + "tuxedo.tmp",
        source_gene_list_dir = TUXEDO_DIR +  "/09-excel",
        source_counts_dir =  TUXEDO_DIR +  "/06-cummerbund",
        source_plots_dir =  TUXEDO_DIR +  "/06-cummerbund",
    shell:
        """rm -rf {output.deliverables_dir} {params.tmp_dir}
        mkdir -p {params.tmp_dir}/counts {params.tmp_dir}/plots
        cp `find {params.source_counts_dir} -maxdepth 2 -name '*_repRawCounts.txt'` {params.tmp_dir}/counts
        (for phenotype in `ls {params.source_plots_dir}`; do cp -r {params.source_plots_dir}/${{phenotype}}/Plots {params.tmp_dir}/plots/$phenotype; done)
        cp -r {params.source_gene_list_dir} {params.tmp_dir}/gene_lists
        cp {input.summary_txt} {params.tmp_dir}/gene_lists
        cp {input.summary_xlsx} {params.tmp_dir}/gene_lists
        mv {params.tmp_dir} {output.deliverables_dir} """
