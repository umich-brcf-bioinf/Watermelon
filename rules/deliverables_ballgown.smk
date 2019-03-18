rule deliverables_ballgown:
    input:
        excel = expand(BALLGOWN_DIR + '04-excel/{phenotype_name}/{comparison}.xlsx',
                       zip,
                       phenotype_name=ALL_PHENOTYPE_NAMES,
                       comparison=ALL_COMPARISON_GROUPS),
        summary_txt = BALLGOWN_DIR + '05-summary/ballgown_summary.txt',
        summary_xlsx = BALLGOWN_DIR + '05-summary/ballgown_summary.xlsx',
    output:
        summary_gene_lists = DELIVERABLES_DIR + 'ballgown/gene_lists/ballgown_summary.txt',
    params:
        final_dir = DELIVERABLES_DIR + 'ballgown',
        tmp_dir = DELIVERABLES_DIR + 'ballgown.tmp',
        source_gene_list_dir = BALLGOWN_DIR +  '/04-excel',
        source_counts_dir =  BALLGOWN_DIR + '01-ballgown_diffex/counts',
        source_plots_dir =  BALLGOWN_DIR + '01-ballgown_diffex/plots',
    shell:
        """
        rm -rf {params.final_dir} {params.tmp_dir}
        mkdir -p {params.tmp_dir}
        # mkdir -p {params.tmp_dir}/counts {params.tmp_dir}/plots
        # cp `find {params.source_counts_dir} -maxdepth 2 -name '*_repRawCounts.txt'` {params.tmp_dir}/counts
        # (for phenotype in `ls {params.source_plots_dir}`; do cp -r {params.source_plots_dir}/${{phenotype}}/Plots {params.tmp_dir}/plots/$phenotype; done)
        cp -r {params.source_gene_list_dir} {params.tmp_dir}/gene_lists
        cp {input.summary_txt} {params.tmp_dir}/gene_lists
        cp {input.summary_xlsx} {params.tmp_dir}/gene_lists
        mv {params.tmp_dir} {params.final_dir}
        """
