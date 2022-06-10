# A rule to build the report from the R markdown
rule report_diffex:
    input:
        report_rmd = WORKFLOW_BASEDIR + '/report/report_diffex.Rmd',
        versions = DELIVERABLES_DIR + 'run_info/env_software_versions.yaml',
        diffex_summary = DIFFEX_DIR + 'summary/deseq2_summary.txt',
        diffex_annot = helper.expand_model_contrast_filenames(
                DIFFEX_DIR + 'diffex_{model_name}/{contrast}.annot.txt',
                DESEQ2_CONTRAST_DICT
            ),
        #deseq2_plots_by_phenotype
        diffex_PCA_pngs = expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{phenotype}/PCAplot_{dim}_top{ngenes}.png',
                phenotype = PHENOTYPES,
                dim = ['12','23'],
                ngenes = ['100','500']
            ),
        diffex_boxplot_pngs = expand(DIFFEX_DIR + 'plots_labeled_by_pheno/{phenotype}/BoxPlot_{transformation}.png',
                phenotype = PHENOTYPES,
                transformation=['raw', 'rlog']
            ),
        diffex_heatmap_pngs = expand(DIFFEX_DIR + 'plots_labeled_by_pheno/SampleHeatmap.png'),
        #deseq2_volcano_plots
        diffex_volcanoplot_pngs = helper.expand_model_contrast_filenames(\
            DIFFEX_DIR + 'diffex_{model_name}/volcano_plots/VolcanoPlot_{contrast}.png',
            DESEQ2_CONTRAST_DICT
            )

    output:
        report_md = REPORT_DIR + 'report_draft.md',
        report_html = REPORT_DIR + 'report_draft.html'
    log:
        JOB_LOG_DIR + 'report_from_counts.log'
    container: 'docker://umichbfxcore/wat_diffex:0.4.0'
    params:
        project_name = config['report_info']['project_name'],
        snakemake_rdata = REPORT_DIR + '.report_draft_snakemake.rda',
        report_dir = REPORT_DIR,
        diffex_dir = DIFFEX_DIR,
        add_background = True, # TODO: These could later be moved out to config
        add_bg_samples = True,
        add_prep_description = True,
        custom_sections = '', # Rscript expects a string of comma separated values. All sections: 'methods,results,appendix'
        contrasts = DESEQ2_CONTRAST_DICT,
        phenotypes = PHENOTYPES,
        sample_phenotypes = PHENOTYPE_SAMPLE_LIST,
        diffex_model_info = DIFFEX_MODEL_INFO
    script:
        '../scripts/report_rmd.R'
