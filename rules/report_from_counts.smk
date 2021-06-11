# A rule to build the report from the R markdown
rule report_from_counts:
    input:
        report_rmd = WORKFLOW_BASEDIR + '/report/report.Rmd',
        versions = DELIVERABLES_DIR + 'run_info/env_software_versions.yaml',
        diffex_summary = DIFFEX_DIR + 'deseq2/summary/deseq2_summary.txt',
        diffex_annot = helper.expand_model_contrast_filenames(
                DIFFEX_DIR + 'deseq2/annotated/{model_name}/{contrast}.annot.txt',
                DESEQ2_CONTRAST_DICT
            ),
        #deseq2_plots_by_phenotype
        diffex_PCA_pngs = expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/PCAplot_{dim}_top{ngenes}.png',
                phenotype = PHENOTYPES,
                dim = ['12','23'],
                ngenes = ['100','500']
            ),
        diffex_boxplot_pngs = expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/BoxPlot_{transformation}.png',
                phenotype = PHENOTYPES,
                transformation=['raw', 'rlog']
            ),
        diffex_heatmap_pngs = expand(DIFFEX_DIR + 'deseq2/plots/by_phenotype/{phenotype}/SampleHeatmap.png',
                phenotype = PHENOTYPES
            ),
        #deseq2_comparison_plots
        diffex_volcanoplot_pngs = helper.expand_model_contrast_filenames(\
            DIFFEX_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.png',
            DESEQ2_CONTRAST_DICT
            )

    output:
        report_md = REPORT_DIR + 'report_draft.md',
        report_html = REPORT_DIR + 'report_draft.html'
    log:
        JOB_LOG_DIR + 'report_from_counts.log'
    singularity: 'docker://umichbfxcore/report:0.1.0'
    params:
        project_name = config['report_info']['project_name'],
        snakemake_rdata = REPORT_DIR + '.report_draft_snakemake.rda',
        report_dir = REPORT_DIR,
        diffex_dir = DIFFEX_DIR,
        add_background = True, # TODO: These could later be moved out to config
        add_wetlab = True,
        custom_sections = '', # Rscript expects a string of comma separated values. All sections: 'methods,results,appendix'
        contrasts = DESEQ2_CONTRAST_DICT,
        sample_phenotypes = PHENOTYPE_MANAGER.phenotype_sample_list,
        diffex_model_info = DIFFEX_MODEL_INFO
    script:
        '../scripts/report_rmd.R'
