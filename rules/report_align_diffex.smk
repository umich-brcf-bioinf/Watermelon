# A rule to build the report from the R markdown
rule report_draft:
    input:
        report_rmd = WORKFLOW_BASEDIR + '/report/report.Rmd',
        versions = DELIVERABLES_DIR + "run_info/env_software_versions.yaml",
        multiqc_html = DELIVERABLES_DIR + "alignment/alignment_qc.html",
        multiqc_gen_stats = ALIGNMENT_DIR + "07-qc/alignment_qc_data/multiqc_general_stats.txt",
        multiqc_star = ALIGNMENT_DIR + "07-qc/alignment_qc_data/multiqc_star.txt",
        diffex_summary = DIFFEX_DIR + 'deseq2/summary/deseq2_summary.txt',
        diffex_annot = rnaseq_snakefile_helper.expand_model_contrast_filenames(
                DIFFEX_DIR + 'deseq2/annotated/{model_name}/{contrast}.annot.txt',
                DESEQ2_CONTRAST_DICT
            )'',
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
        diffex_volcanoplot_pngs = rnaseq_snakefile_helper.expand_model_contrast_filenames(\
            DIFFEX_DIR + 'deseq2/plots/comparison_plots/{model_name}/VolcanoPlot_{contrast}.png',
            DESEQ2_CONTRAST_DICT
            )''

    output:
        report_md = REPORT_DIR + 'report_draft.md',
        report_html = REPORT_DIR + 'report_draft.html'
    log:
        REPORT_DIR + '.log/report_draft.log'
    singularity: 'docker://umichbfxcore/report:0.1.0'
    params:
        snakemake_rdata = REPORT_DIR + '.report_draft_snakemake.rda',
        report_dir = REPORT_DIR,
        diffex_dir = DIFFEX_DIR,
        add_custom = False, # TODO: This could later be moved out to config
        contrasts = DESEQ2_CONTRAST_DICT if config.get('diffex') else '', # Empty val here can also be used within the script to exclude diffex sections
        phenotypes = PHENOTYPES,
        diffex_model_info = DIFFEX_MODEL_INFO if config.get('diffex') else ''
    script:
        '../scripts/report_rmd.R'

# # A rule to build the report from the markdown
# # This final report creation will copy into the results/ directory
# rule report_md:
#     input:
#         report_md = report_md,
#         report_html = report_html,
#     output:
#         report_final_html = report_final_html,
#     benchmark:
#         os.path.join(benchmark_dir, 'report_md.txt')
#     log:
#         os.path.join(log_dir, 'report_md.txt')
#     params:
#         snakemake_rdata = os.path.join(snakemake_rdata_dir, 'report_md.rda'),
#     script:
#         '../scripts/report_md.R'
