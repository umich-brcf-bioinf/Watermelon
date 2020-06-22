# A rule to build the report from the R markdown
rule report_draft:
    input:
        report_rmd = WORKFLOW_BASEDIR + '/report/report.Rmd',
        versions = DELIVERABLES_DIR + 'run_info/env_software_versions.yaml',
        multiqc_html = DELIVERABLES_DIR + 'alignment/alignment_qc.html',
        multiqc_gen_stats = ALIGNMENT_DIR + '07-qc/alignment_qc_data/multiqc_general_stats.txt',
        multiqc_star = ALIGNMENT_DIR + '07-qc/alignment_qc_data/multiqc_star.txt'
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
        add_background = False, # TODO: These could later be moved out to config
        add_custom = False,
        add_wetlab = False 
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
