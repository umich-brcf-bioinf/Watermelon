# A rule to build the report from the R markdown
rule report_draft:
    input:
        report_rmd = WORKFLOW_BASEDIR + '/report/report.Rmd',
        versions = DELIVERABLES_DIR + "run_info/env_software_versions.yaml",
    output:
        report_md = REPORT_DIR + 'report_draft.md',
        report_html = REPORT_DIR + 'report_draft.html'
    log:
        REPORT_DIR + '.log/report_draft.log'
    singularity: 'docker://twsaari/report_env'
    params:
        snakemake_rdata = REPORT_DIR + '.report_draft_snakemake.rda',
        report_dir = REPORT_DIR,
        sample_phenotypes = PHENOTYPE_MANAGER.phenotype_sample_list if 'diffex' in config and config['diffex'] else '',

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
