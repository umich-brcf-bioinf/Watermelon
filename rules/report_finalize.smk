# A rule to build the report from the markdown
# This final report creation will copy into the deliverables location
rule report_finalize:
    input:
        report_md = REPORT_DIR + 'report_draft.md',
    output:
        report_final_html = REPORT_DIR + 'report_final.html',
        report_deliverable = DELIVERABLES_DIR + 'report/report_final.html'
    log:
        JOB_LOG_DIR + 'report_finalize.log'
    singularity: 'docker://umichbfxcore/report:0.1.0'
    params:
        snakemake_rdata = REPORT_DIR + '.report_finalize_snakemake.rda',
        report_dir = REPORT_DIR,
        workflow_basedir = WORKFLOW_BASEDIR
    script:
        '../scripts/report_md.R'
