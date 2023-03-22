# A rule to build the report from the R markdown
rule report_align_qc:
    input:
        report_rmd = WORKFLOW_BASEDIR + '/report/report_align_qc.Rmd',
        methods_rmd = WORKFLOW_BASEDIR + '/report/methods_standalone.Rmd',
        methods_fig = WORKFLOW_BASEDIR + '/report/methods_fig.png',
        versions = DELIVERABLES_DIR + 'run_info/env_software_versions.yaml',
        multiqc_html = DELIVERABLES_DIR + 'alignment/alignment_qc.html',
        multiqc_star = ALIGNMENT_DIR + '07-qc/alignment_qc_data/multiqc_star.txt'
    output:
        methods_doc = REPORT_DIR + 'methods.pdf',
        report_md = REPORT_DIR + 'report_draft.md',
        report_html = REPORT_DIR + 'report_draft.html'
    log:
        JOB_LOG_DIR + 'report_align_only.log'
    container: 'docker://umichbfxcore/wat_diffex:0.4.0'
    params:
        project_name = config['report_info']['project_name'],
        snakemake_rdata = REPORT_DIR + '.report_draft_snakemake.rda',
        report_dir = REPORT_DIR,
        multiqc_gen_stats = ALIGNMENT_DIR + '07-qc/alignment_qc_data/multiqc_general_stats.txt',
        mqc_plots_dir = ALIGNMENT_DIR + '07-qc/multiqc_plots/png/'
    script:
        '../scripts/report_rmd.R'
