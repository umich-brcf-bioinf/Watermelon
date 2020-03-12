# A rule to build the report from the R markdown
rule report_rmd:
    input:
        filter_txt = filter_txt,
        filter_rdata = filter_rdata,
        report_rmd = report_rmd,
    output:
        report_md = report_md,
        report_html = report_html,
    benchmark:
        os.path.join(benchmark_dir, 'report_rmd.txt')
    log:
        os.path.join(log_dir, 'report_rmd.txt')
    params:
        snakemake_rdata = os.path.join(snakemake_rdata_dir, 'report_rmd.rda'),
    script:
        '../scripts/report_rmd.R'

# A rule to build the report from the markdown
# This final report creation will copy into the results/ directory
rule report_md:
    input:
        report_md = report_md,
        report_html = report_html,
    output:
        report_final_html = report_final_html,
    benchmark:
        os.path.join(benchmark_dir, 'report_md.txt')
    log:
        os.path.join(log_dir, 'report_md.txt')
    params:
        snakemake_rdata = os.path.join(snakemake_rdata_dir, 'report_md.rda'),
    script:
        '../scripts/report_md.R'
