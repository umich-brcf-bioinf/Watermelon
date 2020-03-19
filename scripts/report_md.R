log = file(snakemake@log[[1]], open = 'wt')
sink(log)
sink(log, type = 'message')

save(snakemake, file = snakemake@params[['snakemake_rdata']])

################################################################################

library(rmarkdown)

project_dir = snakemake@config[['project_dir']]
analysis_dir = snakemake@config[['analysis_dir']]

setwd(file.path(analysis_dir, 'doc'))

rmarkdown::render('report.md', output_format = 'all')

setwd(file.path(project_dir, analysis_dir))
system(sprintf('cp %s %s', snakemake@input[['report_html']], snakemake@output[['report_final_html']]))
