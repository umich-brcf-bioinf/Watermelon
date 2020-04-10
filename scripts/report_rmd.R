log = file(snakemake@log[[1]], open = 'wt')
sink(log)
sink(log, type = 'message')

save(snakemake, file = snakemake@params[['snakemake_rdata']])

################################################################################

library(rmarkdown)
library(tidyverse)
library(kableExtra)
library(knitr)

report_dir = snakemake@params[['report_dir']]
output_prefix = paste0(report_dir, 'report_draft')

# if(!dir.exists(report_dir)) {
#     dir.create(report_dir, recursive = TRUE)
# }

report_rmd = snakemake@input[['report_rmd']]
report_rmd = sub("/ccmb/BioinfCore/ActiveProjects", "/nfs/med-bfx-activeprojects", report_rmd)

# Two levels up from directory where report.Rmd is located, is project dir. e.g. /path/to/project/Watermelon/report/report.Rmd
project_dir = dirname(dirname(dirname(report_rmd)))

# project_dir = getwd() # This is another option to get the project dir. This is where snakemake is called from, should be project dir.

rmarkdown::render(report_rmd, output_format = 'all', output_file = output_prefix, output_dir = report_dir, params = list(project_dir = project_dir))
