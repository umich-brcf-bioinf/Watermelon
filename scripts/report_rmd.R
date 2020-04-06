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

cat(getwd())

# if(!dir.exists(report_dir)) {
#     dir.create(report_dir, recursive = TRUE)
# }

report_rmd = snakemake@input[['report_rmd']]
report_rmd = sub("/ccmb/BioinfCore/ActiveProjects", "/nfs/med-bfx-activeprojects", report_rmd)

#setwd(dirname(report_rmd))

rmarkdown::render(report_rmd, output_format = 'all', output_file = output_prefix, output_dir = report_dir)
