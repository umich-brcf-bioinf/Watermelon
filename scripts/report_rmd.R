log = file(snakemake@log[[1]], open = 'wt')
sink(log)
sink(log, type = 'message')

save(snakemake, file = snakemake@params[['snakemake_rdata']])

################################################################################

library(rmarkdown)
library(tidyverse)
library(kableExtra)
library(knitr)

project_dir = snakemake@config[['project_dir']]
analysis_dir = snakemake@config[['analysis_dir']]

load(snakemake@input[['filter_rdata']])
# Loads filter_summary and targets

setwd(file.path(analysis_dir, 'doc'))

rmarkdown::render('report.Rmd', output_format = 'all')
