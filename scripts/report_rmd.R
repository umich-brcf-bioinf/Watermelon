log = file(snakemake@log[[1]], open = 'wt')
sink(log)
sink(log, type = 'message')

save(snakemake, file = snakemake@params[['snakemake_rdata']])

################################################################################

library(rmarkdown)
library(tidyverse)
library(kableExtra)
library(knitr)
library(openxlsx)

project_name = snakemake@config[['report_info']][['project_name']]
analyst_name = snakemake@config[['report_info']][['analyst_name']]
acknowledgement_text = snakemake@config[['report_info']][['acknowledgement_text']]

if(analyst_name == "Advanced Genomics Core"){
    analyst_email = "agc-datateam@umich.edu"
    doc_author = paste0("UM ", analyst_name)
} else {
    analyst_email = snakemake@config[['email']][['to']]
    doc_author = "UM Bioinformatics Core"
}


report_dir = snakemake@params[['report_dir']]
report_dir = sub("/$", "", report_dir) #Remove trailing / if there is one
output_prefix = file.path(report_dir, 'report_draft')
methods_out = file.path(report_dir, 'methods.pdf')

# if(!dir.exists(report_dir)) {
#     dir.create(report_dir, recursive = TRUE)
# }

report_rmd = snakemake@input[['report_rmd']]

methods_rmd = snakemake@input[['methods_rmd']]
methods_fig = snakemake@input[['methods_fig']]

# Two levels up from directory where report.Rmd is located, is project dir. e.g. /path/to/project/Watermelon/report/report.Rmd
project_dir = dirname(dirname(dirname(report_rmd)))

# project_dir = getwd() # This is another option to get the project dir. This is where snakemake is called from, should be project dir.


################################################################################

# Render report
rmarkdown::render(report_rmd, output_format = 'all', output_file = output_prefix, output_dir = report_dir, params = list(project_dir = project_dir))

if (!is.null(methods_rmd)){
  # Render standalone methods doc
  rmarkdown::render(methods_rmd, output_format = 'pdf_document', output_file = methods_out, output_dir = report_dir, params = list(methods_fig = methods_fig))
}


# Copy bioinformatics.csl and references_WAT.bib alongside draft report - simplifies report finalization step if they're colocated
bfx.csl = file.path(dirname(report_rmd), 'bioinformatics.csl')
refs.wat = file.path(dirname(report_rmd), 'references_WAT.bib')
copystatus = file.copy(bfx.csl, file.path(project_dir, report_dir), overwrite = TRUE)
if(!copystatus){
  stop("Copying bioinformatics.csl to report dir failed.")
}
copystatus = file.copy(refs.wat, file.path(project_dir, report_dir), overwrite = TRUE)
if(!copystatus){
  stop("Copying references_WAT.bib to report dir failed.")
}
