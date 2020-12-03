log = file(snakemake@log[[1]], open = 'wt')
sink(log)
sink(log, type = 'message')

save(snakemake, file = snakemake@params[['snakemake_rdata']])

################################################################################

library(rmarkdown)

# Two levels up from directory where report.Rmd is located, is project dir. e.g. /path/to/project/Watermelon/report/report.Rmd
#project_dir = dirname(dirname(dirname(snakemake@input[['report_md']])))

project_dir = getwd() # This is another option to get the project dir. This is where snakemake is called from, should be project dir. Need more testing in different contexts

report_draft_fullpath = file.path(project_dir, snakemake@input[['report_md']])
report_final_path = file.path(project_dir, snakemake@output[['report_final_html']])

report_dir = file.path(project_dir, snakemake@params[['report_dir']])

rmarkdown::render(report_draft_fullpath, output_file = report_final_path, output_format = 'html_document')

# Copy final report html to deliverables
report_final_deliverable = file.path(project_dir, snakemake@output[['report_deliverable']])
if(file.exists(report_final_deliverable)){
  stop("Final report exists in deliverables folder. Will not overwrite.")
}
copystatus = file.copy(report_final_path, report_final_deliverable, overwrite = FALSE)
if(!copystatus){
  stop("Copying final report to deliverables dir failed.")
}

# Copy linked multiqc html to deliverables
mqc_copy_path = file.path(dirname(report_draft_fullpath), 'multiqc_linked_report.html')
if(file.exists(mqc_copy_path)){ # May not exist if report_from_counts was run to produce draft
  copystatus = file.copy(mqc_copy_path, dirname(report_final_deliverable), overwrite = TRUE)
  if(!copystatus){
    stop("Copying multiqc_linked_report.html to deliverables dir failed.")
  }
}

methods_pdf = file.path(report_dir, 'methods.pdf')
if(file.exists(methods_pdf)){
  copystatus = file.copy(methods_pdf, dirname(report_final_deliverable), overwrite = TRUE)
  if(!copystatus){
    stop("Copying methods.pdf to deliverables dir failed.")
  }
}