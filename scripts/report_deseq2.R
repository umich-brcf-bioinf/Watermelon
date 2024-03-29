##########
# Load libraries

library(rmarkdown)
library(optparse)
library(tidyverse)
library(kableExtra)
library(knitr)

###################
# Define functions



##########
# Setup

option_list = list(
  make_option(c("-c", "--configfile"), action="store", default=NA, type='character', help="Name of config file (required)"),
  make_option(c("-m", "--markdownfile"), action="store", default=NA, type='character', help="R Markdown file to knit (required)"),
  make_option(c("-d", "--project_dir"), action="store", default=getwd(), type='character', help="Project directory. Defaults to current working directory"),
  make_option(c("-r", "--rdatafile"), action="store", default=NA, type='character', help="Rdata file containing R objects (from the analysis) to include in the report (required)")
)

opt = parse_args(OptionParser(option_list=option_list))

for(cli_option in c("configfile", "markdownfile", "rdatafile")){
  if is.na(opt[[cli_option]]){
    stop(paste0("Required argument --", cli_option, " is missing. For help, see --help"))
  }
}

# Load in the config
config = yaml.load_file(opt$configfile)

# Set up directories and file naming
DIFFEX_DIR = config[['dirs']][['diffex_results']]
PLOTS_DIR = file.path(DIFFEX_DIR, 'plots_labeled_by_pheno')
COUNTS_DIR = file.path(DIFFEX_DIR, 'counts')
SUMMARY_DIR = file.path(DIFFEX_DIR, 'summary')
REPORT_DIR = config[['dirs']][['report']]
REPORT_OUT_HTML = file.path(REPORT_DIR, 'report_final.html')
DELIVERABLES_DIR = config[['dirs']][['deliverables']]

# Create needed output directories?

# Set variables for knitting
project_dir = opt$project_dir

# Load the Rdata from the analysis
load(opt$rdatafile)

#################
# Knitting Section

rmarkdown::render(opt$markdownfile, output_file = REPORT_OUT_HTML, output_format = 'html_document')

# Copy final report html to deliverables
report_final_deliverable = file.path(project_dir, DELIVERABLES_DIR, "report_final.html")
if(file.exists(report_final_deliverable)){
  stop("Final report exists in deliverables folder. Will not overwrite.")
}
copystatus = file.copy(report_final_path, report_final_deliverable, overwrite = FALSE)
if(!copystatus){
  stop("Copying final report to deliverables dir failed.")
}
