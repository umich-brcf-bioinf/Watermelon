###################
# Load universally needed libraries

library(optparse)
library(reticulate)
library(yaml)

###################
# Define functions

named.filepaths.from.dir = function(rsem.dir, file.extension) {
  fnames = list.files(rsem.dir)
  genes.fnames = fnames[grepl(file.extension, fnames)]
  sample.names = sub(file.extension, "", genes.fnames, fixed=T)
  fpathlist = file.path(rsem.dir, genes.fnames)
  names(fpathlist) = sample.names
  return(fpathlist)
}

make_ncbilink = function(id_str){
  if(is.na(id_str)){
    return(id_str)
  } else{
    url = paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", id_str)
    return(url)
  }
}

make_diffex_model_contrast_info_dfs = function(diffex_config){
  not_model_names = c('count_min_cutoff')
  model_names = names(diffex_config)[! names(diffex_config) %in% not_model_names]
  info_df_rows = list()
  for (m in model_names) {
    curr_info_row = list(
        model = diffex_config[[m]][['DESeq2']][['design']],
        linear_fold_change = diffex_config[[m]][['linear_fold_change']],
        adjustedPValue = diffex_config[[m]][['adjustedPValue']]
    )
    if(is.null(curr_info_row[['linear_fold_change']])) {
      curr_info_row[['linear_fold_change']] = NULL  # Looks weird but that's R for you. This removes its slot from list if it's NULL
    }
    info_df_rows[[m]] = curr_info_row
  }
  info_df = bind_rows(info_df_rows, .id = "model_name")
  return(info_df)
}

###################
# Setup
if(interactive()){
  # If running interactively, set opts manually here
  opt = list()
  opt$configfile = NA  # Configfile always required. Set this to config file path
  opt$markdownfile = NA  # Required if knitting or finalizing
  opt$rdatafile = NA  # Required if knitting and no_analysis == TRUE
  # Can modify opts below if needed
  opt$project_dir = getwd()
  opt$threads = 8
  opt$import_from = 'count_matrix'
  opt$no_analysis = FALSE
  opt$no_knit = FALSE
  opt$report_finalize = FALSE
} else {
  option_list = list(
    make_option(c("-c", "--configfile"), action="store", default=NA, type='character',
                help="Name of config file (required)."),
    make_option(c("-d", "--project_dir"), action="store", default=getwd(), type='character', help="Project directory. Defaults to current working directory"),
    make_option(c("-i", "--import_from"), action="store", default='count_matrix', type='character',
                help="Import from count matrix or from rsem via tximport. Options 'count_matrix' or 'tximport'. Default 'count_matrix'"),
    make_option(c("-m", "--markdownfile"), action="store", default=NA, type='character', help="R Markdown file to knit. Required if knitting or finalizing report."),
    make_option("--no_analysis", action="store_true", default=FALSE, help="Do not run analysis code."),
    make_option("--no_knit", action="store_true", default=FALSE, help="Do not knit report document."),
    make_option("--report_finalize", action="store_true", default=FALSE, help="Finalize the report (supplied by markdownfile)."),
    make_option(c("-r", "--rdatafile"), action="store", default=NA, type='character',
                help="Rdata file containing R objects (from the analysis) to include in the report. Required for, and can only be used when combining with --no_analysis."),
    make_option(c("-t", "--threads"), action="store", default=8, type='integer', help="Number of threads to use")
  )

  opt = parse_args(OptionParser(option_list=option_list))
}

if(is.na(opt$configfile)){
  stop("Required argument --configfile is missing. For help, see --help")
}

# Load in the config
config = yaml.load_file(opt$configfile)

# Some hard-coded defaults (almost never change)
SAMPLE_COLUMN = 'sample'
WAT_DIR = "Watermelon" # Expect to use copy of Watermelon in the project's folder
DIFFEX_DIR = config[['dirs']][['diffex_results']]
PLOTS_DIR = file.path(DIFFEX_DIR, 'plots_labeled_by_pheno')
COUNTS_DIR = file.path(DIFFEX_DIR, 'counts')
SUMMARY_DIR = file.path(DIFFEX_DIR, 'summary')
RUN_INFO_DIR = file.path(DIFFEX_DIR, 'run_info')
SCRIPTS_DIR = file.path(WAT_DIR, 'scripts')
REPORT_SRC_DIR = file.path(WAT_DIR, 'report')
REPORT_OUT_DIR = config[['dirs']][['report']]

if(opt$report_finalize) {
  if(is.na(opt$markdownfile)) {
    stop("Required argument --markdownfile is missing. For help, see --help")
  }
  if(opt$no_analysis || opt$no_knit) {
    warning("--report_finalize is always independent of analysis and knitting, so both --no_analysis and --no_knit are implied. Arguments --no_analysis and --no_knit will be ignored")
  }
  opt$no_analysis = TRUE
  opt$no_knit = TRUE

  # Note: If changing output_file or output_dir, also adjust deliv_rows entry to match. They aren't synced like analysis files are
  rmarkdown::render(opt$markdownfile, output_file = 'report_final.html', output_dir = REPORT_OUT_DIR, output_format = 'html_document')

}

if(!opt$no_analysis){

  ###################
  # Analysis section

  library(BiocParallel)
  library(data.table)
  library(DESeq2)
  library(openxlsx)
  library(tximport)

  # In the end, the rows are combined to create deliv_df
  deliv_rows = list()
  # Start with outputs that are same for every project
  deliv_rows[['count_data_rda']] = list(obj_name = 'count.data.list', file_name = file.path(COUNTS_DIR, 'count_data.rda'), deliverable = TRUE)
  deliv_rows[['counts_raw']] = list(obj_name = 'counts_raw', file_name = file.path(COUNTS_DIR, 'deseq2_raw_counts.txt'), deliverable = TRUE)
  deliv_rows[['counts_norm']] = list(obj_name = 'counts_norm', file_name = file.path(COUNTS_DIR, 'deseq2_depth_normalized_counts.txt'), deliverable = TRUE)
  deliv_rows[['counts_rlog']] = list(obj_name = 'counts_rlog', file_name = file.path(COUNTS_DIR, 'deseq2_rlog_normalized_counts.txt'), deliverable = TRUE)
  deliv_rows[['counts_rlog_annot']] = list(obj_name = 'counts_rlog_annot', file_name = file.path(COUNTS_DIR, 'deseq2_rlog_normalized_counts.annot.txt'), deliverable = TRUE)
  deliv_rows[['sample_heatmap_pdf']] = list(obj_name = 'sample_heatmap', file_name = file.path(PLOTS_DIR, 'SampleHeatmap.pdf'), deliverable = TRUE)
  deliv_rows[['sample_heatmap_png']] = list(obj_name = 'sample_heatmap', file_name = file.path(PLOTS_DIR, 'SampleHeatmap.png'), deliverable = TRUE)
  deliv_rows[['topVar_heatmap_pdf']] = list(obj_name = 'topVar_heatmap', file_name = file.path(PLOTS_DIR, 'Heatmap_TopVar.pdf'), deliverable = TRUE)
  deliv_rows[['topVar_heatmap_png']] = list(obj_name = 'topVar_heatmap', file_name = file.path(PLOTS_DIR, 'Heatmap_TopVar.png'), deliverable = TRUE)
  deliv_rows[['topExp_heatmap_pdf']] = list(obj_name = 'topExp_heatmap', file_name = file.path(PLOTS_DIR, 'Heatmap_TopExp.pdf'), deliverable = TRUE)
  deliv_rows[['topExp_heatmap_png']] = list(obj_name = 'topExp_heatmap', file_name = file.path(PLOTS_DIR, 'Heatmap_TopExp.png'), deliverable = TRUE)
  deliv_rows[['summary_txt']] = list(obj_name = 'summary_df', file_name = file.path(SUMMARY_DIR, 'deseq2_summary.txt'), deliverable = TRUE)
  deliv_rows[['summary_xlsx']] = list(obj_name = 'summary_df', file_name = file.path(SUMMARY_DIR, 'deseq2_summary.xlsx'), deliverable = TRUE)
  deliv_rows[['report_draft_md']] = list(file_name = file.path(REPORT_OUT_DIR, 'report_draft.md'), deliverable = FALSE)
  deliv_rows[['report_draft_html']] = list(file_name = file.path(REPORT_OUT_DIR, 'report_draft.html'), deliverable = FALSE)
  # Note: If changing the following file_names, must also adjust in report_finalize section. They aren't necessarily synced-up like analysis files are
  deliv_rows[['report_final_html']] = list(file_name = file.path(REPORT_OUT_DIR, 'report_final.html'), deliverable = TRUE)
  deliv_rows[['sw_versions']] = list(obj_name = 'sw_versions', file_name = file.path(RUN_INFO_DIR, 'env_software_versions.yaml'), deliverable = TRUE)

  # Create needed output directories - Note: more are created below during plotting steps
  dir.create(COUNTS_DIR, recursive=T, mode="775")
  dir.create(PLOTS_DIR, mode="775")
  dir.create(RUN_INFO_DIR, mode="775")
  dir.create(SUMMARY_DIR, mode="775")

  ###################
  # Counts

  # Get phenotype matrix
  sample.info.file = config[['samplesheet']]
  pdata_full = read.csv(sample.info.file, comment.char = "#", colClasses=c('character'), stringsAsFactors = FALSE)
  pdata_full$input_glob = NULL # Drop if exists - used for align_qc but not useful in diffex

  # Import data
  message('Importing count data')
  if(opt$import_from == 'tximport'){
      message('Using tximport on rsem data')
      rsem_dir = config[['rsem_dir']]
      gene.files.list = named.filepaths.from.dir(rsem_dir, ".genes.results")
      # tximport will have the correct order only if the files list matches pdata_full
      gene.files.list = gene.files.list[order(match(names(gene.files.list), pdata_full[[SAMPLE_COLUMN]]), na.last = NA)]
      txi.rsem.gene.results = tximport(gene.files.list, type = "rsem", txIn = F, txOut = F)
      #Some genes have length zero (what does this even mean?), causing issues with creating DESeqDataSet
      #Mike Love recommends changing these from 0 to 1; https://support.bioconductor.org/p/84304/#84368
      txi.rsem.gene.results$length[txi.rsem.gene.results$length == 0] = 1
      # Create DESeqDataSet from tximport object
      # Note use of no design (i.e. = ~ 1) here since we only care about counts. For diffex, will need to use actual design
      counts_dds = DESeqDataSetFromTximport(txi = txi.rsem.gene.results, colData = pdata_full, design = ~ 1)
  } else if(opt$import_from == 'count_matrix'){
      message('Using count matrix')
      # Load count data and subset based on samplesheet - Note check.names=F prevents conversion of dashes to dots (e.g. samplename 'Sample_716-AS-1')
      counts_full = read.table(config[['count_matrix']], header=TRUE, check.names=FALSE, row.names=1, sep="\t")
      samples.list = pdata_full[[SAMPLE_COLUMN]]
      stopifnot("All samples in the sample sheet must be present in the count matrix"= all(samples.list %in% names(counts_full)))
      counts_full = counts_full[,samples.list]
      # Create DESeqDataSet from matrix
      # Note use of no design (i.e. = ~ 1) here since we only care about counts. For diffex, will need to use actual design
      counts_dds = DESeqDataSetFromMatrix(countData = counts_full, colData = pdata_full, design = ~ 1)
  }

  # Define annotation_df - used for annotating rlog counts and each deseq2 comparison later
  annotation_tsv = config[['references']][['annotation_tsv']]
  annotation_df = read.delim(annotation_tsv, stringsAsFactors = FALSE)
  #The mapping table should not have duplicate values in the index
  #The onus is on whoever creates the mapping table to do it correctly
  #There are instances of biomaRt queries returning multiple results
  #These are cases where a query ID matches equally well to multiple
  #items from another database, so biomaRt returns them all. Current solution
  #is to have the conflicting attributes separated by commas

  # Remove redundant external_gene_name column if it's identical to gene_id
  # TODO: Can probably remove this if sticking with ENSEMBL's GTFs
  if(all(annotation_df$gene_id == annotation_df$external_gene_name)){
    annotation_df$external_gene_name = NULL
  }

  # Filter lowly expressed genes
  threshold = config[['diffex']][['count_min_cutoff']]
  counts_dds = counts_dds[ rowSums(counts(counts_dds)) > threshold, ]

  # Extract counts
  counts_raw = counts(counts_dds) #raw
  # TWS - need to estimate size factors before normalization can be done
  dds.withSF = estimateSizeFactors(counts_dds)
  counts_norm = counts(dds.withSF, normalized = TRUE) #raw/lib size factors
  # Log normed
  counts_rlog = rlog(counts_dds, blind = FALSE)

  #Save the dataset and extracted counts to rdata files
  count.data.list = c('counts_raw', 'counts_norm', 'counts_rlog')
  if(exists('txi.rsem.gene.results')){
      count.data.list = c(count.data.list, 'txi.rsem.gene.results')
  } else if(exists('counts_full')){
      count.data.list = c(count.data.list, 'counts_full')
  }
  save(list = count.data.list, file = deliv_rows[['count_data_rda']][['file_name']])


  # Write counts tables
  message('Writing count files')
  raw_counts_df = as.data.frame(counts_raw)
  data.table::setDT(raw_counts_df, keep.rownames = TRUE)[] #set rownames to valid column
  setnames(raw_counts_df, 'rn', 'id')
  write.table(x = raw_counts_df, file = deliv_rows[['counts_raw']][['file_name']],
      append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

  counts_norm_df = as.data.frame(counts_norm)
  data.table::setDT(counts_norm_df, keep.rownames = TRUE)[] #set rownames to valid column
  setnames(counts_norm_df, 'rn', 'id')
  write.table(x = counts_norm_df, file = deliv_rows[['counts_norm']][['file_name']],
      append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

  rld_df = as.data.frame(assay(counts_rlog))
  data.table::setDT(rld_df, keep.rownames = TRUE)[] #set rownames to valid column
  setnames(rld_df, 'rn', 'id')
  write.table(x = rld_df, file = deliv_rows[['counts_rlog']][['file_name']],
      append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

  rld_annot_df = merge(annotation_df, rld_df, by.x='gene_id', by.y='id', all.x=FALSE, all.y=TRUE, sort=TRUE)
  write.table(x = rld_annot_df, file = deliv_rows[['counts_rlog_annot']][['file_name']],
      append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

  #####################
  # Plots by Phenotype

  # Use sample column also as rownames for plotting
  rownames(pdata_full) = pdata_full$sample

  # Source the plotting functions
  source(file.path(SCRIPTS_DIR, "deseq2_plotting_fxns.R"))

  # Set up variables and data
  phenotypes = colnames(pdata_full)[grep(SAMPLE_COLUMN, colnames(pdata_full), invert=TRUE)] # Everything colname from samplesheet except SAMPLE_COLUMN
  for (p in phenotypes) {
    dir.create(file.path(PLOTS_DIR, p), mode="775")
  }

  mat = as.matrix(assay(counts_rlog))

  boxplot_title = 'Rlog normalized counts'
  boxplot_y_lab = 'log2(counts)'

  # Plotting

  message(sprintf('Plotting sample heatmap for %s', paste(phenotypes, collapse = ', ')))
  sample_heatmap = plot_sample_correlation_heatmap(mat = mat, pdata = pdata_full, factors = phenotypes, out_basename = 'SampleHeatmap')


  message(sprintf('Plotting top variably expressed genes heatmap for %s', paste(phenotypes, collapse = ', ')))
  topVar_heatmap = plot_top_variably_expressed_heatmap(mat = mat, pdata = pdata_full, factors = phenotypes, top_n = 500, out_basename = 'Heatmap_TopVar')

  message(sprintf('Plotting top expressed genes heatmap for %s', paste(phenotypes, collapse = ', ')))
  topExp_heatmap = plot_top_expressed_heatmap(mat = mat, pdata = pdata_full, factors = phenotypes, top_n = 500, out_basename = 'Heatmap_TopExp')

  # Plots for each of the phenotypes
  rlog_boxplot_list = list()
  raw_boxplot_list = list()
  pca_12_top500_list = list()
  pca_12_top500_labeled_list = list()
  pca_23_top500_list = list()
  pca_23_top500_labeled_list = list()
  pca_12_top100_list = list()
  pca_12_top100_labeled_list = list()
  pca_23_top100_list = list()
  pca_23_top100_labeled_list = list()
  scree_top500_list = list()
  scree_top100_list = list()


  for(phenotype in phenotypes) {
    
    message(sprintf('Plotting boxplots for %s', phenotype))
    boxplot_basepath = file.path(PLOTS_DIR, phenotype, 'BoxPlot_rlog')
    log2_boxplot = plot_boxplot(mat = mat, pdata = pdata_full, factor_name = phenotype, title = boxplot_title, y_label = boxplot_y_lab, out_basepath = boxplot_basepath)
    rlog_boxplot_list[[phenotype]] = log2_boxplot # In addition to writing output figure files, add to list of plots
    # Add newly created plots to deliverables
    deliv_rows[[paste0('boxplot_rlog_', phenotype, '_pdf')]] = list(obj_name = paste0('rlog_boxplot_list[[\'', phenotype, '\']]'), file_name = paste0(boxplot_basepath, '.pdf'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('boxplot_rlog_', phenotype, '_png')]] = list(obj_name = paste0('rlog_boxplot_list[[\'', phenotype, '\']]'), file_name = paste0(boxplot_basepath, '.png'), phenotype = phenotype, deliverable = TRUE)

    boxplot_basepath = file.path(PLOTS_DIR, phenotype, 'BoxPlot_raw')
    raw_boxplot = plot_boxplot(mat = log2(counts_raw), pdata = pdata_full, factor_name = phenotype, title = 'Non-normalized counts', y_label = boxplot_y_lab, out_basepath = boxplot_basepath)
    raw_boxplot_list[[phenotype]] = raw_boxplot # In addition to writing output figure files, add to list of plots
    # Add newly created plots to deliverables
    deliv_rows[[paste0('boxplot_raw_', phenotype, '_pdf')]] = list(obj_name = paste0('raw_boxplot_list[[\'', phenotype, '\']]'), file_name = paste0(boxplot_basepath, '.pdf'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('boxplot_raw_', phenotype, '_png')]] = list(obj_name = paste0('raw_boxplot_list[[\'', phenotype, '\']]'), file_name = paste0(boxplot_basepath, '.png'), phenotype = phenotype, deliverable = TRUE)
    
    # PCA top 500
    message(sprintf('Plotting PCA for %s in dim 1 and 2, top 500', phenotype))
    pca_basepath = file.path(PLOTS_DIR, phenotype, 'PCAplot_12_top500')
    pca_result_12 = compute_PCA(mat = mat, pdata = pdata_full, factor_name = phenotype, top_n = 500, dims = c('PC1','PC2'))
    log2_pca_12 = plot_PCA(compute_PCA_result = pca_result_12, out_basepath = pca_basepath)
    log2_pca_12_lab = plot_PCA(compute_PCA_result = pca_result_12, out_basepath = paste0(pca_basepath, '_labeled'), label_samples = TRUE)
    pca_12_top500_list[[phenotype]] = log2_pca_12 # In addition to writing output figure files, add to list of plots
    pca_12_top500_labeled_list[[phenotype]] = log2_pca_12_lab # In addition to writing output figure files, add to list of plots
    # Add newly created plots to deliverables
    deliv_rows[[paste0('pca_12_top500_', phenotype, '_pdf')]] = list(obj_name = paste0('pca_12_top500_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '.pdf'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('pca_12_top500_', phenotype, '_png')]] = list(obj_name = paste0('pca_12_top500_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '.png'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('pca_12_top500_', phenotype, '_labeled_pdf')]] = list(obj_name = paste0('pca_12_top500_labeled_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '_labeled.pdf'), phenotype = phenotype, deliverable = FALSE)
    deliv_rows[[paste0('pca_12_top500_', phenotype, '_labeled_png')]] = list(obj_name = paste0('pca_12_top500_labeled_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '_labeled.png'), phenotype = phenotype, deliverable = FALSE)
    message(sprintf('Plotting PCA for %s in dim 2 and 3, top 500', phenotype))
    pca_basepath = file.path(PLOTS_DIR, phenotype, 'PCAplot_23_top500')
    pca_result_23 = compute_PCA(mat = mat, pdata = pdata_full, factor_name = phenotype, top_n = 500, dims = c('PC2','PC3'))
    log2_pca_23 = plot_PCA(compute_PCA_result = pca_result_23, out_basepath = pca_basepath)
    log2_pca_23_lab = plot_PCA(compute_PCA_result = pca_result_23, out_basepath = paste0(pca_basepath, 'labeled'), label_samples = TRUE)
    pca_23_top500_list[[phenotype]] = log2_pca_23 # In addition to writing output figure files, add to list of plots
    pca_23_top500_labeled_list[[phenotype]] = log2_pca_23_lab # In addition to writing output figure files, add to list of plots
    # Add newly created plots to deliverables
    deliv_rows[[paste0('pca_23_top500_', phenotype, '_pdf')]] = list(obj_name = paste0('pca_23_top500_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '.pdf'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('pca_23_top500_', phenotype, '_png')]] = list(obj_name = paste0('pca_23_top500_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '.png'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('pca_23_top500_', phenotype, '_labeled_pdf')]] = list(obj_name = paste0('pca_23_top500_labeled_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '_labeled.pdf'), phenotype = phenotype, deliverable = FALSE)
    deliv_rows[[paste0('pca_23_top500_', phenotype, '_labeled_png')]] = list(obj_name = paste0('pca_23_top500_labeled_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '_labeled.png'), phenotype = phenotype, deliverable = FALSE)
    
    message(sprintf('Plotting scree for %s, top 500', phenotype))
    scree_basepath = file.path(PLOTS_DIR, phenotype, 'ScreePlot_top500')
    scree_plot = plot_scree(compute_PCA_result = pca_result_12, out_basepath = scree_basepath)
    scree_top500_list[[phenotype]] = scree_plot # In addition to writing output figure files, add to list of plots
    # Add newly created plots to deliverables
    deliv_rows[[paste0('scree_top500_', phenotype, '_pdf')]] = list(obj_name = paste0('scree_top500_list[[\'', phenotype, '\']]'), file_name = paste0(scree_basepath, '.pdf'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('scree_top500_', phenotype, '_png')]] = list(obj_name = paste0('scree_top500_list[[\'', phenotype, '\']]'), file_name = paste0(scree_basepath, '.png'), phenotype = phenotype, deliverable = TRUE)
    
    # PCA top 100
    message(sprintf('Plotting PCA for %s in dim 1 and 2, top 100', phenotype))
    pca_basepath = file.path(PLOTS_DIR, phenotype, 'PCAplot_12_top100')
    pca_result_12 = compute_PCA(mat = mat, pdata = pdata_full, factor_name = phenotype, top_n = 100, dims = c('PC1','PC2'))
    log2_pca_12 = plot_PCA(compute_PCA_result = pca_result_12, out_basepath = pca_basepath)
    log2_pca_12_lab = plot_PCA(compute_PCA_result = pca_result_12, out_basepath = paste0(pca_basepath, '_labeled'), label_samples = TRUE)
    pca_12_top100_list[[phenotype]] = log2_pca_12 # In addition to writing output figure files, add to list of plots
    pca_12_top100_labeled_list[[phenotype]] = log2_pca_12_lab # In addition to writing output figure files, add to list of plots
    # Add newly created plots to deliverables
    deliv_rows[[paste0('pca_12_top100_', phenotype, '_pdf')]] = list(obj_name = paste0('pca_12_top100_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '.pdf'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('pca_12_top100_', phenotype, '_png')]] = list(obj_name = paste0('pca_12_top100_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '.png'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('pca_12_top100_', phenotype, '_labeled_pdf')]] = list(obj_name = paste0('pca_12_top100_labeled_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '_labeled.pdf'), phenotype = phenotype, deliverable = FALSE)
    deliv_rows[[paste0('pca_12_top100_', phenotype, '_labeled_png')]] = list(obj_name = paste0('pca_12_top100_labeled_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '_labeled.png'), phenotype = phenotype, deliverable = FALSE)
    
    message(sprintf('Plotting PCA for %s in dim 2 and 3, top 100', phenotype))
    pca_basepath = file.path(PLOTS_DIR, phenotype, 'PCAplot_23_top100')
    pca_result_23 = compute_PCA(mat = mat, pdata = pdata_full, factor_name = phenotype, top_n = 100, dims = c('PC2','PC3'))
    log2_pca_23 = plot_PCA(compute_PCA_result = pca_result_23, out_basepath = pca_basepath)
    log2_pca_23_lab = plot_PCA(compute_PCA_result = pca_result_23, out_basepath = paste0(pca_basepath, '_labeled'), label_samples = TRUE)
    pca_23_top100_list[[phenotype]] = log2_pca_23 # In addition to writing output figure files, add to list of plots
    pca_23_top100_labeled_list[[phenotype]] = log2_pca_23_lab # In addition to writing output figure files, add to list of plots
    # Add newly created plots to deliverables
    deliv_rows[[paste0('pca_23_top100_', phenotype, '_pdf')]] = list(obj_name = paste0('pca_23_top100_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '.pdf'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('pca_23_top100_', phenotype, '_png')]] = list(obj_name = paste0('pca_23_top100_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '.png'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('pca_23_top100_', phenotype, '_labeled_pdf')]] = list(obj_name = paste0('pca_23_top100_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '_labeled.pdf'), phenotype = phenotype, deliverable = FALSE)
    deliv_rows[[paste0('pca_23_top100_', phenotype, '_labeled_png')]] = list(obj_name = paste0('pca_23_top100_list[[\'', phenotype, '\']]'), file_name = paste0(pca_basepath, '_labeled.png'), phenotype = phenotype, deliverable = FALSE)

    message(sprintf('Plotting scree for %s, top 100', phenotype))
    scree_basepath = file.path(PLOTS_DIR, phenotype, 'ScreePlot_top500')
    scree_plot = plot_scree(compute_PCA_result = pca_result_12, out_basepath = scree_basepath)
    scree_top100_list[[phenotype]] = scree_plot # In addition to writing output figure files, add to list of plots
    # Add newly created plots to deliverables
    deliv_rows[[paste0('scree_top100_', phenotype, '_pdf')]] = list(obj_name = paste0('scree_top100_list[[\'', phenotype, '\']]'), file_name = paste0(scree_basepath, '.pdf'), phenotype = phenotype, deliverable = TRUE)
    deliv_rows[[paste0('scree_top100_', phenotype, '_png')]] = list(obj_name = paste0('scree_top100_list[[\'', phenotype, '\']]'), file_name = paste0(scree_basepath, '.png'), phenotype = phenotype, deliverable = TRUE)
  }

  ##########################
  # Differential Expression

  model_names = grep("count_min_cutoff", names(config$diffex), value=TRUE, invert=TRUE)

  # Will have rows of model_name, comparison, total_count, count_diff_expressed, count_annotated, percent_annotated
  summary_rows = list()
  volcano_plot_list = list()
  pdata_subset_list = list()

  for (model_name in model_names) {
    # Create directory structure for the results of this model
    if(length(config[['diffex']][[model_name]][['contrasts']]) == 0) {
      dir.create(file.path(DIFFEX_DIR, sprintf("diffex_%s", model_name)), mode="775")  # If contrasts not specified, don't create volcano plot dir
    } else {
      dir.create(file.path(DIFFEX_DIR, sprintf("diffex_%s/volcano_plots", model_name)), recursive=TRUE, mode="775")
    }
    #######################################
    # DESeq2 Initialization for Each Model
    design = config[['diffex']][[model_name]][['DESeq2']][['design']]
    deseq.params = config[['diffex']][[model_name]][['DESeq2']][['DESeq']]
    feature_subset = config[['diffex']][[model_name]][['subset']]
    if (!is.null(feature_subset) && feature_subset != "" && feature_subset != "all") {
      feature_subset_split = unlist(strsplit(feature_subset, '::'))
      feature_label = feature_subset_split[1]
      feature_value = feature_subset_split[2]
      pdata_keeprows = which(pdata_full[[feature_label]] == feature_value)
      pdata = pdata_full[pdata_keeprows,]
      pdata_subset_list[[model_name]] = pdata
      counts = counts_full[,as.character(pdata$sample)]
      subset_key = paste0('counts_raw_subset_', model_name)
      deliv_rows[[subset_key]] = list(obj_name = paste0('pdata_subset_list[[\'', model_name, '\']]'), file_name = file.path(COUNTS_DIR, paste0('deseq2_raw_counts_subset_', model_name, '.txt')), model_name = model_name, deliverable = FALSE)
      write.table(counts,file=deliv_rows[[subset_key]][['file_name']], sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
    } else {
      pdata = pdata_full
      counts = counts_full
    }
    
    # Create DESeqDataSet and filter lowly expressed genes
    message('Initializing DESeq2 result')
    if(opt$import_from == 'tximport'){
      dds = DESeqDataSetFromTximport(txi = txi.rsem.gene.results, colData = pdata, design = as.formula(design))
    } else if(opt$import_from == 'count_matrix'){
      dds = DESeqDataSetFromMatrix(countData = counts, colData = pdata, design = as.formula(design))
    }
    
    # Filter lowly expressed genes
    threshold = config[['diffex']][['count_min_cutoff']]
    dds = dds[ rowSums(counts(dds)) > threshold, ]
    
    # Parse the params for DESeq call, if it can't be converted return it as string
    deseq.params.parsed = lapply(deseq.params, function(x) {
      tryCatch(eval(parse(text=x)),
              error=function(e){
                message(paste0("DESeq2 config parameter '", x, "' not parsable to known R type. Leaving as string"))
                x=as.character(x)
              })
    })
    # Add dds object to params list
    deseq.params.parsed[['object']] = dds
    # Set parallel to TRUE
    deseq.params.parsed[['parallel']] = TRUE
    
    # Set up multithreading
    multicore_param = MulticoreParam(workers = opt$threads)
    register(multicore_param, default=TRUE)
    
    # Call to DESeq
    dds = do.call(DESeq, deseq.params.parsed)
    
    #################################
    # Setup Getting Results for Each Model
    fdr_cutoff = as.numeric(config[['diffex']][[model_name]][['adjustedPValue']])
    fc_cutoff_linear = config[['diffex']][[model_name]][['linear_fold_change']]
    if(is.null(fc_cutoff_linear)) {
      fc_cutoff = NULL  # If linear_fold_change not set, use NULL here and that triggers diff behavior below
    } else {
      fc_cutoff = log2(as.numeric(fc_cutoff_linear))
    }

    # Can specify a list of factor::level pairs to set reference level in dds object
    set_ref_levels = config[['diffex']][[model_name]][['set_ref_levels']]
    if(!is.null(set_ref_levels)) {
      for(lvl_set in set_ref_levels) {
        set_parts = strsplit(lvl_set, '::')[[1]]
        set_factor = set_parts[[1]]
        set_level = set_parts[[2]]
        dds[[set_factor]] = relevel(dds[[set_factor]], ref = set_level)
      }
    }

    factor_name = config[['diffex']][[model_name]][['DESeq2']][['factor_name']]

    # Print the resultsNames of the dataset, for easier debugging
    message("resultsNames() of dds:")
    message(paste(resultsNames(dds), collapse=" "))
    
    # If contrasts is empty, instead of iterating over a set of contrasts, we'll have some different behavior triggered by
    # a contrast name of 'results' instead of 'foo_v_bar'. This will be used in this analysis script and while knitting the report
    if(length(config[['diffex']][[model_name]][['contrasts']]) == 0) {
      config[['diffex']][[model_name]][['contrasts']] = 'results'
    }
    for (contrast_name in config[['diffex']][[model_name]][['contrasts']]) {
      if(contrast_name != "results") {
        cont_split = unlist(strsplit(contrast_name, "_v_"))
        #Define the test and reference name from this
        test_name = cont_split[1] ; reference_name = cont_split[2]
        message(sprintf('Testing %s: %s vs %s', factor_name, test_name, reference_name))
      }
      
      deseq_config_keys = names(config[['diffex']][[model_name]][['DESeq2']])
      results.params = config[['diffex']][[model_name]][['DESeq2']][['results']]
      lfcShrink.params = config[['diffex']][[model_name]][['DESeq2']][['lfcShrink']]
      use_lfcShrink = FALSE # Default call is to results
      if(!is.null(results.params) && !is.null(lfcShrink.params)) {
        stop("Cannot use both 'results' and 'lfcShrink' to define results. Must choose one.")
      } else if( !('results' %in% deseq_config_keys) && is.null(lfcShrink.params)) {
        # ^ Allow results to be stated but empty (to use default args). On the other hand, lfsShrink can't be run with no args,
        # so missing or empty lfcShrink are treated the same
        stop("Must have config section 'results' or 'lfcShrink' to define results")
      } else if(!is.null(lfcShrink.params)) {
        # Move lfcShrink.params into results.params variable to avoid code duplication
        results.params = lfcShrink.params
        use_lfcShrink = TRUE # Set this true to change the object & fxn call when needed
      }
      
      # Parse the params for results {DESeq2} call, if it can't be converted return it as string
      results.params.parsed = lapply(results.params, function(x) {
        tryCatch(eval(parse(text=x)),
                error=function(e){
                  message(paste0("DESeq2 config parameter '", x, "' not parsable to known R type. Leaving as string"))
                  x=as.character(x)
                })
      })
      
      # Set parallel to TRUE
      results.params.parsed[['parallel']] = TRUE
      
      message('Pulling results from DESeq2 object')
      if(use_lfcShrink) {
        # function definition has 'dds' instead of 'object'
        results.params.parsed[['dds']] = dds
        # Grab the results with a call to lfcShrink()
        cat("Using lfcShrink")
        res = do.call(lfcShrink, results.params.parsed)
      } else {
        # Add dds object to params list
        results.params.parsed[['object']] = dds
        # Grab the results with a call to results()
        res = do.call(results, results.params.parsed)
      }
      
      
      # Order by adjusted p value
      res = res[order(res$padj),]
      
      # If using lfcShrink(), a 'stat' column is not generated. Create an empty one
      if (! 'stat' %in% colnames(res)) {res$stat = '.'}
      # Add gene_id column and move it to first column
      res$gene_id = rownames(res)
      res = res[,c('gene_id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')]
      if(contrast_name != 'results'){
        # Add Condition and Control columns
        res$Condition = test_name
        res$Control = reference_name
      }
      
      #make DEG calls and select DEGs
      res$Call = "NO"
      if(is.null(fc_cutoff)) {
        res$Call[which(res$padj <= fdr_cutoff)] = 'YES'  # If fc_cutoff is NULL, call DE without considering fold change
        res$log2FoldChange = NULL  # Also remove log2FoldChange column from results
      } else {
        res$Call[which(res$padj <= fdr_cutoff & abs(res$log2FoldChange) >= fc_cutoff)] = 'YES'
      }

      res = res[order(-rank(res$Call), res$pvalue), ]
      
      #write to individual tab-delimited txt files
      message('Writing diffex genes')
      res_basepath = file.path(DIFFEX_DIR, sprintf("diffex_%s/%s", model_name, contrast_name))
      write.table(
        x = res,
        file = paste0(res_basepath, ".txt"),
        append = FALSE, sep = '\t', na = '.', row.names = FALSE, quote = FALSE)
      # Also add to deliverables
      deliv_rows[[paste("results", model_name, contrast_name, sep="_")]] = list(obj_name = 'res', file_name = paste0(res_basepath, '.txt'), model_name = model_name, contrast_name = contrast_name, deliverable = TRUE)
      
      ##################
      # Add annotations
      
      annotated_results = merge(as.data.frame(res), annotation_df, by='gene_id', all.x=TRUE, all.y=FALSE, sort=FALSE)
      keep_colnames = c('gene_id', 'entrezgene_id', 'external_gene_name', 'description', 'baseMean','log2FoldChange','lfcSE','stat','pvalue','padj', 'Condition', 'Control', 'Call')
      if(contrast_name == 'results'){
        keep_colnames = keep_colnames[! keep_colnames %in% c('Condition', 'Control')]  # Remove Condition and Control columns
      }
      if(is.null(fc_cutoff)){
        keep_colnames = keep_colnames[keep_colnames != 'log2FoldChange']  # Remove log2FoldChange column
      }

      annotated_results = annotated_results[, keep_colnames]

      # Append row of details for summary DF, which will be written to file at end of analysis script
      curr_summary_row = list(
        model_name = sub("^model_", "", model_name),
        comparison = contrast_name,
        total_count = nrow(res),
        count_diff_expressed = sum(annotated_results$Call == 'YES'),
        count_annotated = sum(!is.na(annotated_results$entrezgene_id))
      )
      summary_rows[[length(summary_rows) + 1]] = curr_summary_row
      
      # Write the annotated gene list - plain text
      message('Writing annotated diffex genes')
      write.table(
        x = annotated_results,
        file = paste0(res_basepath, '.annot.txt'),
        append = FALSE, sep = '\t', na = '.', row.names = FALSE, quote = FALSE)
      # Also add to deliverables
      deliv_rows[[paste("results_annot", model_name, contrast_name, sep="_")]] = list(obj_name = 'annotated_results', file_name = paste0(res_basepath, '.annot.txt'), model_name = model_name, contrast_name = contrast_name, deliverable = TRUE)
      
      # Write the annotated gene list - xlsx
      message('Writing annotated diffex genes xlsx')
      # Load glossary for insertion into the workbook
      glossary = read.delim(file.path(SCRIPTS_DIR, "deseq2_glossary.txt"), stringsAsFactors = FALSE)
      # Create an appropriate xlsx workbook
      diffex_wb = createWorkbook()
      # Add annotated results as worksheet
      addWorksheet(diffex_wb, 'diffex_genes')
      # Make ncbi links from the Entrezgene_ids
      ncbi_links = sapply(annotated_results$entrezgene_id, make_ncbilink)
      annotated_results$ncbi_link = ncbi_links
      names(annotated_results$ncbi_link) = annotated_results$entrezgene_id # This doesn't seem to work as described
      class(annotated_results$ncbi_link) = "hyperlink"
      # Add the annotated results to the diffex_genes worksheet
      writeData(diffex_wb, 'diffex_genes', annotated_results)
      # Add glossary as worksheet
      addWorksheet(diffex_wb, 'glossary')
      writeData(diffex_wb, 'glossary', glossary)
      # Write the workbook to an xlsx file
      saveWorkbook(
        diffex_wb,
        paste0(res_basepath, '.annot.xlsx'),
        overwrite = TRUE)
      # Also add to deliverables
      deliv_rows[[paste("results_xlsx", model_name, contrast_name, sep="_")]] = list(obj_name = 'diffex_wb', file_name = paste0(res_basepath, '.annot.xlsx'), model_name = model_name, contrast_name = contrast_name, deliverable = TRUE)
      
      if(contrast_name != 'results'){
        ######################
        # Create volcano plot
        message(sprintf('Plotting volcano plot for %s %s', factor_name, contrast_name))
        
        volcano_basepath = file.path(DIFFEX_DIR, sprintf("diffex_%s/volcano_plots/VolcanoPlot_%s", model_name, contrast_name))
        
        volcano_plot = plot_volcano(
          de_list = annotated_results,
          method = "deseq2",
          exp_name = test_name,
          con_name = reference_name,
          fdr_cutoff = fdr_cutoff,
          logfc_cutoff = fc_cutoff,
          out_basepath = volcano_basepath)
        # In addition to the pdf and png listed above, also add the plot to a list of volcano plots
        volcano_key = paste("model", model_name, contrast_name, "volcano", sep="_")
        volcano_plot_list[[volcano_key]] <- volcano_plot # We may want to alter naming based on how it's used in report Rmd
        # Add newly created plots to deliverables
        deliv_rows[[paste0(volcano_key, '_pdf')]] = list(obj_name = paste0('volcano_plot_list[[\'', volcano_key, '\']]'), file_name = paste0(volcano_basepath, '.pdf'), model_name = model_name, contrast_name = contrast_name, deliverable = TRUE)
        deliv_rows[[paste0(volcano_key, '_png')]] = list(obj_name = paste0('volcano_plot_list[[\'', volcano_key, '\']]'), file_name = paste0(volcano_basepath, '.png'), model_name = model_name, contrast_name = contrast_name, deliverable = TRUE)
      }
    } # End iteration over contrasts for each model
  } # End iteration over models


  # Write summary files
  summary_df = bind_rows(summary_rows)
  summary_df$percent_annotated = round((summary_df$count_annotated / summary_df$total_count * 100), digits = 2)

  message('Writing summary table')
  write.table(
    x = summary_df,
    file=deliv_rows[['summary_txt']][['file_name']],
    append = FALSE, sep = '\t', na = '.', row.names = FALSE, quote = FALSE)

  message('Writing summary xlsx')
  # Load glossary for insertion into the workbook
  glossary = read.delim(file.path(SCRIPTS_DIR, "deseq2_glossary.txt"), stringsAsFactors = FALSE)
  # Create an appropriate xlsx workbook
  summary_wb = createWorkbook()
  addWorksheet(summary_wb, 'Sheet1')
  # Add the annotated results to the diffex_genes worksheet
  writeData(summary_wb, 'Sheet1', summary_df)
  # Write the workbook to an xlsx file
  saveWorkbook(
    summary_wb,
    deliv_rows[['summary_xlsx']][['file_name']],
    overwrite = TRUE)

  # Create info_df object from diffex section of config, to be later used in report
  info_df = make_diffex_model_contrast_info_dfs(config[['diffex']])

  # Create deliv_df object from deliv_rows list, to be later used when moving deliverables
  deliv_df = bind_rows(deliv_rows)

  # Save all objects from analysis to an Rdata file
  # Exclude `opt` to prevent unintentional override later
  env_to_save = ls()[ls() != 'opt']
  save(list=env_to_save, file=file.path(DIFFEX_DIR, "deseq2_analysis.Rdata"))
}

if(!opt$no_knit){

  ###################
  # Knitting Section

  library(rmarkdown)
  library(tidyverse)
  library(kableExtra)
  library(knitr)
  library(yaml)

  if(is.na(opt$markdownfile)) {stop("Required argument --markdownfile is missing. For help, see --help")}
  if(opt$no_analysis) {
    if(is.na(opt$rdatafile)) {
      stop("Required argument --rdatafile is missing. For help, see --help")
    } else {
      message(paste0("Loading ", opt$rdatafile, "..."))
      load(opt$rdatafile)
    }
  } else {
    if(!is.na(opt$rdatafile)) {stop("--rdatafile should only be used with --no_analysis")}
  }

  if(!dir.exists(REPORT_OUT_DIR)) {dir.create(REPORT_OUT_DIR, mode = "775")}

  #For recreating deliverables_run_info functionality
  #Do this before knitting report
  py_run_string("import os")
  py_run_string("WORKFLOW_BASEDIR = os.path.abspath('Watermelon')")
  source_python(file.path(WAT_DIR, 'version_info.smk'))
  sw_versions = VER_INFO[c('watermelon', 'WAT_diffex')]
  con = file(deliv_rows[['sw_versions']][['file_name']])
  writeLines(c(
    "#This file contains software version information in the format:",
    "#environment:",
    "#    software_package:",
    "#        software_version"
  ),con)
  writeLines(yaml::as.yaml(sw_versions), con)
  close(con)

  #Reporting set-up
  project_name = config[['report_info']][['project_name']]
  analyst_name = config[['report_info']][['analyst_name']]
  acknowledgement_text = config[['report_info']][['acknowledgement_text']]

  if(analyst_name == "Advanced Genomics Core"){
      analyst_email = "agc-datateam@umich.edu"
      doc_author = paste0("UM ", analyst_name)
  } else {
      analyst_email = config[['email']][['to']]
      doc_author = "UM Bioinformatics Core"
  }

  project_dir = opt$project_dir

  ################################################################################

  diffex_annot_file = '%s/diffex_%s/%s.annot.txt'
  diffex_summary_file = '%s/summary/deseq2_summary.txt'
  diffex_volcano_file = '%s/diffex_%s/volcano_plots/VolcanoPlot_%s.png'
  qc_boxplot_file = '%s/plots_labeled_by_pheno/%s/BoxPlot_%s.png'
  qc_heatmap_file = '%s/plots_labeled_by_pheno/SampleHeatmap.png'
  qc_pca_file = '%s/plots_labeled_by_pheno/%s/PCAplot_12_%s.png'

  ################################################################################

  output_prefix = fs::path_ext_remove(basename(deliv_rows[['report_draft_md']][['file_name']]))
  rmarkdown::render(opt$markdownfile, output_dir = REPORT_OUT_DIR, output_format = 'all', output_file = output_prefix)

  # Copy bioinformatics.csl and references_WAT.bib alongside draft report - simplifies report finalization step if they're colocated
  bfx.csl = file.path(REPORT_SRC_DIR, 'bioinformatics.csl')
  refs.wat = file.path(REPORT_SRC_DIR, 'references_WAT.bib')

  copystatus = file.copy(bfx.csl, file.path(project_dir, REPORT_OUT_DIR), overwrite = TRUE)
  if(!copystatus){
    stop("Copying bioinformatics.csl to report dir failed.")
  }
  copystatus = file.copy(refs.wat, file.path(project_dir, REPORT_OUT_DIR), overwrite = TRUE)
  if(!copystatus){
    stop("Copying references_WAT.bib to report dir failed.")
  }

  # Create deliverables_list.txt - something to feed into rsync command
  deliv_vec = deliv_df %>% filter(deliverable == TRUE) %>% pull(file_name)
  # Add some /./'s for rsync's dest path creation (re-create everything after /./)
  deliv_vec = sub("diffex_results/", "diffex_results/./", deliv_vec)
  deliv_vec = sub("/report", "/./report", deliv_vec)
  writeLines(deliv_vec, "deliverables_list.txt")

}