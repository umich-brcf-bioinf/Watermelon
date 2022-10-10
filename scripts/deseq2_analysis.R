##########
# Load libraries

library(BiocParallel)
library(data.table)
library(DESeq2)
library(openxlsx)
library(optparse)
library(tximport)
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

##########
# Setup

option_list = list(
  make_option(c("-c", "--configfile"), action="store", default=NA, type='character',
              help="Name of config file"),
  make_option(c("-i", "--import_from"), action="store", default='count_matrix', type='character',
              help="Import from count matrix or from rsem via tximport. Options 'count_matrix' or 'tximport'. Default 'count_matrix'"),
  make_option(c("-t", "--threads"), action="store", default=8, type='integer', help="Number of threads to use")
)

opt = parse_args(OptionParser(option_list=option_list))

# Load in the config
config = yaml.load_file(opt$configfile)

# Some hard-coded defaults (almost never change)
SAMPLE_COLUMN = 'sample'
DIFFEX_DIR = config[['dirs']][['diffex_results']]
PLOTS_DIR = file.path(DIFFEX_DIR, 'plots_labeled_by_pheno')
COUNTS_DIR = file.path(DIFFEX_DIR, 'counts')
SUMMARY_DIR = file.path(DIFFEX_DIR, 'summary')
SCRIPTS_DIR = "Watermelon/scripts" # Expect to use copy of Watermelon in the project's folder
# Not a comprehensive list - only those that are the same for every project
OUTPUT_FILES = list('count_data_rda' = file.path(COUNTS_DIR, 'count_data.rda'),
                    'raw' = file.path(COUNTS_DIR, 'deseq2_raw_counts.txt'),
                    'norm' = file.path(COUNTS_DIR, 'deseq2_depth_normalized_counts.txt'),
                    'rlog' = file.path(COUNTS_DIR, 'deseq2_rlog_normalized_counts.txt'))

# Create needed output directories - Note: more are created below during plotting steps
dir.create(COUNTS_DIR, recursive=T, mode="775")
dir.create(PLOTS_DIR, mode="775")
dir.create(SUMMARY_DIR, mode="775")

#################
# Counts Section

# Get phenotype matrix
sample.info.file = config[['samplesheet']]
pdata = read.csv(sample.info.file, comment.char = "#", colClasses=c('character'), stringsAsFactors = FALSE)
pdata$input_glob = NULL # Drop if exists - used for align_qc but not useful in diffex

# Import data
message('Importing count data')
if(opt$import_from == 'tximport'){
    message('Using tximport on rsem data')
    rsem_dir = config[['rsem_dir']]
    gene.files.list = named.filepaths.from.dir(rsem_dir, ".genes.results")
    # tximport will have the correct order only if the files list matches pdata
    gene.files.list = gene.files.list[order(match(names(gene.files.list), pdata[[SAMPLE_COLUMN]]), na.last = NA)]
    txi.rsem.gene.results = tximport(gene.files.list, type = "rsem", txIn = F, txOut = F)
    #Some genes have length zero (what does this even mean?), causing issues with creating DESeqDataSet
    #Mike Love recommends changing these from 0 to 1; https://support.bioconductor.org/p/84304/#84368
    txi.rsem.gene.results$length[txi.rsem.gene.results$length == 0] = 1
    # Create DESeqDataSet from tximport object
    # Note use of no design (i.e. = ~ 1) here since we only care about counts. For diffex, will need to use actual design
    counts_dds = DESeqDataSetFromTximport(txi = txi.rsem.gene.results, colData = pdata, design = ~ 1)
} else if(opt$import_from == 'count_matrix'){
    message('Using count matrix')
    # Load count data and subset based on samplesheet - Note check.names=F prevents conversion of dashes to dots (e.g. samplename 'Sample_716-AS-1')
    counts = read.table(config[['count_matrix']], header=TRUE, check.names=FALSE, row.names=1, sep="\t",stringsAsFactors=FALSE)
    samples.list = pdata[,SAMPLE_COLUMN]
    stopifnot("Sample names differ between sample sheet and count matrix"= identical(sort(names(counts)),sort(samples.list)))
    counts = counts[,samples.list]
    # Create DESeqDataSet from matrix
    # Note use of no design (i.e. = ~ 1) here since we only care about counts. For diffex, will need to use actual design
    counts_dds = DESeqDataSetFromMatrix(countData = counts, colData = pdata, design = ~ 1)
}

# Filter lowly expressed genes
threshold = config[['diffex']][['count_min_cutoff']]
counts_dds = counts_dds[ rowSums(counts(counts_dds)) > threshold, ]

# Extract counts
raw_counts = counts(counts_dds) #raw
# TWS - need to estimate size factors before normalization can be done
dds.withSF = estimateSizeFactors(counts_dds)
norm_counts = counts(dds.withSF, normalized = TRUE) #raw/lib size factors
# Log normed
rld = rlog(counts_dds, blind = FALSE)

#Save the dataset and extracted counts to rdata files
count.data.list = c('raw_counts', 'norm_counts', 'rld')
if(exists('txi.rsem.gene.results')){
    count.data.list = c(count.data.list, 'txi.rsem.gene.results')
} else if(exists('counts')){
    count.data.list = c(count.data.list, 'counts')
}
save(list = count.data.list, file=OUTPUT_FILES[['count_data_rda']])


# Write counts tables
message('Writing count files')
raw_counts_df = as.data.frame(raw_counts)
data.table::setDT(raw_counts_df, keep.rownames = TRUE)[] #set rownames to valid column
setnames(raw_counts_df, 'rn', 'id')
write.table(x = raw_counts_df, file=OUTPUT_FILES[['raw']],
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

norm_counts_df = as.data.frame(norm_counts)
data.table::setDT(norm_counts_df, keep.rownames = TRUE)[] #set rownames to valid column
setnames(norm_counts_df, 'rn', 'id')
write.table(x = norm_counts_df, file=OUTPUT_FILES[['norm']],
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

rld_df = as.data.frame(assay(rld))
data.table::setDT(rld_df, keep.rownames = TRUE)[] #set rownames to valid column
setnames(rld_df, 'rn', 'id')
write.table(x = rld_df, file=OUTPUT_FILES[['rlog']],
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

#############################
# Plots by Phenotype Section

# Source the plotting functions
source(file.path(SCRIPTS_DIR, "deseq2_plotting_fxns.R"))
# Note: PLOTS_DIR must be defined for these functions to work

# Set up variables and data
phenotypes = colnames(pdata)[grep(SAMPLE_COLUMN, colnames(pdata), invert=TRUE)] # Everything colname from samplesheet except SAMPLE_COLUMN
for (p in phenotypes) {
  dir.create(file.path(PLOTS_DIR, p), mode="775")
}

mat = as.matrix(assay(rld))

boxplot_title = 'Rlog normalized counts'
boxplot_y_lab = 'log2(counts)'

# Plotting

message(sprintf('Plotting sample heatmap for %s', paste(phenotypes, collapse = ', ')))
log2_heatmap = plot_sample_correlation_heatmap(mat = mat, pdata = pdata, factors = phenotypes, out_basename = 'SampleHeatmap')

message(sprintf('Plotting top variably expressed genes heatmap for %s', paste(phenotypes, collapse = ', ')))
plot_top_variably_expressed_heatmap(mat = mat, pdata = pdata, factors = phenotypes, top_n = 500, out_basename = 'Heatmap_TopVar')

message(sprintf('Plotting top expressed genes heatmap for %s', paste(phenotypes, collapse = ', ')))
plot_top_expressed_heatmap(mat = mat, pdata = pdata, factors = phenotypes, top_n = 500, out_basename = 'Heatmap_TopExp')

# Plots for each of the phenotypes
for(phenotype in phenotypes) {
  
  message(sprintf('Plotting boxplots for %s', phenotype))
  log2_boxplot = plot_boxplot(mat = mat, pdata = pdata, factor_name = phenotype, title = boxplot_title, y_label = boxplot_y_lab, out_basename = 'BoxPlot_rlog')
  raw_boxplot = plot_boxplot(mat = log2(raw_counts), pdata = pdata, factor_name = phenotype, title = 'Non-normalized counts', y_label = boxplot_y_lab, out_basename = 'BoxPlot_raw')
  
  # PCA top 500
  message(sprintf('Plotting PCA for %s in dim 1 and 2, top 500', phenotype))
  pca_result_12 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 500, dims = c('PC1','PC2'))
  log2_pca_12 = plot_PCA(compute_PCA_result = pca_result_12, out_basename = 'PCAplot_12_top500')
  
  message(sprintf('Plotting PCA for %s in dim 2 and 3, top 500', phenotype))
  pca_result_23 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 500, dims = c('PC2','PC3'))
  log2_pca_23 = plot_PCA(compute_PCA_result = pca_result_23, out_basename = 'PCAplot_23_top500')
  
  message('Plotting scree, top 500')
  scree_plot = plot_scree(compute_PCA_result = pca_result_12, out_basename = 'ScreePlot_top500')
  
  # PCA top 100
  message(sprintf('Plotting PCA for %s in dim 1 and 2, top 100', phenotype))
  pca_result_12 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 100, dims = c('PC1','PC2'))
  log2_pca_12 = plot_PCA(compute_PCA_result = pca_result_12, out_basename = 'PCAplot_12_top100')
  
  message(sprintf('Plotting PCA for %s in dim 2 and 3, top 100', phenotype))
  pca_result_23 = compute_PCA(mat = mat, pdata = pdata, factor_name = phenotype, top_n = 100, dims = c('PC2','PC3'))
  log2_pca_23 = plot_PCA(compute_PCA_result = pca_result_23, out_basename = 'PCAplot_23_top100')
  
  message('Plotting scree, top 100')
  scree_plot = plot_scree(compute_PCA_result = pca_result_12, out_basename = 'ScreePlot_top100')
  
}

##################################
# Differential Expression Section

model_names = grep("count_min_cutoff", names(config$diffex), value=TRUE, invert=TRUE)

# Will have rows of model_name, comparison, total_count, count_diff_expressed, count_annotated, percent_annotated
summary_rows = list()

for (model_name in model_names) {
  # Create directory structure for the results of this model
  dir.create(file.path(DIFFEX_DIR, sprintf("diffex_%s/volcano_plots", model_name)), recursive=TRUE, mode="775")
  #######################################
  # DESeq2 Initialization for Each Model
  design = config[['diffex']][[model_name]][['DESeq2']][['design']]
  deseq2.params = config[['diffex']][[model_name]][['DESeq2']][['DESeq2']]
  
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
  deseq2.params.parsed = lapply(deseq2.params, function(x) {
    tryCatch(eval(parse(text=x)),
             error=function(e){
               message(paste0("DESeq2 config parameter '", x, "' not parsable to known R type. Leaving as string"))
               x=as.character(x)
             })
  })
  # Add dds object to params list
  deseq2.params.parsed[['object']] = dds
  # Set parallel to TRUE
  deseq2.params.parsed[['parallel']] = TRUE
  
  # Set up multithreading
  multicore_param = MulticoreParam(workers = opt$threads)
  register(multicore_param, default=TRUE)
  
  # Call to DESeq
  dds = do.call(DESeq, deseq2.params.parsed)
  
  ####################################################
  # Getting Results for Each Model

  factor_name = config[['diffex']][[model_name]][['DESeq2']][['factor_name']]
  
  for (contrast_name in config[['diffex']][[model_name]][['contrasts']]) {
    base_filename = contrast_name #Assign base filename using contrast name
    cont_split = unlist(strsplit(contrast_name, "_v_"))
    #Define the test and reference name from this
    test_name = cont_split[1] ; reference_name = cont_split[2]
    
    # Establish cutoffs
    fdr_cutoff = as.numeric(config[['diffex']][[model_name]][['adjustedPValue']])
    fc_cutoff = log2(as.numeric(config[['diffex']][[model_name]][['linear_fold_change']]))
    
    message(sprintf('Testing %s: %s vs %s', factor_name, test_name, reference_name))
    
    # Print the resultsNames of the dataset, for easier debugging
    message("resultsNames() of dds:")
    message(paste(resultsNames(dds), collapse=" "))
    
    results.params = config[['diffex']][[model_name]][['DESeq2']][['results']]
    lfcShrink.params = config[['diffex']][[model_name]][['DESeq2']][['lfcShrink']]
    use_lfcShrink = FALSE # Default call is to results
    if(is.null(results.params) && is.null(lfcShrink.params)) {
      stop("Must have config section 'results' or 'lfcShrink' to define results")
    } else if(!is.null(results.params) && !is.null(lfcShrink.params)) {
      stop("Cannot use both 'results' and 'lfcShrink' to define results. Must choose one.")
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
    # Add Condition and Control columns
    res$Condition = test_name
    res$Control = reference_name
    
    #make DEG calls and select DEGs
    res$Call = "NO"
    res$Call[which(res$padj <= fdr_cutoff & abs(res$log2FoldChange) >= fc_cutoff)] = 'YES'
    res = res[order(-rank(res$Call), res$pvalue), ]
    
    #write to individual tab-delimited txt files
    message('Writing diffex genes')
    write.table(
      x = res,
      file=file.path(DIFFEX_DIR, sprintf("diffex_%s/%s.txt", model_name, contrast_name)),
      append = FALSE, sep = '\t', na = '.', row.names = FALSE, quote = FALSE)
    
    ##################
    # Add annotations
    annotation_tsv = config[['references']][['annotation_tsv']]
    annotation_df = read.delim(annotation_tsv, stringsAsFactors = FALSE)
    #The mapping table should not have duplicate values in the index
    #The onus is on whoever creates the mapping table to do it correctly
    #There are instances of biomaRt queries returning multiple results
    #These are cases where a query ID matches equally well to multiple
    #items from another database, so biomaRt returns them all. Current solution
    #is to have the conflicting attributes separated by commas
    
    #Remove redundant external_gene_name column if it's identical to gene_id
    if(all(annotation_df$gene_id == annotation_df$external_gene_name)){
      annotation_df$external_gene_name = NULL
    }
    
    annotated_results = merge(as.data.frame(res), annotation_df, by='gene_id', all.x=TRUE, all.y=FALSE, sort=FALSE)
    annotated_results = annotated_results[,c('gene_id', 'entrezgene_id', 'external_gene_name', 'description', 'baseMean','log2FoldChange','lfcSE','stat','pvalue','padj', 'Condition', 'Control', 'Call')]

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
      file=file.path(DIFFEX_DIR, sprintf("diffex_%s/%s.annot.txt", model_name, contrast_name)),
      append = FALSE, sep = '\t', na = '.', row.names = FALSE, quote = FALSE)
    
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
      file.path(DIFFEX_DIR, sprintf("diffex_%s/%s.annot.xlsx", model_name, contrast_name)),
      overwrite = TRUE)
    
    ######################
    # Create volcano plot
    message(sprintf('Plotting volcano plot for %s %s', factor_name, contrast_name))
    
    out_pdf = file.path(DIFFEX_DIR, sprintf("diffex_%s/volcano_plots/VolcanoPlot_%s.pdf", model_name, contrast_name))
    out_png = file.path(DIFFEX_DIR, sprintf("diffex_%s/volcano_plots/VolcanoPlot_%s.png", model_name, contrast_name))
    
    volcano_plot = plot_volcano(
      de_list = annotated_results,
      method = "deseq2",
      exp_name = test_name,
      con_name = reference_name,
      fdr_cutoff = fdr_cutoff,
      logfc_cutoff = fc_cutoff,
      out_filepath_pdf = out_pdf,
      out_filepath_png = out_png)
  } # End iteration over contrasts for each model
} # End iteration over models


# Write summary files
summary_df = bind_rows(summary_rows)
summary_df$percent_annotated = round((summary_df$count_annotated / summary_df$total_count * 100), digits = 2)

message('Writing summary table')
write.table(
  x = summary_df,
  file=file.path(SUMMARY_DIR, "deseq2_summary.txt"),
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
  file.path(SUMMARY_DIR, "deseq2_summary.xlsx"),
  overwrite = TRUE)
