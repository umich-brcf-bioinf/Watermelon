##########
# Set up logging and save snakemake S4 object (for debugging or running manually)
log = file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')
save(snakemake, file = snakemake@params[['snakemake_rdata']])

#load("/nfs/med-bfx-activeprojects/trsaari/sandbox/20190725_test_Delono_RS1/analysis_test_Delano_RS1/diffex_results/deseq2/gene_lists/phenotype.CellState.treatment/DIO.WCLP.none_v_DIO.DCLP.none_snakemake.rda") #TWS DEBUG

#Isolate conda environment: https://github.com/conda-forge/r-base-feedstock/issues/37
#If we move away from conda in the future, we may want to remove this
.libPaths(R.home("library"))

##########
# Load libraries

lib.vector = c("BiocParallel", "data.table", "DESeq2", "GGally", "ggfortify", "ggplot2",
               "ggrepel", "RColorBrewer", "stringr", "tidyr", "openxlsx")
foo = suppressMessages(lapply(lib.vector, library, character.only=T, warn.conflicts=F, quietly=T))

##########
# Define functions

make_ncbilink = function(id_str){
  if(is.na(id_str)){
    return(id_str)
  } else{
    url = paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", id_str)
    return(url)
  }
}


plot_volcano = function(de_list, method = c('ballgown', 'deseq2'), exp_name, con_name, fdr_cutoff, logfc_cutoff, out_filepath_pdf, out_filepath_png) {

  method = match.arg(method)

  # Determine the correct column names on the basis of deseq2 or ballgown
  if(method == 'deseq2') {
    log2fc = 'log2FoldChange'
    pval = 'pvalue'
    padj = 'padj'
    id = 'gene_id'
    de_call = 'Call'
  } else {
    log2fc = 'log2fc'
    pval = 'pval'
    padj = 'qval'
    id = 'gene_id'
    de_call = 'diff_exp'
  }

  # Add direction column
  de_list$direction = 'NS'
  de_list$direction[de_list[, de_call] == 'YES' & de_list[, log2fc] <= 0] = 'Down'
  de_list$direction[de_list[, de_call] == 'YES' & de_list[, log2fc] > 0] = 'Up'

  # Determine direction labels (with number per)
  direction_table = table(de_list$direction)

  if('Up' %in% names(direction_table)) {
    up_label = sprintf('Up: %s', direction_table['Up'])
  } else {
    up_label = 'Up: 0'
  }

  if('Down' %in% names(direction_table)) {
    down_label = sprintf('Down: %s', direction_table['Down'])
  } else {
    down_label = 'Down: 0'
  }

  if('NS' %in% names(direction_table)) {
    ns_label = sprintf('NS: %s', direction_table['NS'])
  } else {
    ns_label = 'NS: 0'
  }

  de_list$direction_count = factor(
    de_list$direction,
    levels = c('Up', 'Down', 'NS'),
    labels = c(up_label, down_label, ns_label))

  if(all( !( c('Up', 'Downl') %in% names(direction_table) ) )) {
    warning(sprintf('No genes were DE at fdr < %s and |log2fc| > %s. Consider a different threshold.', fdr_cutoff, logfc_cutoff))
  }

  # Transform qval to -log10 scale
  de_list$log10qval = -log10(de_list[, padj])

  # Add top 10 Up and 10 Down gene labels
  # de_list is assumed to be ordered by Call/diff_exp and then qvalue from deseq2_diffex.R and ballgown_diffex.R
  top = rbind(
    head(subset(de_list, direction == 'Up'), 10),
    head(subset(de_list, direction == 'Down'), 10))
  top$label = top$external_gene_name
  de_list = merge(x = de_list, y = top[, c(id,'label')], by = id, all.x = TRUE, sort = FALSE)

  # Volcano Plot
  volcano_plot = ggplot(de_list, aes_string(x = log2fc, y = 'log10qval', color = 'direction_count', alpha = 0.5)) +
    scale_alpha(guide = 'none') +
    geom_point(size = 1) +
    scale_color_manual(name = '', values=c('#B31B21', '#1465AC', 'darkgray')) +
    geom_vline(xintercept = c(0, -1*logfc_cutoff, logfc_cutoff), linetype = c(1, 2, 2), color = c('black', 'black', 'black')) +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2, color = 'black') +
    labs(
      title = sprintf('%s_v_%s', exp_name, con_name),
      x = 'Log2 fold-change',
      y = '-Log10 adjusted q-value') +
    theme_classic()
  # Add gene symbol labels
  if(!all(is.na(de_list$label))) {
    volcano_plot = volcano_plot + ggrepel::geom_label_repel(label = de_list$label, force = 3, segment.alpha = 0.4)
  }
  ggsave(filename = out_filepath_pdf, plot = volcano_plot, height = 8, width = 8, dpi = 300)
  ggsave(filename = out_filepath_png, plot = volcano_plot, height = 8, width = 8, dpi = 300)

  return(volcano_plot)
}

####################################################
# Set up for calling differentially expressed genes

model_name = snakemake@wildcards[['model_name']]
factor_name = snakemake@config[['diffex']][[model_name]][['DESeq2']][['factor_name']]

contrast_name = snakemake@wildcards[['contrast']]
base_filename = contrast_name #Assign base filename using contrast wildcard
cont_split = unlist(strsplit(contrast_name, "_v_"))
#Define the test and reference name from this
test_name = cont_split[1] ; reference_name = cont_split[2]

# Establish cutoffs
fdr_cutoff = as.numeric(snakemake@config[['diffex']][[model_name]][['adjustedPValue']])
fc_cutoff = log2(as.numeric(snakemake@config[['diffex']][[model_name]][['linear_fold_change']]))

# Set up multithreading
multicore_param = MulticoreParam(workers = snakemake@resources[['cpus']])
register(multicore_param, default=TRUE)

# Get phenotype matrix
sample.info.file = snakemake@config[['samplesheet']]
pdata = read.csv(sample.info.file, comment.char = "#")
pdata$input_dir = NULL
pdata$input_glob = NULL # Drop either if they exist - transition from dir to glob

message(sprintf('Testing %s: %s vs %s', factor_name, test_name, reference_name))

# Load DESeq dataset, generated via deseq2_init into variable dds
load(snakemake@input[['rda']])


results.params = snakemake@config[['diffex']][[model_name]][['DESeq2']][['results']]
# Parse the params for results {DESeq2} call, if it can't be converted return it as string
results.params.parsed = lapply(results.params, function(x) {
  tryCatch(eval(parse(text=x)),
  error=function(e){
    message(paste0("DESeq2 config parameter '", x, "' not parsable to known R type. Leaving as string"))
    x=as.character(x)
  })
})

# Add dds object to params list
results.params.parsed[['object']] = dds
# Set parallel to TRUE
results.params.parsed[['parallel']] = TRUE

message('Pulling results from DESeq2 object')
# Grab the results
res = do.call(results, results.params.parsed)

# Order by adjusted p value
res = res[order(res$padj),]

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
  file=snakemake@output[['gene_list']],
  append = FALSE, sep = '\t', na = '.', row.names = FALSE, quote = FALSE)

##################
# Add annotations
annotation_tsv = snakemake@config[['references']][['annotation_tsv']]
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

# Write the annotated gene list - plain text
message('Writing annotated diffex genes')
write.table(
  x = annotated_results,
  file=snakemake@output[['annot_results']],
  append = FALSE, sep = '\t', na = '.', row.names = FALSE, quote = FALSE)

# Write the annotated gene list - xlsx
message('Writing annotated diffex genes xlsx')
# Load glossary for insertion into the workbook
glossary = read.delim(snakemake@input[['glossary']], stringsAsFactors = FALSE)
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
saveWorkbook(diffex_wb, snakemake@output[['annot_results_xlsx']], overwrite = TRUE)

######################
# Create volcano plot
message(sprintf('Plotting volcano plot for %s %s', factor_name, contrast_name))

out_pdf = snakemake@output[['volcano_plot_pdf']]
out_png = snakemake@output[['volcano_plot_png']]

volcano_plot = plot_volcano(
  de_list = annotated_results,
  method = "deseq2",
  exp_name = test_name,
  con_name = reference_name,
  fdr_cutoff = fdr_cutoff,
  logfc_cutoff = fc_cutoff,
  out_filepath_pdf = out_pdf,
  out_filepath_png = out_png)
