##########
# Set up logging and save snakemake S4 object (for debugging or running manually)
log = file(snakemake@log[[1]], open='wt')
sink(log, split=TRUE)
save(snakemake, file = snakemake@params[['snakemake_rdata']])

#Isolate conda environment: https://github.com/conda-forge/r-base-feedstock/issues/37
#If we move away from conda in the future, we may want to remove this
.libPaths(R.home("library"))

##########
# Load libraries
suppressMessages(library(data.table, warn.conflicts=F, quietly=T))
suppressMessages(library(tximport, warn.conflicts=F, quietly=T))
suppressMessages(library(DESeq2, warn.conflicts=F, quietly=T))

##########
# Define functions
named.filepaths.from.dir <- function(rsem.dir, file.extension) {
  fnames <- list.files(rsem.dir)
  genes.fnames <- fnames[grepl(file.extension, fnames)]
  sample.names <- sub(file.extension, "", genes.fnames, fixed=T)
  fpathlist <- file.path(rsem.dir, genes.fnames)
  names(fpathlist) <- sample.names
  return(fpathlist)
}

##########
# Main

# Get phenotype matrix
sample.info.file = snakemake@config[['sample_description_file']]
pdata = read.csv(sample.info.file, comment.char = "#")

# Import data
message('Importing rsem data')
rsem_dir = snakemake@params[['rsem_dir']]
gene.files.list <- named.filepaths.from.dir(rsem_dir, ".genes.results")
txi.rsem.gene.results <- tximport(gene.files.list, type = "rsem", txIn = F, txOut = F)
#Some genes have length zero (what does this even mean?), causing issues with creating DESeqDataSet
#Mike Love recommends changing these from 0 to 1; https://support.bioconductor.org/p/84304/#84368
txi.rsem.gene.results$length[txi.rsem.gene.results$length == 0] <- 1

# Create DESeqDataSet and filter lowly expressed genes
# Import data from rsem
# Note use of no design (i.e. = ~ 1) here since we only care about counts. For diffex, will need to use actual design
dds = DESeqDataSetFromTximport(txi = txi.rsem.gene.results, colData = pdata, design = ~ 1)
# Filter lowly expressed genes QUESTION: Should this be more stringent?
dds = dds[ rowSums(counts(dds)) > 1, ]

# Extract counts
raw_counts = counts(dds) #raw
# TWS - need to estimate size factors before normalization can be done
dds.withSF = estimateSizeFactors(dds)
norm_counts = counts(dds.withSF, normalized = TRUE) #raw/lib size factors
# Log normed
rld = rlog(dds, blind = FALSE)

#Save the dataset and extracted counts to rdata files
save(txi.rsem.gene.results, file=snakemake@output[['txi']])
counts.tables.list = c('raw_counts', 'norm_counts', 'rld')
save(list = counts.tables.list, file=snakemake@output[['count_tables_rda']])


# Write counts tables
message('Writing count files')
raw_counts_df = as.data.frame(raw_counts)
data.table::setDT(raw_counts_df, keep.rownames = TRUE)[] #set rownames to valid column
setnames(raw_counts_df, 'rn', 'id')
write.table(x = raw_counts_df, file=snakemake@output[['raw']],
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

norm_counts_df = as.data.frame(norm_counts)
data.table::setDT(norm_counts_df, keep.rownames = TRUE)[] #set rownames to valid column
setnames(norm_counts_df, 'rn', 'id')
write.table(x = norm_counts_df, file=snakemake@output[['norm']],
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

rld_df = as.data.frame(assay(rld))
data.table::setDT(rld_df, keep.rownames = TRUE)[] #set rownames to valid column
setnames(rld_df, 'rn', 'id')
write.table(x = rld_df, file=snakemake@output[['rlog']],
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
