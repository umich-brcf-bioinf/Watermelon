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
suppressMessages(library(BiocParallel, warn.conflicts=F, quietly=T))
suppressMessages(library(DESeq2, warn.conflicts=F, quietly=T))

##########
# Main

# Get deseq params from config
factor_name = snakemake@wildcards[['factor_name']]
deseq2.params = snakemake@config[['diffex']][[factor_name]][['DESeq2']][['DESeq2']]
# Get phenotype matrix
sample.info.file = snakemake@config[['sample_description_file']]
pdata = read.csv(sample.info.file, comment.char = "#")

# Load count data (imported in previous rule via tximport to variable txi.rsem.gene.results)
load(snakemake@input[['data_import']])
# Create DESeqDataSet and filter lowly expressed genes
message('Initializing DESeq2 result')
# Create dataset from tximport object
dds = DESeqDataSetFromTximport(txi = txi.rsem.gene.results, colData = pdata, design = as.formula(deseq2.params$full))
# Filter lowly expressed genes - QUESTION: Should this be more stringent?
dds = dds[ rowSums(counts(dds)) > 1, ]

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
multicore_param = MulticoreParam(workers = snakemake@threads)
register(multicore_param, default=TRUE)

# Call to DESeq
dds = do.call(DESeq, deseq2.params.parsed)

save(dds, file=snakemake@output[['rda']])
