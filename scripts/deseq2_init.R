##########
# Set up logging and save snakemake S4 object (for debugging or running manually)
log = file(snakemake@log[[1]], open = 'wt')
sink(log)
sink(log, type = 'message')
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
model_name = snakemake@wildcards[['model_name']]
design = snakemake@config[['diffex']][[model_name]][['DESeq2']][['design']]
deseq2.params = snakemake@config[['diffex']][[model_name]][['DESeq2']][['DESeq2']]
# Get phenotype matrix
sample.info.file = snakemake@config[['samplesheet']]
pdata = read.csv(sample.info.file, comment.char = "#")

# Load count data (imported in previous rule via tximport to variable txi.rsem.gene.results)
load(snakemake@input[['data_import']])
# Create DESeqDataSet and filter lowly expressed genes
message('Initializing DESeq2 result')
if(exists('txi.rsem.gene.results')){
    # Create dataset from tximport object if it was imported above
    dds = DESeqDataSetFromTximport(txi = txi.rsem.gene.results, colData = pdata, design = as.formula(design))
} else if(exists('counts')){
    # Create dataset from count matrix if it was imported above
    dds = DESeqDataSetFromMatrix(countData = counts, colData = pdata, design = as.formula(design))
}

# Filter lowly expressed genes
threshold = snakemake@config[['diffex']][['count_min_cutoff']]
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
multicore_param = MulticoreParam(workers = snakemake@threads)
register(multicore_param, default=TRUE)

# Call to DESeq
dds = do.call(DESeq, deseq2.params.parsed)

save(dds, file=snakemake@output[['rda']])
