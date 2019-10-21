##########
# Set up logging and save snakemake S4 object (for debugging or running manually)
log = file(snakemake@log[[1]], open='wt')
sink(log, split=TRUE)
save(snakemake, file = snakemake@params[['snakemake_rdata']])

#load("/nfs/med-bfx-activeprojects/trsaari/sandbox/20190725_test_Delono_RS1/analysis_test_Delano_RS1/diffex_results/deseq2/gene_lists/phenotype.CellState.treatment/DIO.WCLP.none_v_DIO.DCLP.none_snakemake.rda") #TWS DEBUG

#Isolate conda environment: https://github.com/conda-forge/r-base-feedstock/issues/37
#If we move away from conda in the future, we may want to remove this
.libPaths(R.home("library"))

##########
# Load libraries

lib.vector = c("BiocParallel", "data.table", "DESeq2")
foo = suppressMessages(lapply(lib.vector, library, character.only=T, warn.conflicts=F, quietly=T))

##########
# Main

model_name = snakemake@wildcards[['model_name']]
factorName = snakemake@config[['diffex']][[model_name]][['DESeq2']][['factor_name']]

contrast = snakemake@wildcards[['contrast']]
base.file.name = contrast #Assign base filename using contrast wildcard
cont.split = unlist(strsplit(contrast, "_v_"))
#Define the test and reference name from this
testName = cont.split[1] ; referenceName = cont.split[2]

# Establish cutoffs
fdr_cutoff = as.numeric(snakemake@config[['diffex']][['adjustedPValue']])
fc_cutoff = log2(as.numeric(snakemake@config[['diffex']][['linear_fold_change']]))

# Set up multithreading
multicore_param = MulticoreParam(workers = snakemake@threads)
register(multicore_param, default=TRUE)

# Get phenotype matrix
sample.info.file = snakemake@config[['sample_description_file']]
pdata = read.csv(sample.info.file, comment.char = "#")

message(sprintf('Testing %s: %s vs %s', factorName, testName, referenceName))

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

# Grab the results
res = do.call(results, results.params.parsed)

# Order by adjusted p value
res = res[order(res$padj),]

diffexData = as.data.frame(res)
#set rownames to valid column
data.table::setDT(diffexData, keep.rownames = TRUE)[] #TWS - Is data.table really useful here at all?

#rename first column
colnames(diffexData) = c('id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')
diffexData$Condition = testName
diffexData$Control = referenceName

#make DEG calls and select DEGs
diffexData$Call = rep('NO', nrow(diffexData))
diffexData$Call[which(diffexData$padj <= fdr_cutoff & abs(diffexData$log2FoldChange) >= fc_cutoff)] = 'YES'
diffexData = diffexData[order(-rank(diffexData$Call), diffexData$pvalue), ]


#write to individual tab-delimited txt files
message('Writing diffex genes')
write.table(
  x = diffexData,
  file=snakemake@output[['gene_list']],
  append = FALSE, sep = '\t', na = 'NA', row.names = FALSE, quote = FALSE)
