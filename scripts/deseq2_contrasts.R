##########
# Set up logging and save snakemake S4 object (for debugging or running manually)
log = file(snakemake@log[[1]], open='wt')
sink(log, split=TRUE)
save(snakemake, file = snakemake@params[['snakemake_rdata']])

#load("/nfs/med-bfx-activeprojects/trsaari/example_output_Watermelon/diffex_results/pheno_gender/DESeq2/pheno.Gend_DM.fem_DM.male_snakemake.rda") #TWS DEBUG

##########
# Load libraries
library(BiocParallel)
library(data.table)
library(DESeq2)

##########
# Main

model_name = snakemake@wildcards[['model_name']]
contrast = snakemake@wildcards[['contrast']]

# Establish cutoffs
fdr_cutoff = as.numeric(snakemake@config[['diffex']][['adjustedPValue']])
fc_cutoff = log2(as.numeric(snakemake@config[['diffex']][['fold_change']]))

# Set up multithreading
multicore_param = MulticoreParam(workers = snakemake@threads)
register(multicore_param, default=TRUE)


# Get phenotype matrix
sample.info.file = snakemake@config[['sample_description_file']]
pdata = read.csv(sample.info.file)

# Parse contrast string
conparts = unlist(strsplit(contrast, "_", fixed=T))
factorName = conparts[1] ; testName = conparts[2] ; referenceName = conparts[3]
# Give _v_ style basename

message(sprintf('Testing %s: %s vs %s', factorName, testName, referenceName))

# Collect test samples in a list
testSamples = subset(x=pdata, subset=pdata[, factorName] == testName)
# Collect reference samples in a list
referenceSamples = subset(x=pdata, subset=pdata[, factorName] == referenceName)

# TWS - The above is imitating the previous approach, but are we really only interested in the # samples in each (for listValues?)
# resultsNames of dds are "pheno.GendDM.fem"   "pheno.GendDM.male"  "pheno.GendNDM.fem"  "pheno.GendNDM.male"
# So I'll set the contrast up as below, i.e. using paste0

# Load DESeq dataset (dds)
load(snakemake@input[['rda']])

# Grab the results
res = results(
  dds,
  contrast = list(paste0(factorName, testName), paste0(factorName, referenceName)),
  parallel = T,
  listValues = c(1/nrow(testSamples), -1/nrow(referenceSamples))
)
res = res[order(res$padj),]

diffexData = as.data.frame(res)
#set rownames to valid column
data.table::setDT(diffexData, keep.rownames = TRUE)[] #TWS - Is data.table really useful here at all?

# TWS LEFT OFF HERE

#rename first column
colnames(diffexData) = c('id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')
diffexData$Condition = testName
diffexData$Control = referenceName

#make DEG calls and select DEGs
diffexData$Call = rep('NO', nrow(diffexData))
diffexData$Call[which(diffexData$padj <= pval & abs(diffexData$log2FoldChange) >= log2(fc))] = 'YES'
diffexData = diffexData[order(-rank(diffexData$Call), diffexData$pvalue), ]

# Add the individual diffexData to the collection
de_results[[factorName]][[baseName]] = list(gene_list = diffexData)

#write to individual tab-delimited txt files
message('Writing diffex genes')
write.table(
  x = diffexData,
  file=paste0(dir_diffex,'/',baseName, '.txt'),
  append = FALSE, sep = '\t', na = 'NA', row.names = FALSE, quote = FALSE)
}
