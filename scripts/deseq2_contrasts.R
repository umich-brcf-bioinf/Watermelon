##########
# Set up logging and save snakemake S4 object (for debugging or running manually)
log = file(snakemake@log[[1]], open='wt')
sink(log, split=TRUE)
save(snakemake, file = snakemake@params[['snakemake_rdata']])

#load("/nfs/med-bfx-activeprojects/trsaari/example_output_Watermelon/diffex_results/deseq2/gene_lists/pheno_gender/DM.fem_v_NDM.fem_snakemake.rda") #TWS DEBUG

##########
# Load libraries

lib.vector = c("BiocParallel", "data.table", "DESeq2")
foo = suppressMessages(lapply(lib.vector, library, character.only=T, warn.conflicts=F, quietly=T))

##########
# Main

factor_name = snakemake@wildcards[['factor_name']]

#Match contrasts wildcard with contrast from comparison config
#all of this is just to ensure that the correct factorName is used, which can only be pulled from the config
#Wildcard portion
contrast.wc = snakemake@wildcards[['contrast']]
base.file.name = contrast.wc #Assign base filename using contrast wildcard
cont.wc.split = unlist(strsplit(contrast.wc, "_v_"))
cont.wc.rejoined = paste(cont.wc.split[[1]], cont.wc.split[[2]], sep="^")
#Config portion
contrast.conf = snakemake@config[['diffex']][[factor_name]][['DESeq2']][['results']][['contrasts']]
cont.match = contrast.conf[which(grepl(cont.wc.rejoined, contrast.conf, fixed=T))] #Find the match
#Define the variables from this
conparts = unlist(strsplit(cont.match, "^", fixed=T))
factorName = conparts[1] ; testName = conparts[2] ; referenceName = conparts[3]

# Establish cutoffs
fdr_cutoff = as.numeric(snakemake@config[['diffex']][['adjustedPValue']])
fc_cutoff = log2(as.numeric(snakemake@config[['diffex']][['fold_change']]))

# Set up multithreading
multicore_param = MulticoreParam(workers = snakemake@threads)
register(multicore_param, default=TRUE)

# Get phenotype matrix
sample.info.file = snakemake@config[['sample_description_file']]
pdata = read.csv(sample.info.file)

message(sprintf('Testing %s: %s vs %s', factorName, testName, referenceName))

# Load DESeq dataset, generated via deseq2_init into variable dds
load(snakemake@input[['rda']])

# Grab the results
res = results(
  dds,
  contrast = list(paste0(factorName, testName), paste0(factorName, referenceName)),
  parallel = T
)
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

# Add the individual diffexData to the collection
# de_results[[factorName]][[baseName]] = list(gene_list = diffexData)

# TWS - Will need to see if/how subsequent steps depend on de_results structure
# Will have to decide if it makes sense to just read in the table files or to load in the different rdata files
# They're going to be separate with the new paradigm anyways.

#write to individual tab-delimited txt files
message('Writing diffex genes')
write.table(
  x = diffexData,
  file=snakemake@output[['gene_list']],
  append = FALSE, sep = '\t', na = 'NA', row.names = FALSE, quote = FALSE)
