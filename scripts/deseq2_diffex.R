log = file(snakemake@log[[1]], open='wt')
sink(log)
sink(log, type='message')
save(snakemake, file = snakemake@params[['snakemake_rdata']])

#load("/nfs/med-bfx-activeprojects/trsaari/example_output_Watermelon/diffex_results/pheno_gender/DESeq2/pheno.Gend_DM.male_NDM.male_snakemake.rda") #TWS DEBUG

# library(optparse)
#
# option_list = list(
#     make_option('--rsem_dir', type='character', help='[Required] Path to dir containing rsem results'),
#     make_option('--config_file', type='character', help='[Required] Path to configfile'),
#     make_option('--threads', type = 'integer', default = 4, help='Number of threads to use.')
# )
# opt = parse_args(OptionParser(option_list=option_list))

rsem_dir = snakemake@params[['rsem_dir']]
#rsem_dir = paste0("/nfs/med-bfx-activeprojects/trsaari/example_output_Watermelon/", rsem_dir) #TWS DEBUG
#config_file = opt$config_file
threads = snakemake@threads

#######################################

library(BiocParallel)
library(data.table)
library(tximport)
library(DESeq2)
#library(readr)
#library(stringr)
#library(yaml)

#######################################

multicore_param = MulticoreParam(workers = threads)
register(multicore_param, default=TRUE)

#######################################

#Function to generate a named list of filepaths
#For loading these data via tximport
named.filepaths.from.dir <- function(rsem.dir, file.extension) {
  fnames <- list.files(rsem.dir)
  genes.fnames <- fnames[grepl(file.extension, fnames)]
  sample.names <- sub(file.extension, "", genes.fnames, fixed=T)
  fpathlist <- file.path(rsem.dir, genes.fnames)
  names(fpathlist) <- sample.names
  return(fpathlist)
}

#Function to generate DESeq expression from diffex config parameters
# DESeq.call.from.params = function(obj, paramslist) {
#   param.kv = list()
#   for(i in seq_along(paramslist)) { append(param.kv, paste(names(paramslist)[i], paramslist[[i]], sep="=")) }
#   return(paste("DESeq", "obj"
# }
#test.call <- parse(text="DESeq(object=dds, full=~0 + pheno, reduced=~0, betaPrior=False, test='LRT')")
#foo <- eval(test.call)

# #######################################
# # yaml parsing
#
# yaml = read_yaml(config_file)
# diffex_yaml = read_yaml("/nfs/med-bfx-activeprojects/trsaari/Watermelon/config/example_comparisons.yaml") #TWS DEBUG
#
# # Establish directories
# diffex_dir = yaml$dirs$diffex_output
# results_dir = sprintf('%s/deseq2/01-deseq2_diffex', diffex_dir)
# counts_dir = sprintf('%s/counts', results_dir)
# gene_lists_dir = sprintf('%s/gene_lists', results_dir)

# Get comparison info
factor_name = snakemake@wildcards[['factor_name']]
deseq2.params = snakemake@config[['diffex']][[factor_name]][['DESeq2']][['DESeq2']]

# Get phenotype matrix
sample.info.file = snakemake@config[['sample_description_file']]
pdata = read.csv(sample.info.file)

# # Establish comparisons
# comparisons = handle_comparisons(diffex_yaml)
# comparison_types = names(comparisons)

# Establish cutoffs
fdr_cutoff = as.numeric(snakemake@config[['diffex']][['adjustedPValue']])
fc_cutoff = log2(as.numeric(snakemake@config[['diffex']][['linear_fold_change']]))

#######################################

# Import count data
gene.files.list <- named.filepaths.from.dir(rsem_dir, ".genes.results")
txi.rsem.gene.results <- tximport(gene.files.list, type = "rsem", txIn = F, txOut = F)
#Some genes have length zero (what does this even mean?), causing issues with creating DESeqDataSet
#Mike Love recommends changing these from 0 to 1; https://support.bioconductor.org/p/84304/#84368
txi.rsem.gene.results$length[txi.rsem.gene.results$length == 0] <- 1

#######################################
# Create DESeqDataSet from matrix and filter lowly expressed genes

message('Initializing DESeq2 result')
#Import data from rsem
dds = DESeqDataSetFromTximport(txi = txi.rsem.gene.results, colData = pdata, design = as.formula(deseq2.params$full))
# QUESTION: Should this be more stringent?
# TWS
dds = dds[ rowSums(counts(dds)) > 1, ]

# QUESTION: Why do this here before DESeq()?
rld = rlog(dds, blind = FALSE)

#Parse the params for DESeq call, if it can't be converted return it as string
deseq2.params.parsed = lapply(deseq2.params, function(x) {
  tryCatch(eval(parse(text=x)),
  error=function(e){
    message(paste0(x, " cannot be parsed. Leaving as string"))
    x=as.character(x)
  })
})
#Add dds object to params list
deseq2.params.parsed[['object']] = dds

# TWS - betaPrior cannot be used with a model that has an intercept. Don't think this should be generalized...
# message('Calculating dispersion')
# #Normalize and calculate dispersions
# dds = DESeq(dds, betaPrior = F, parallel = TRUE)


#Call to DESeq
dds = do.call(DESeq, deseq2.params.parsed)

quit()

# TWS - This should be moved out to a separate rule/script which will be executed once for all included samples
# ###################
# # Extract counts
# raw_counts = counts(dds) #raw
# # TWS - need to estimate size factors before this can be done
# dds.withSF = estimateSizeFactors(dds)
# norm_counts = counts(dds.withSF, normalized = TRUE) #raw/lib size factors
#
# ###################
# # Write counts
# message('Writing count files')
# raw_counts_df = as.data.frame(raw_counts)
# data.table::setDT(raw_counts_df, keep.rownames = TRUE)[] #set rownames to valid column
# setnames(raw_counts_df, 'rn', 'id')
# write.table(x = raw_counts_df, file=paste(countsDir,'raw_counts.txt', sep = '/'),
#     append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
#
# norm_counts_df = as.data.frame(norm_counts)
# data.table::setDT(norm_counts_df, keep.rownames = TRUE)[] #set rownames to valid column
# setnames(norm_counts_df, 'rn', 'id')
# write.table(x = norm_counts_df, file=paste(countsDir,'depth_normalized_counts.txt', sep = '/'),
#     append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
#
# rld_df = as.data.frame(assay(rld))
# data.table::setDT(rld_df, keep.rownames = TRUE)[] #set rownames to valid column
# setnames(rld_df, 'rn', 'id')
# write.table(x = rld_df, file=paste(countsDir,'rlog_normalized_counts.txt', sep = '/'),
#     append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

#######################################
# Diffex analysis in loop

contrastData = read.table(file = "/nfs/med-bfx-activeprojects/trsaari/example_output_Watermelon/contrasts.txt", header=TRUE, sep = '\t', stringsAsFactors = FALSE, strip.white = TRUE) #TWS DEBUG
colData <- read.table(file = "/nfs/med-bfx-activeprojects/trsaari/example_output_Watermelon/sample_metadata.txt" , header = TRUE, sep = '\t', stringsAsFactors = FALSE, na.strings = '', check.names = FALSE) #TWS DEBUG

# Container for the results
de_results = split(contrastData, contrastData$factor)
de_results = lapply(de_results, function(de){split(de, de$base_file_name)})

for (i in 1:nrow(contrastData)){
i = 7 #TWS DEBUG
    #dir assignment
    dir_diffex = paste(diffexDir, contrastData$factor[i], sep = '/')

    #collect references, etc
    factorName = contrastData$factor[i]
    testName = contrastData$test_level[i]
    referenceName = contrastData$reference_level[i]
    baseName = contrastData$base_file_name[i]

    message(sprintf('Testing %s: %s vs %s', factorName, testName, referenceName))
    #diffex analysis use generic name, subset DEGs

    #collect test samples in a list, from groups column
    testSamples = subset(
        x = colData,
        subset = colData[, factorName] == testName & !duplicated(colData$combinatoric_group),
        select = combinatoric_group)
    # put 'combinatoric_group' at beginning of name
    testSamples$combinatoric_group = sub('^', 'combinatoric_group', testSamples$combinatoric_group)
    # replace '^' with '.' to match name in deseq obj
    testSamples$combinatoric_group = gsub('^', '.', testSamples$combinatoric_group, fixed = TRUE)
    # replace '-' with '.' to match name in deseq obj
    testSamples$combinatoric_group = gsub('-', '.', testSamples$combinatoric_group, fixed = TRUE)


    #collect reference samples in a list
    referenceSamples = subset(
        x = colData,
        subset = colData[,factorName] == referenceName & !duplicated(colData$combinatoric_group),
        select = combinatoric_group)
    # put 'combinatoric_group' at beginning of name
    referenceSamples$combinatoric_group = sub('^', 'combinatoric_group', referenceSamples$combinatoric_group)
    # replace '^' with '.' to match name in deseq obj
    referenceSamples$combinatoric_group = gsub('^', '.', referenceSamples$combinatoric_group, fixed = TRUE)
    # replace '-' with '.' to match name in deseq obj
    referenceSamples$combinatoric_group = gsub('-', '.', referenceSamples$combinatoric_group, fixed = TRUE)


    # Grab the results
    res = results(
        dds,
        parallel = TRUE,
        contrast = list(c(unlist(testSamples)), c(unlist(referenceSamples))),
        listValues = c(1/nrow(testSamples), -1/nrow(referenceSamples)))
    res = res[order(res$padj),]


    diffexData = as.data.frame(res)
    #set rownames to valid column
    data.table::setDT(diffexData, keep.rownames = TRUE)[]
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

save_list = c('dds', 'rld', 'raw_counts', 'norm_counts', 'de_results')
deseq2_rda_file = sprintf('%s/deseq2_data.rda', out_dir)
message(sprintf('Saving deseq RData to %s', deseq2_rda_file))
save(list = save_list, file = deseq2_rda_file)

cat('deseq2_diffex.R done.\n')
