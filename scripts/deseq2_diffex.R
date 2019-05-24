library(optparse)

option_list = list(
    make_option('--count_file', type='character', help='[Required] Path to stringtie counts matrix'),
    make_option('--config_file', type='character', help='[Required] Path to configfile'),
    make_option('--threads', type = 'integer', default = 4, help='Number of threads to use.')
)
opt = parse_args(OptionParser(option_list=option_list))

#count_file = opt$count_file
#config_file = opt$config_file
#threads = opt$threads

results_paths_file = "example_output_Watermelon/gene_results_files.txt" #TWS DEBUG
config_file = "Watermelon/config/config_190521_test.yaml" #TWS DEBUG
threads = 4 #TWS DEBUG

#######################################

library(BiocParallel)
library(data.table)
library(tximport)
library(DESeq2)
library(readr)
library(stringr)
library(yaml)

#######################################

multicore_param = MulticoreParam(workers = threads)
register(multicore_param, default=TRUE)

#######################################

versus_split = function(str) {
    tmp = unlist(strsplit(str, '_v_', fixed = TRUE))
    tmp = str_trim(tmp)
    return(tmp)
}

handle_comparisons = function(yaml) {
    comps = lapply(names(yaml$comparisons),
        function(name){
            comparison = yaml$comparisons[[name]]
            tmp = data.frame(t(sapply(comparison, versus_split, USE.NAMES = FALSE)))
            colnames(tmp) = c('exp', 'con')
            tmp$out_name = apply(tmp, 1, paste, collapse = '_v_')
            tmp$factor = name
            return(tmp)
        }
    )
    comps = Reduce(rbind, comps)
    return(comps)
}

#######################################
# yaml parsing

yaml = read_yaml(config_file)
diffex_yaml = read_yaml("Watermelon/config/example_comparisons.yaml") #TWS DEBUG

# Establish directories
diffex_dir = yaml$dirs$diffex_output
results_dir = sprintf('%s/deseq2/01-deseq2_diffex', diffex_dir)
counts_dir = sprintf('%s/counts', results_dir)
gene_lists_dir = sprintf('%s/gene_lists', results_dir)

# Get phenotyp matrix
pdata = read_csv("Watermelon/config/example_sample_description_subset.csv") #TWS DEBUG

# Establish comparisons
comparisons = handle_comparisons(diffex_yaml)
comparison_types = names(comparisons)

# Establish cutoffs
fdr_cutoff = diffex_yaml$deseq2_adjustedPValue
fc_cutoff = log2(diffex_yaml$fold_change)

#######################################

# Import count data
files.table <- read.table(results_paths_file, sep = '\t', col.names = c("sample", "filepath"))
files.list <- file.path(files.table$filepath)
names(files.list) <- files$sample

txi.rsem.gene.results <- tximport(files.list, type = "rsem", txIn = F, txOut = F)
#Some genes have length zero (what does this even mean?), causing issues with creating DESeqDataSet
#Mike Love recommends changing these from 0 to 1; https://support.bioconductor.org/p/84304/#84368
txi.rsem.gene.results$length[txi.rsem.gene.results$length == 0] <- 1

#######################################
# Create DESeqDataSet from matrix and filter lowly expressed genes

message('Initializing DESeq2 result')
dds = DESeqDataSetFromTximport(txi = txi.rsem.gene.results, colData = pdata, design = ~0 + pheno) #TWS DEBUG - combinatoric_group doesn't exist
# QUESTION: Should this be more stringent?
dds = dds[ rowSums(counts(dds)) > 1, ]

# QUESTION: Why do this here before DESeq()?
rld = rlog(dds, blind = FALSE)

message('Calculating dispersion')
#Normalize and calculate dispersions
dds = DESeq(dds, betaPrior = TRUE, parallel = TRUE)

###################
# Extract counts
raw_counts = counts(dds) #raw
norm_counts = counts(dds, normalized = TRUE) #raw/lib size factors

###################
# Write counts
message('Writing count files')
raw_counts_df = as.data.frame(raw_counts)
data.table::setDT(raw_counts_df, keep.rownames = TRUE)[] #set rownames to valid column
setnames(raw_counts_df, 'rn', 'id')
write.table(x = raw_counts_df, file=paste(countsDir,'raw_counts.txt', sep = '/'),
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

norm_counts_df = as.data.frame(norm_counts)
data.table::setDT(norm_counts_df, keep.rownames = TRUE)[] #set rownames to valid column
setnames(norm_counts_df, 'rn', 'id')
write.table(x = norm_counts_df, file=paste(countsDir,'depth_normalized_counts.txt', sep = '/'),
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

rld_df = as.data.frame(assay(rld))
data.table::setDT(rld_df, keep.rownames = TRUE)[] #set rownames to valid column
setnames(rld_df, 'rn', 'id')
write.table(x = rld_df, file=paste(countsDir,'rlog_normalized_counts.txt', sep = '/'),
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

#######################################
# Diffex analysis in loop

contrastData = read.table(file = "example_output_Watermelon/contrasts.txt", header=TRUE, sep = '\t', stringsAsFactors = FALSE, strip.white = TRUE) #TWS DEBUG
colData <- read.table(file = "example_output_Watermelon/sample_metadata.txt" , header = TRUE, sep = '\t', stringsAsFactors = FALSE, na.strings = '', check.names = FALSE) #TWS DEBUG

# Container for the results
de_results = split(contrastData, contrastData$factor)
de_results = lapply(de_results, function(de){split(de, de$base_file_name)})

for (i in 1:nrow(contrastData)){

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
