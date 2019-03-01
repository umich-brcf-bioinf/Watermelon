#!/usr/bin/env Rscript

#######
# Perform differential expression analysis using DESeq2
# Can be single or multi-factor, and multi-level (will use combined factor approach)
#
# Required Input files or options:
#  1) Raw count matrix
#     - tab-delimited
#     - columns are samples
#     - rows are genes
#     - no extra spaces at the end of any row
#  2) Metadata file
#     - tab-delimited between 1st and 2nd column, and ^ delimited between all other columns (factors)
#     - Name column followed by 1 column for each factor (ie, Name, Gender, Diet, Cell)
#     - Header row followed by 1 row for each sample (ie, Name, Sample_1, Sample_2, Sample_3)
#     - Row and column together give matrix with sample names and the specific level for each factor for that specific sample
#  3) Contrast file
#     - tab-delimited
#     - 1 row per contrast to be made
#     - Columns in order:
#       - factor = the factor (variable) of interest for that contrast
#       - test_level = the level of the factor to be considered the numerator, or test level
#       - reference_level = the level of the factor to be considered the denominator, or reference level
#       - file_name = the output file name for that contrast
#  4) Fold-change cutoff for consideration as differentially expressed
#  5) Adjusted p-value cutoff for consideration as differentially expressed
#
#  ALL OPTIONS ARE REQUIRED. NO DEFAULTS ARE SET.
#
# Output files:
#  1) Tab-delimited text files for each comparison (between levels) for this single factor.
#     - Sorted by adjusted pvalue
#     - A file containing the full differential expression dataset is created
#     - A file containing only data for differentially expressed genes is created
#  2) Excel file for each comparison (as above) names by the comparison.
#     - Sorted by adjusted pvalue
#     - A sheet containing the full differential expression dataset is created
#     - A sheet containing only data for differentially expressed genes is created
#######

DEFAULT_THREADS = 4

print_options <- function(opt) {
  opt_df <- data.frame(matrix(unlist(opt), byrow=T),stringsAsFactors=FALSE)
  row.names(opt_df)<-names(opt)
  names(opt_df) <- NULL
  message('options:')
  opt_df
}

library(optparse)

#create parser object
option_list = list(
  make_option(c('-c', '--countDataFile'), type = 'character', default = NULL,
              help = 'Name of file containing tab-delimited RAW counts.'),
  make_option(c('-m', '--metaDataFile'), type = 'character', default = NULL,
              help = 'Name of file containing tab-delimited sample information, ie. factors, levels, etc.'),
  make_option(c('-f', '--contrastFile'), type = 'character', default = NULL,
              help = 'Factor to compare from metaDataFile, must be a column name from metaDataFile'),
  make_option(c('-o', '--outDir'), type = 'character', default = 'deseq_output',
              help = 'Output directory to write all generated subdirectories.'),
  make_option(c('--foldChange'), type = 'character', default = '1.5',
              help = 'Absolute numeric value above which a gene is considered for differential expression.'),
  make_option(c('--adjustedPValue'), type = 'character', default = '0.05',
              help = 'Adjusted p-value below which a genes is considered for differential expression.'),
  make_option(c('--threads'), type = 'integer', default = DEFAULT_THREADS,
              help = 'Number of parallel processes (CPUs) to use.'),
  make_option(c('--countsDir'), type = 'character', default = 'counts',
              help = 'Output directory for count data.'),
  make_option(c('--geneListsDir'), type = 'character', default = 'gene_lists',
              help = 'Output directory for differential expression data.')
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#check for necessary files
if (is.null(opt$countDataFile)) {
  print_help(opt_parser)
  stop('A tab-delimited raw counts file must be supplied (countDataFile).', call.=FALSE)
} else if (is.null(opt$metaDataFile)) {
  print_help(opt_parser)
  stop('A tab-delimited sample information file must be supplied (metaDataFile).', call.=FALSE)
} else if (is.null(opt$contrastFile)) {
  print_help(opt_parser)
  stop('The tab-delimited file specifying contrasts to be performed must supplied (contrastFile).', call.=FALSE)
} else
  message('Inputs detected: countDataFile, metaDataFile, contrastFile. Continue...')

print_options(opt)

message('loading libraries begins')
suppressPackageStartupMessages(library('DESeq2', character.only=TRUE))
suppressPackageStartupMessages(library('data.table', character.only=TRUE))
suppressPackageStartupMessages(library('BiocParallel', character.only=TRUE))

multicore_param <- MulticoreParam(workers=opt$threads)
register(multicore_param, default=TRUE)

message('loading libraries complete')

#convert fc and padj options to numeric values
fc <- as.numeric(opt$foldChange)
pval <- as.numeric(opt$adjustedPValue)

#read in countDataFile, convert to integers if not already converted
countData <- read.table(file=opt$countDataFile, header = TRUE, sep = ',', row.names = 1, strip.white = TRUE, quote = '', stringsAsFactors = FALSE)
countData <- round(countData)

#read in metaDataFile
colData <- read.table(file =opt$metaDataFile , header = TRUE, sep = '\t', stringsAsFactors = FALSE, na.strings = '', check.names = FALSE)
colData[is.na(colData)] <- 'other' # repliace NA values with 'other'

#read in contrastFile
contrastData <- read.table(file =opt$contrastFile, header=TRUE, sep = '\t', stringsAsFactors = FALSE, strip.white = TRUE)

####
# Create output directories
####
#create directories
cat('creating parent output directories\n')

#create output dir
dir.create(path = opt$outDir, recursive = TRUE) #creates relative
out_dir <- normalizePath(opt$outDir)

#create output subdirs
# counts directory
dir.create(path = paste(out_dir, opt$countsDir, sep = '/'), showWarnings = TRUE, recursive = TRUE, mode = '0777')
# gene_lists directory
dir.create(path = paste(out_dir, opt$geneListsDir, sep = '/'), showWarnings = TRUE, recursive = TRUE, mode = '0777')

#get subdir paths
countsDir <- normalizePath(paste(out_dir, opt$countsDir, sep = '/'))
diffexDir <- normalizePath(paste(out_dir, opt$geneListsDir, sep ='/'))

#create directories
cat('Directories created...\noutput directory base:', out_dir,'\ncounts directory base:', countsDir, '\ngene lists directory base:', diffexDir)

cat('creating contrast-specific (comparison-specific) output directories\n')
for (i in 1:nrow(contrastData)){
  #dir creation
  dir_diffex <- paste(diffexDir, contrastData$factor[i], sep = '/')
  cat('creating ', dir_diffex, '\n')
  dir.create(dir_diffex, showWarnings = FALSE, recursive = TRUE, mode = '0777') # gene_lists/contrast-specific directory
}

####
# Create dds, normalize and produce PCA plot, dispersion, and heatmap for checking
####

message('Initializing DESeq2 result')
#create DESeqDataSet from matrix and filter lowly expressed genes
# want to account for all factors/interactions
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ combinatoric_group)
# NOTE: Should this be more stringent?
dds <- dds[ rowSums(counts(dds)) > 1, ]

rld <- rlog(dds, blind = FALSE)

message('Calculating dispersion')
#Normalize and calculate dispersions
dds <- DESeq(dds, betaPrior = TRUE, parallel = TRUE)

rawCounts <- counts(dds) #raw
normCounts <- counts(dds, normalized = TRUE) #raw/lib size factors

message('Writing count files')
rawCountsDf <- as.data.frame(rawCounts)
data.table::setDT(rawCountsDf, keep.rownames = TRUE)[] #set rownames to valid column
setnames(rawCountsDf, 'rn', 'id')
write.table(x = rawCountsDf, file=paste(countsDir,'raw_counts.txt', sep = '/'),
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

normCountsDf <- as.data.frame(normCounts)
data.table::setDT(normCountsDf, keep.rownames = TRUE)[] #set rownames to valid column
setnames(normCountsDf, 'rn', 'id')
write.table(x = normCountsDf, file=paste(countsDir,'depth_normalized_counts.txt', sep = '/'),
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

rldDf <- as.data.frame(assay(rld))
data.table::setDT(rldDf, keep.rownames = TRUE)[] #set rownames to valid column
setnames(rldDf, 'rn', 'id')
write.table(x = rldDf, file=paste(countsDir,'rlog_normalized_counts.txt', sep = '/'),
    append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


####
# Diffex analysis in loop
####

# Container for the results
de_results = split(contrastData, contrastData$factor)
de_results = lapply(de_results, function(de){split(de, de$base_file_name)})

for (i in 1:nrow(contrastData)){

    #dir assignment
    dir_diffex <- paste(diffexDir, contrastData$factor[i], sep = '/')

    #collect references, etc
    factorName <- contrastData$factor[i]
    testName <- contrastData$test_level[i]
    referenceName <- contrastData$reference_level[i]
    baseName <- contrastData$base_file_name[i]

    message(sprintf('Testing %s: %s vs %s', factorName, testName, referenceName))
    #diffex analysis use generic name, subset DEGs

    #collect test samples in a list, from groups column
    testSamples <- subset(
        x = colData,
        subset = colData[, factorName] == testName & !duplicated(colData$combinatoric_group),
        select = combinatoric_group)
    # put 'combinatoric_group' at beginning of name
    testSamples$combinatoric_group <- sub('^', 'combinatoric_group', testSamples$combinatoric_group)
    # replace '^' with '.' to match name in deseq obj
    testSamples$combinatoric_group <- gsub('^', '.', testSamples$combinatoric_group, fixed = TRUE)
    # replace '-' with '.' to match name in deseq obj
    testSamples$combinatoric_group <- gsub('-', '.', testSamples$combinatoric_group, fixed = TRUE)


    #collect reference samples in a list
    referenceSamples <- subset(
        x = colData,
        subset = colData[,factorName] == referenceName & !duplicated(colData$combinatoric_group),
        select = combinatoric_group)
    # put 'combinatoric_group' at beginning of name
    referenceSamples$combinatoric_group <- sub('^', 'combinatoric_group', referenceSamples$combinatoric_group)
    # replace '^' with '.' to match name in deseq obj
    referenceSamples$combinatoric_group <- gsub('^', '.', referenceSamples$combinatoric_group, fixed = TRUE)
    # replace '-' with '.' to match name in deseq obj
    referenceSamples$combinatoric_group <- gsub('-', '.', referenceSamples$combinatoric_group, fixed = TRUE)


    # Grab the results
    res <- results(
        dds,
        parallel = TRUE,
        contrast = list(c(unlist(testSamples)), c(unlist(referenceSamples))),
        listValues = c(1/nrow(testSamples), -1/nrow(referenceSamples)))
    res <- res[order(res$padj),]


    diffexData <- as.data.frame(res)
    #set rownames to valid column
    data.table::setDT(diffexData, keep.rownames = TRUE)[]
    #rename first column
    colnames(diffexData) <- c('id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')
    diffexData$Condition <- testName
    diffexData$Control <- referenceName

    #make DEG calls and select DEGs
    diffexData$Call <- rep('NO', nrow(diffexData))
    diffexData$Call[which(diffexData$padj <= pval & abs(diffexData$log2FoldChange) >= log2(fc))] = 'YES'
    diffexData <- diffexData[order(-rank(diffexData$Call), diffexData$pvalue), ]

    # Add the individual diffexData to the collection
    de_results[[factorName]][[baseName]] = list(gene_list = diffexData)

    #write to individual tab-delimited txt files
    message('Writing diffex genes')
    write.table(
        x = diffexData,
        file=paste0(dir_diffex,'/',baseName, '.txt'),
        append = FALSE, sep = '\t', na = 'NA', row.names = FALSE, quote = FALSE)
}

save_list = c('dds', 'rld', 'rawCounts', 'normCounts', 'de_results')
deseq2_rda_file = sprintf('%s/deseq2_data.rda', out_dir)
message(sprintf('Saving deseq RData to %s', deseq2_rda_file))
save(list = save_list, file = deseq2_rda_file)

cat('deseq2_diffex.R done.\n')
