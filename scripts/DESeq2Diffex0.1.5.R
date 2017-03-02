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

library(optparse)
#create parser object
option_list = list(
  make_option(c("-c", "--countDataFile"), type = "character", default = NULL,
              help = "Name of file containing tab-delimited RAW counts."),
  make_option(c("-m", "--metaDataFile"), type = "character", default = NULL,
              help = "Name of file containing tab-delimited sample information, ie. factors, levels, etc."),
  make_option(c("-f", "--contrastFile"), type = "character", default = NULL,
              help = "Factor to compare from metaDataFile, must be a column name from metaDataFile"),
  make_option(c("--foldChange"), type = "character", default = NULL,
              help = "Absolute numeric value above which a gene is considered for differential expression."),
  make_option(c("--adjustedPValue"), type = "character", default = NULL,
              help = "Adjusted p-value below which a genes is considered for differential expression.")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#check for necessary files
if (is.null(opt$countDataFile)) {
  print_help(opt_parser)
  stop("A tab-delimited raw counts file must be supplied (countDataFile).", call.=FALSE)
} else if (is.null(opt$metaDataFile)) {
  print_help(opt_parser)
  stop("A tab-delimited sample information file must be supplied (metaDataFile).", call.=FALSE)
} else if (is.null(opt$contrastFile)) {
  print_help(opt_parser)
  stop("The tab-delimited file specifying constats to be performed must supplied (contrastFile).", call.=FALSE)
} else if (is.null(opt$foldChange)) {
  print_help(opt_parser)
  stop("The fold-change cutoff is not specified (foldChange)", call.=FALSE)
} else if (is.null(opt$adjustedPValue)) {
  print_help(opt_parser)
} else
  print("countDataFile, metaDataFile, contrastFile, foldChange, and adjustedPValue detected. Continue.")

#convert fc and padj options to numeric values
fc <- as.numeric(opt$foldChange)
pval <- as.numeric(opt$adjustedPValue)

####
#Load some libraries, read in and parse files, create directory for output
####

#load
library(DESeq2)
library(data.table)
library(BiocParallel)
register(MulticoreParam(6))
options(java.parameters = "-Xmx16000m")
library(xlsx)
library(geneplotter)
library(calibrate)
library(GGally)

#read in countDataFile, convert to integers if not already converted
countData <- read.table(file=opt$countDataFile, header = TRUE, sep = "\t", row.names = 1, strip.white = TRUE, stringsAsFactors = FALSE)
countData <- round(countData)

#read in metaDataFile
colData <- read.table(file = opt$metaDataFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#read in contrastFile
contrastData <- read.table(file = opt$contrastFile, header=TRUE, sep = "\t", stringsAsFactors = FALSE, strip.white = TRUE)
setDT(contrastData)[, c("filePath","fileName") := tstrsplit(file_name, "/", type.convert = TRUE, fixed = TRUE)] #split fileName column to path and name.

#create directories for plots and normalized data
dir.create("plots", showWarnings = TRUE, recursive = FALSE, mode = "0777") #create directory for output
dir.create("normalizedData", showWarnings = TRUE, recursive = FALSE, mode = "0777") #create directory for output

####
# Create dds, normalize and produce PCA plot, dispersion, and heatmap for checking
####

#create DESeqDataSet from matrix and filter lowly expressed genes
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~combinatoricGroups) #want to account for all factors/interactions
dds <- dds[ rowSums(counts(dds)) > 1, ]

#rlog normalization and PCA plot
library(ggplot2)
rld <- rlog(dds, blind = FALSE)
pdf(file = paste0('plots','/','PCA.pdf'), onefile = TRUE)
for (i in names(colData[1:length(colData)])) {
  p <- plotPCA(rld, intgroup = i) + ggtitle(label = as.character(i)) + theme(plot.title = element_text(hjust = 0.5))
  plot(p)
}
dev.off()

#Normalize and calculate dispersions
dds <- DESeq(dds, betaPrior = TRUE, parallel = TRUE)

#get raw counts values for later
rawCounts <- counts(dds)
normCounts <- counts(dds, normalized = TRUE)
idx.nz <- apply(rawCounts, 1, function(x) { all(x > 0)})

#plot dispersions
pdf(file = paste0('plots','/','Dispersion.pdf'), onefile = TRUE)
disp <- plotDispEsts(dds)
dev.off()

#heatmap of normalized data, sample distibution matrix
library(pheatmap)
library(RColorBrewer)
sampleDists <- dist(t(assay(rld))) #rlog normalized data
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$combinatoricGroups, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
pdf(file = paste0('plots','/','CorrelHeatmap.pdf'), onefile = TRUE)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

#heatmap with top 20 expressed genes, rlog normalized data
select <- order(rowMeans(assay(rld)), decreasing=TRUE)[1:30]
df <- data.frame(Group = colData(rld)[,c("combinatoricGroups")], row.names = rownames(colData(dds)))
pdf(file = paste0('plots','/','Top30Heatmap.pdf'), onefile = TRUE)
pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df, fontsize = 8, las = 2)
dev.off()

#boxplot of non-normalized and normalized data
pdf(file = paste0('plots','/','Boxplot.pdf'), onefile = TRUE)
names <- as.factor(colData$combinatoricGroups)
boxplot(log2(rawCounts), range = 0, las = 2, main = "Non-normalized", ylab = "Log2 expression", col= names, fontsize = 8)
boxplot(assay(rld), range = 0, las =2, main = "Rlog normalized", ylab = "Log2 normalized expression", col = names, fontsize = 8)
dev.off()

#raw and normalized count density, removing rows with 0 values
pdf(file = paste0('plots','/','Density.pdf'), onefile = TRUE)
multidensity(log10(rawCounts[idx.nz,]), xlab="mean counts", xlim=c(0,10))#raw counts
multidensity(log10(normCounts[idx.nz,]), xlab="mean counts", xlim=c(0, 10))#normalized (counts/size factors)
dev.off()

#correlation plot between samples --- needs to be revisited because sometime too many samples
pdf(file = paste0('plots','/','RLEmatrix.pdf'), onefile = TRUE)
ggpairs(data = normalizedData[2:length(ncol(normalizedData))], upper = list(continuous = wrap(ggally_cor, use = "pairwise.complete.obs", method = "spearman")), title = "Relative Log2 Expression Correlation Matrix")
dev.off()

#write out normalized values, rlog normalized
normalizedData <- as.data.frame(assay(rld)) #get normalized data from deseq object
setDT(normalizedData, keep.rownames = TRUE)[] #set rownames to valid column
write.table(x = normalizedData, file=paste0('normalizedData',"/","normalizedExpData.txt"), append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE) #write file

####
# Diffex analysis in loop
####

for (i in 1:nrow(contrastData)){
  #newDir creation
  newDir <- paste(as.character(contrastData$filePath[i])) #create directory with date and time
  dir.create(newDir, showWarnings = TRUE, recursive = FALSE, mode = "0777") #create directory for output
  
  #collect references, etc
  referenceName <- contrastData$reference_level[i] 
  testName <- contrastData$test_level[i]
  factorName <- contrastData$factor[i]
  
  #diffex analysis use generic name, subset DEGs
  testSamples <- subset(x = colData, subset = colData[,factorName] == testName & !duplicated(colData$combinatoricGroups), select = combinatoricGroups) #collect test samples in a list, from groups column
  testSamples$combinatoricGroups <- sub("^", "combinatoricGroups", testSamples$combinatoricGroups)
  testSamples$combinatoricGroups <- gsub("^", ".", testSamples$combinatoricGroups, fixed = TRUE)
  referenceSamples <- subset(x = colData, subset = colData[,factorName] == referenceName & !duplicated(colData$combinatoricGroups), select = combinatoricGroups) #collect reference samples in a list
  referenceSamples$combinatoricGroups <- sub("^", "combinatoricGroups", referenceSamples$combinatoricGroups)
  referenceSamples$combinatoricGroups <- gsub("^", ".", referenceSamples$combinatoricGroups, fixed = TRUE)
  res <- results(dds, parallel = TRUE, contrast = list(c(unlist(testSamples)), c(unlist(referenceSamples))), listValues = c(1/nrow(testSamples), -1/nrow(referenceSamples)))
  res <- res[order(res$padj),] #order by adjusted p-value
  diffexData <- as.data.frame(res)
  setDT(diffexData, keep.rownames = TRUE)[] #set rownames to valid column
  colnames(diffexData) <- c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj") #rename first column
  diffexData$Condition <- testName
  diffexData$Control <- referenceName
  
  #diffex call
  for (i in 1:nrow(diffexData)) {
    if (abs(diffexData$log2FoldChange) >= log2(opt$foldChange) & diffexData$padj <= opt$adjustedPValue)
      diffexData$Call[i] <- "YES" #Yes if differentiall expressed
    else
      exp.data$dot[i] <- "NO" #red if significant but not in top 20 dataframe
  }
  
  #select DEGs
  diffexDataDEG <- subset(x = diffexData, subset = abs(diffexData$log2FoldChange) >= log2(opt$foldChange) & diffexData$padj <= opt$adjustedPValue)
  
  #create MA plot
  mainTitle <- list(bquote(paste( "Log"[2], " Fold-Change vs. Mean of normalized expression:")), bquote(paste(~ .(contrastData$fileName[1]))))
  pdf(file = paste0('plots','/','MAplot_',contrastData$fileName[i],'.pdf'), onefile = TRUE)
  plotMA(res)
  mtext(do.call(expression, mainTitle), side=3, line = 1:0)
  dev.off()
  
  #create volcano plot
  mainTitle <- list(bquote(paste( "-Log"[10], "adjusted P-value vs. Log"[2], "Fold-Change:")), bquote(paste(~ .(contrastData$fileName[1]))))
  pdf(file = paste0('plots','/','Volcanoplot_',contrastData$fileName[i],'.pdf'), onefile = TRUE)
  with(diffexData, plot(log2FoldChange, -log10(padj), pch=20)) #plain volcano plot
  with(subset(diffexData, padj<0.05), points(log2FoldChange,  -log10(padj), pch=20, col="red")) #padj <0.05, red
  with(subset(diffexData, abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="orange")) # |log2FC| > 1, orange
  with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="green")) # both, green
  mtext(do.call(expression, mainTitle), side=3, line = 1:0)
  dev.off()
  
  #write to individual tab-delimited txt files, named "diffExpData.[comparison].txt" and a comma-delimited format.
  write.table(x = diffexData, file=paste0(newDir,"/","diffExpData.",contrastData$fileName[i], ".txt"), append = FALSE, sep = "\t", na = 'NA', row.names = FALSE, quote = FALSE)
  write.table(x = diffexData, file=paste0(newDir,"/","diffExpData.",contrastData$fileName[i], ".csv"), append = FALSE, sep = ",", na = 'NA', row.names = FALSE, quote = FALSE)

  #write to an excel file, named "diffExpData.[comparison].txt"
  write.xlsx2(x = diffexData, file=paste0(newDir,"/","diffExpData.",contrastData$fileName[i], ".xlsx"), append = TRUE, sheetName = paste(contrastData$fileName[1], sep = ""), col.names = TRUE, row.names = FALSE)
  write.xlsx2(x = diffexDataDEG, file=paste0(newDir,"/","diffExpData.",contrastData$fileName[i], ".xlsx"), append = TRUE, sheetName = paste(contrastData$fileName[1], "DEG",sep = ""), col.names = TRUE, row.names = FALSE)
}


