#!/usr/bin/env Rscript

####
#Load some libraries, read in and parse files, create directory for output
####

#load
library(optparse)
library(DESeq2)
library(data.table)
library(BiocParallel)
register(MulticoreParam(6))
options(java.parameters = "-Xmx16000m")
library(xlsx)
library(genefilter)
library(geneplotter)
library(ggfortify)
library(ggplot2)
library(calibrate)
library(GGally)
library(reshape2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

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

#read in countDataFile, convert to integers if not already converted
countData <- read.table(file=opt$countDataFile, header = TRUE, sep = "\t", row.names = 1, strip.white = TRUE, stringsAsFactors = FALSE)
countData <- round(countData)

#read in metaDataFile
colData <- read.table(file =opt$metaDataFile , header = TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = '')
colData[is.na(colData)] <- "other" # repliace NA values with 'other'

#read in contrastFile
contrastData <- read.table(file = opt$contrastFile, header=TRUE, sep = "\t", stringsAsFactors = FALSE, strip.white = TRUE)

#create directories for plots and normalized data
dir.create("plots", showWarnings = TRUE, recursive = FALSE, mode = "0777") #create directory for output
dir.create("normalizedData", showWarnings = TRUE, recursive = FALSE, mode = "0777") #create directory for output

####
# Create dds, normalize and produce PCA plot, dispersion, and heatmap for checking
####

#create DESeqDataSet from matrix and filter lowly expressed genes
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~combinatoric_group) #want to account for all factors/interactions
dds <- dds[ rowSums(counts(dds)) > 1, ]

#rlog normalization and PCA plot, shape based on groups, color bas
rld <- rlog(dds, blind = FALSE)
pdf(file = paste0('plots/PCA.pdf'), onefile = TRUE)

#PCA plot for all samples
p.all <- plotPCA(rld, intgroup = 'Name')
CombinatoricGroup <- factor(colData$combinatoric_group)
SampleName <- factor(colData$Name)
gp <- ggplot(p.all$data, aes(x = PC1, y = PC2, color = SampleName, shape = CombinatoricGroup)) + scale_shape_manual(values=1:nlevels(CombinatoricGroup), name = "Combinatoric Group") + geom_point(size=2) + ggtitle(label = as.character('All samples')) + theme(plot.title = element_text(hjust = 0.5)) + guides(colour=guide_legend(nrow=12, title = "Sample"), legend.key = element_rect(size = 1), legend.key.size = unit(0, 'cm')) + theme_classic(base_size = 10) + theme(legend.margin=margin(t = 0, unit='mm')) 
plot(gp)

#get replicate df and sample df
replicateColData <- colData[,grep(".rep", colnames(colData))]
colnames(replicateColData) <- gsub('.rep$','', colnames(replicateColData))
sampleColData <- colData[,-grep(".rep", colnames(colData))]
sampleColData <- sampleColData[,2:ncol(sampleColData)]

#PCA plot for all contrasts
z <- 1
for (i in names(sampleColData[1:length(sampleColData)])) {
  p <- plotPCA(rld, intgroup = i) #get PCA components
  Group <- factor(unlist(sampleColData[z]))
  Replicates <- factor(unlist(replicateColData[z]))
  g <- ggplot(p$data, aes(x = PC1, y = PC2, color = Replicates, shape = Group)) + scale_shape_manual(values=1:nlevels(Group)) + geom_point(size=2) + ggtitle(label = as.character(i)) + theme(plot.title = element_text(hjust = 0.5)) + theme_classic()
  plot(g)
  z <- z + 1
}
dev.off()

#Normalize and calculate dispersions
dds <- DESeq(dds, betaPrior = TRUE, parallel = TRUE)

#get raw counts values for later
rawCounts <- counts(dds) #raw
normCounts <- counts(dds, normalized = TRUE) #raw/lib size factors
idx.nz <- apply(rawCounts, 1, function(x) { all(x > 0)})

#plot dispersions
pdf(file = paste0('plots/Dispersion.pdf'), onefile = FALSE)
disp <- plotDispEsts(dds)
dev.off()

#heatmap of normalized data, sample distibution matrix
sampleDists <- dist(t(assay(rld))) #rlog normalized data
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$combinatoric_group, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pdf(file = paste0('plots/CorrelHeatmap.pdf'), onefile = FALSE)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

#heatmap with top 500 variant or expressed genes, rlog normalized data
colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
select <- order(rowVars(assay(rld)), decreasing=TRUE)[1:500]
df <- data.frame(Group = colData(rld)[,c("combinatoric_group")], row.names = rownames(colData(dds)))
pdf(file = paste0('plots/TopVarHeatmap.pdf'), onefile = FALSE)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, fontsize = 7, las = 2, fontsize_row = 7, color = colors, main = '500 Top Variably Expressed Genes Heatmap')
dev.off()

select <- order(rowMeans(assay(rld)), decreasing=TRUE)[1:500]
df <- data.frame(Group = colData(rld)[,c("combinatoric_group")], row.names = rownames(colData(dds)))
pdf(file = paste0('plots/TopExpHeatmap.pdf'), onefile = FALSE)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, fontsize = 7, las = 2, fontsize_row = 7, color = colors, main = '500 Top Expressed Genes Heatmap')
dev.off()

#boxplot of non-normalized and normalized data
pdf(file = paste0('plots/Boxplot.pdf'), onefile = TRUE)
rawCountsDf <- as.data.frame(rawCounts)
df <- melt(log2(rawCountsDf), variable.name = "Samples", value.name = "count") # reshape the matrix 
df$Condition <- colData$combinatoric_group[match(df$Samples,colData$Name)]
ggplot(df, aes(x = df$Samples, y = count, fill = Condition)) + geom_boxplot(notch = TRUE, outlier.shape = NA) + ggtitle('Non-normalized Counts') + xlab("") + ylab(expression(paste(Log[2]," counts"))) + theme_classic() + theme(axis.text.x  = element_text(angle=90, vjust=0.5))

normCountsDf <- as.data.frame(normCounts)
dfn <- melt(log2(normCountsDf), variable.name = "Samples", value.name = "count") # reshape the matrix 
dfn$Condition <- colData$combinatoric_group[match(dfn$Samples,colData$Name)]
ggplot(dfn, aes(x = dfn$Samples, y = count, fill = Condition)) + geom_boxplot(notch = TRUE, outlier.shape = NA) + ggtitle('Depth-normalized Counts') + xlab("") + ylab(expression(paste(Log[2]," depth-normalized counts"))) + theme_classic() + theme(axis.text.x  = element_text(angle=90, vjust=0.5))

rldDf <- as.data.frame(assay(rld))
dfr <- melt(rldDf, variable.name = "Samples", value.name = "count") # reshape the matrix 
dfr$Condition <- colData$combinatoric_group[match(dfr$Samples,colData$Name)]
ggplot(dfr, aes(x = dfr$Samples, y = count, fill = Condition)) + geom_boxplot(notch = TRUE, outlier.shape = NA) + ggtitle('Rlog-normalized Counts') + xlab("") + ylab(expression(paste(Regularized-Log[2]," normalized counts"))) + theme_classic() + theme(axis.text.x  = element_text(angle=90, vjust=0.5))
dev.off()

#raw and normalized count density, removing rows with 0 values
pdf(file = paste0('plots/Density.pdf'), onefile = TRUE)
df <- as.data.frame(log2(rawCounts[idx.nz,])) # raw counts, removed 0s
df <- melt(df, variable.name = "Samples", value.name = "count") # reshape the matrix 
ggplot(df, aes(x = count, colour = Samples)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.5, size = 0.25)  +
  theme(legend.position = "right") + ylab('Density') + xlab(expression(paste(Log[2]," counts"))) + ggtitle('Non-normalized Counts') + theme_classic()

dfn <- as.data.frame(log2(normCounts[idx.nz,])) #normalized counts (counts/size factors)
dfn <- melt(dfn, variable.name = "Samples", value.name = "count") # reshape the matrix 
ggplot(dfn, aes(x = count, colour = Samples)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.5, size = 0.25)  +
  theme(legend.position = "right") + ylab('Density') + xlab(expression(paste(Log[2]," depth-normalized counts"))) + ggtitle('Depth-normalized Counts') + theme_classic()

dfr <- as.data.frame(rldDf[idx.nz,]) #normalized counts (counts/size factors)
dfr <- melt(dfr, variable.name = "Samples", value.name = "count") # reshape the matrix 
ggplot(dfr, aes(x = count, colour = Samples)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.5, size = 0.25)  +
  theme(legend.position = "right") + ylab('Density') + xlab(expression(paste(Regularized-Log[2]," normalized counts"))) + ggtitle('Rlog-normalized Counts') + theme_classic()
dev.off()

#correlation plot between samples 
uniqGroups <- unique(colData$combinatoric_group) # identify unique groups
pdf(file = paste0('plots/RLEmatrix.pdf'), onefile = TRUE)
for (i in 1:length(uniqGroups)){
  g <- colData[colData$combinatoric_group %in% uniqGroups[i], 'Name'] # collect sample names in those groups
  gc <- rldDf[, colnames(rldDf) %in% g] # pull out columns for those groups
  gcp <- ggpairs(data = gc[1:ncol(gc)], upper = list(continuous = wrap(ggally_cor, use = "pairwise.complete.obs", method = "spearman")), title = as.character(uniqGroups[i])) + theme_bw()
  print(gcp)
}
ggpairs(data = rldDf[1:ncol(rldDf)],upper = list(continuous = wrap(ggally_cor, use = "pairwise.complete.obs", method = "spearman")), title = 'All samples') + theme_bw()
dev.off()

#write out normalized values, rlog normalized
setDT(rldDf, keep.rownames = TRUE)[] #set rownames to valid column
write.table(x = rldDf, file=paste0('normalizedData/RlogNormalizedExpData.txt'), append = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE) #write file

####
# Diffex analysis in loop
####

for (i in 1:nrow(contrastData)){
  #newDir creation
  newDir <- paste0('./',as.character(contrastData$directory_name[i])) #create directory with date and time
  dir.create(newDir, showWarnings = TRUE, recursive = TRUE, mode = "0777") #create directory for output
  
  #collect references, etc
  referenceName <- contrastData$reference_level[i] 
  testName <- contrastData$test_level[i]
  factorName <- contrastData$factor[i]
  
  #diffex analysis use generic name, subset DEGs
  testSamples <- subset(x = colData, subset = colData[,factorName] == testName & !duplicated(colData$combinatoric_group), select = combinatoric_group) #collect test samples in a list, from groups column
  testSamples$combinatoric_group <- sub("^", "combinatoric_group", testSamples$combinatoric_group)
  testSamples$combinatoric_group <- gsub("^", ".", testSamples$combinatoric_group, fixed = TRUE)
  referenceSamples <- subset(x = colData, subset = colData[,factorName] == referenceName & !duplicated(colData$combinatoric_group), select = combinatoric_group) #collect reference samples in a list
  referenceSamples$combinatoric_group <- sub("^", "combinatoric_group", referenceSamples$combinatoric_group)
  referenceSamples$combinatoric_group <- gsub("^", ".", referenceSamples$combinatoric_group, fixed = TRUE)
  res <- results(dds, parallel = TRUE, contrast = list(c(unlist(testSamples)), c(unlist(referenceSamples))), listValues = c(1/nrow(testSamples), -1/nrow(referenceSamples)))
  res <- res[order(res$padj),] #order by adjusted p-value
  diffexData <- as.data.frame(res)
  setDT(diffexData, keep.rownames = TRUE)[] #set rownames to valid column
  colnames(diffexData) <- c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj") #rename first column
  diffexData$Condition <- testName
  diffexData$Control <- referenceName
  
  #### MA plot and volcano plot
  #assign an calculate values based upon logFC and padj
  df <- diffexData
  df$dot <- rep(3, nrow(df))
  df$dot[which(df$padj <= pval & df$log2FoldChange < 0 & abs(df$log2FoldChange) >= log2(fc))] = 2
  df$dot[which(df$padj <= pval & df$log2FoldChange > 0 & abs(df$log2FoldChange) >= log2(fc))] = 1
  df$sig <- df$dot
  
  #take top 10 up, down, then combine, assign label
  top <- rbind(head(subset(df, df$dot == 1), 10),head(subset(df, df$dot == 2), 10))
  df$label <- rep('', nrow(df))
  df$label[which(df$id %in% top$id)] = df$id
  
  #count the number of significan up and down genes, assign value for legend
  df$dot <- factor(df$dot,labels = c( paste0("Up: ", sum(df$dot == 1)),paste0("Down: ", sum(df$dot == 2)),"NS"))
  
  #MA plot
  pdf(file = paste0('plots/MAplot_',contrastData$base_file_name[i],'.pdf'), onefile = FALSE)
  p <- ggplot(df, aes(x = log2(baseMean+1), y = log2FoldChange)) + geom_point(aes(color = df$dot), size = 1) + theme_classic() + xlab(expression(paste(Log[2]," mean normalized expression"))) + ylab(expression(paste(Log[2]," fold-change")))
  p <- p + scale_color_manual(name = '', values=c("#B31B21", "#1465AC", "darkgray"))
  p <- p + scale_x_continuous(breaks=seq(0, max(log2(df$baseMean+1)), 2)) + geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2), color = c("black", "black", "black"))
  p <- p + geom_label_repel(label = df$label, force = 3, segment.alpha = 0.4) + ggtitle(as.character(contrastData$base_file_name[i]))
  print(p)
  dev.off()
  
  #Volcano plot
  pdf(file = paste0('plots/Volcano_',contrastData$base_file_name[i],'.pdf'), onefile = FALSE)
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) + geom_point(aes(color = df$dot), size = 1) + theme_classic() + xlab(expression(paste(Log[2]," fold-change"))) + ylab(expression(paste(-Log[10]," adjusted p-value")))
  p <- p + scale_color_manual(name = '', values=c("#B31B21", "#1465AC", "darkgray"))
  p <- p + geom_vline(xintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2), color = c("black", "black", "black")) + geom_hline(yintercept = -log10(pval), linetype = 2, color = "black")
  p <- p + geom_label_repel(label = df$label, force = 3, segment.alpha = 0.4) + ggtitle(as.character(contrastData$base_file_name[i]))
  print(p)
  dev.off()
  
  #make DEG calls and select DEGs
  diffexData$Call <- rep("NO", nrow(df))
  diffexData$Call[which(diffexData$padj <= pval & abs(diffexData$log2FoldChange) >= log2(fc))] = 'YES'
  diffexDataDEG <- subset(x = diffexData, subset = diffexData$Call == 'YES')
  
  #write to individual tab-delimited txt files, named "diffExpData.[comparison].txt" and a comma-delimited format.
  write.table(x = diffexData, file=paste0('./',newDir,"/diffExpData.",contrastData$base_file_name[i], ".txt"), append = FALSE, sep = "\t", na = 'NA', row.names = FALSE, quote = FALSE)
  write.table(x = diffexData, file=paste0('./',newDir,"/diffExpData.",contrastData$base_file_name[i], ".csv"), append = FALSE, sep = ",", na = 'NA', row.names = FALSE, quote = FALSE)
  
  #write to an excel file, named "diffExpData.[comparison].txt"
  write.xlsx2(x = diffexData, file=paste0('./',newDir,"/diffExpData.",contrastData$base_file_name[i], ".xlsx"), append = TRUE, sheetName = paste(contrastData$base_file_name[i], sep = ""), col.names = TRUE, row.names = FALSE)
  write.xlsx2(x = diffexDataDEG, file=paste0('./',newDir,"/diffExpData.",contrastData$base_file_name[i], ".xlsx"), append = TRUE, sheetName = paste(contrastData$base_file_name[i], "DEG",sep = ""), col.names = TRUE, row.names = FALSE)
}
