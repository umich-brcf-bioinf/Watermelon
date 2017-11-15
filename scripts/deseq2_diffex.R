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

DEFAULT_JAVA_MEMORY_IN_GB = 8
DEFAULT_PANDOC_MEMORY_IN_GB = 8
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
              help = 'Output directory for differential expression data.'),
  make_option(c('--plotsDir'), type = 'character', default = 'plots',
              help = 'Output directory for plots.'),
  make_option(c('--javaMemoryInGb'), type = 'integer', default = DEFAULT_JAVA_MEMORY_IN_GB,
              help = 'Max memory (Gb) allocated to Java VM (-Xmx).'),
  make_option(c('--pandocMemoryInGb'), type = 'integer', default = DEFAULT_PANDOC_MEMORY_IN_GB,
              help = 'Max memory (Gb) allocated to pandoc.')
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

options(java.parameters = paste0('-Xmx',opt$javaMemoryInGb,'g'))
options(pandoc.stack.size=paste0(opt$pandocMemoryInGb,'g'))
threads <- opt$threads

print_options(opt)

message('loading libraries begins')
suppressPackageStartupMessages(library('plotly', character.only=TRUE))
suppressPackageStartupMessages(library('DESeq2', character.only=TRUE))
suppressPackageStartupMessages(library('data.table', character.only=TRUE))
suppressPackageStartupMessages(library('BiocParallel', character.only=TRUE))
suppressPackageStartupMessages(library('xlsx', character.only=TRUE))
suppressPackageStartupMessages(library('genefilter', character.only=TRUE))
suppressPackageStartupMessages(library('geneplotter', character.only=TRUE))
suppressPackageStartupMessages(library('ggfortify', character.only=TRUE))
suppressPackageStartupMessages(library('ggplot2', character.only=TRUE))
suppressPackageStartupMessages(library('calibrate', character.only=TRUE))
suppressPackageStartupMessages(library('GGally', character.only=TRUE))
suppressPackageStartupMessages(library('reshape2', character.only=TRUE))
suppressPackageStartupMessages(library('ggrepel', character.only=TRUE))
suppressPackageStartupMessages(library('pheatmap', character.only=TRUE))
suppressPackageStartupMessages(library('RColorBrewer', character.only=TRUE))

multicore_param <- MulticoreParam(workers=opt$threads)
register(multicore_param, default=TRUE)

message('loading libraries complete')

#convert fc and padj options to numeric values
fc <- as.numeric(opt$foldChange)
pval <- as.numeric(opt$adjustedPValue)

#read in countDataFile, convert to integers if not already converted
countData <- read.table(file=opt$countDataFile, header = TRUE, sep = '\t', row.names = 1, strip.white = TRUE, stringsAsFactors = FALSE)
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
dir.create(path = paste(out_dir, opt$countsDir, sep = '/'), showWarnings = TRUE, recursive = TRUE, mode = '0777') # counts directory
dir.create(path = paste(out_dir, opt$geneListsDir, sep = '/'), showWarnings = TRUE, recursive = TRUE, mode = '0777') # gene_lists directory
dir.create(path = paste(out_dir, opt$plotsDir, 'summary_plots', sep = '/'), showWarnings = TRUE, recursive = TRUE, mode = '0777') # plots/summary_plots directory
dir.create(path = paste(out_dir, opt$plotsDir, 'comparison_plots', sep = '/'), showWarnings = TRUE, recursive = TRUE, mode = '0777') # plots/comparison_plots directory

#get subdir paths
countsDir <- normalizePath(paste(out_dir, opt$countsDir, sep = '/'))
diffexDir <- normalizePath(paste(out_dir, opt$geneListsDir, sep ='/'))
plotsDir_summary <- normalizePath(paste(out_dir, opt$plotsDir,'summary_plots',  sep = '/'))
plotsDir_comparison <- normalizePath(paste(out_dir, opt$plotsDir, 'comparison_plots', sep = '/'))

#create directories
cat('Directories created...\noutput directory base:', out_dir,'\nplots summary directory:', plotsDir_summary, '\nplots comparison directory:', plotsDir_summary, '\ncounts directory base:', countsDir, '\ngene lists directory base:', diffexDir)

cat('creating contrast-specific (comparison-specific) output directories\n')
for (i in 1:nrow(contrastData)){
  #dir creation
  dir_diffex <- paste(diffexDir, contrastData$factor[i], sep = '/') 
  dir_plots <- paste(plotsDir_comparison, contrastData$factor[i], sep = '/') 
  cat('creating ', dir_diffex, '\n')
  dir.create(dir_diffex, showWarnings = FALSE, recursive = TRUE, mode = '0777') # gene_lists/contrast-specific directory
  cat('creating ', dir_plots, '\n')
  dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE, mode = '0777') # plots/comparison_plots/contrast-specific directory
}

####
# Create dds, normalize and produce PCA plot, dispersion, and heatmap for checking
####

cat('initializing DESeq2 result\n')
#create DESeqDataSet from matrix and filter lowly expressed genes
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~combinatoric_group) #want to account for all factors/interactions
dds <- dds[ rowSums(counts(dds)) > 1, ]

cat('building PCA plots\n')
#rlog normalization and PCA plot
rld <- rlog(dds, blind = FALSE)
pdf(file = paste(plotsDir_summary,'PCAplot_All.pdf', sep = '/'), onefile = TRUE)

message('calculating dispersion')
#Normalize and calculate dispersions
dds <- DESeq(dds, betaPrior = TRUE, parallel = TRUE)

#get raw counts values for later
rawCounts <- counts(dds) #raw
normCounts <- counts(dds, normalized = TRUE) #raw/lib size factors
idx.nz <- apply(rawCounts, 1, function(x) { all(x > 0)})

#PCA plot for Rlog-Normalized counts for all samples
p.all <- plotPCA(rld, intgroup = 'sample_name')
CombinatoricGroup <- factor(colData$combinatoric_group)
SampleName <- factor(colData$sample_name)
gp <- ggplot(p.all$data, aes(x = PC1, y = PC2, color = SampleName, shape = CombinatoricGroup)) + xlab(p.all$labels[2]) + ylab(p.all$labels[1]) + scale_shape_manual(values=1:nlevels(CombinatoricGroup), name = 'Combinatoric Group') + geom_point(size=2) + ggtitle(label = as.character('All samples Rlog-Normalized')) + theme(plot.title = element_text(hjust = 0.5)) + guides(colour=guide_legend(nrow=12, title = 'Sample'), legend.key = element_rect(size = 1), legend.key.size = unit(0, 'cm')) + theme_classic(base_size = 10) + theme(legend.margin=margin(t = 0, unit='mm'))
plot(gp)

#PCA for Depth-Normalized counts for all samples
DNC <- SummarizedExperiment(log2(counts(dds, normalized = TRUE)), colData=colData(dds))
p.DNC <- plotPCA(DESeqTransform(DNC), intgroup = 'sample_name')
CombinatoricGroup <- factor(colData$combinatoric_group)
SampleName <- factor(colData$sample_name)
gpDNC <- ggplot(p.DNC$data, aes(x = PC1, y = PC2, color = SampleName, shape = CombinatoricGroup)) + xlab(p.DNC$labels[2]) + ylab(p.DNC$labels[1]) + scale_shape_manual(values=1:nlevels(CombinatoricGroup), name = 'Combinatoric Group') + geom_point(size=2) + ggtitle(label = as.character('All samples Depth-Normalized')) + theme(plot.title = element_text(hjust = 0.5)) + guides(colour=guide_legend(nrow=12, title = 'Sample'), legend.key = element_rect(size = 1), legend.key.size = unit(0, 'cm')) + theme_classic(base_size = 10) + theme(legend.margin=margin(t = 0, unit='mm'))
plot(gpDNC)

#PCA for Raw counts for all samples
RC <- SummarizedExperiment(log2(counts(dds, normalized = FALSE)), colData=colData(dds))
p.RC <- plotPCA(DESeqTransform(RC), intgroup = 'sample_name')
CombinatoricGroup <- factor(colData$combinatoric_group)
SampleName <- factor(colData$sample_name)
gpRC <- ggplot(p.RC$data, aes(x = PC1, y = PC2, color = SampleName, shape = CombinatoricGroup)) + xlab(p.RC$labels[2]) + ylab(p.RC$labels[1]) + scale_shape_manual(values=1:nlevels(CombinatoricGroup), name = 'Combinatoric Group') + geom_point(size=2) + ggtitle(label = as.character('All samples Raw')) + theme(plot.title = element_text(hjust = 0.5)) + guides(colour=guide_legend(nrow=12, title = 'Sample'), legend.key = element_rect(size = 1), legend.key.size = unit(0, 'cm')) + theme_classic(base_size = 10) + theme(legend.margin=margin(t = 0, unit='mm'))
plot(gpRC)
dev.off()

#get replicate df and sample df
replicateColData <- colData[,grep('.replicate', colnames(colData))]
colnames(replicateColData) <- gsub('.replicate$','', colnames(replicateColData))
sampleColData <- colData[,-grep('.replicate', colnames(colData))]
sampleColData <- sampleColData[,2:ncol(sampleColData)]

#PCA plot for all contrasts
z <- 1
for (i in names(sampleColData[1:length(sampleColData)])) {
  if (i != 'combinatoric_group'){
    pca_dir <- paste(plotsDir_comparison, i, sep='/')
    dir.create(pca_dir, recursive=TRUE)
    pca_filename <- paste(pca_dir, 'PCAplot.pdf', sep = '/') 
    cat(paste0('building plot: [', pca_filename, ']\n'))
    pdf(file = pca_filename)
    p <- plotPCA(rld, intgroup = i) #get PCA components
    Group <- factor(unlist(sampleColData[z]))
    Replicates <- factor(unlist(replicateColData[z]))
    g <- ggplot(p$data, aes(x = PC1, y = PC2, color = Replicates, shape = Group)) + xlab(p$labels[2]) + ylab(p$labels[1]) + scale_shape_manual(values=1:nlevels(Group)) + geom_point(size=2) + ggtitle(label = as.character(i)) + theme(plot.title = element_text(hjust = 0.5)) + theme_classic()
    plot(g)
    z <- z + 1
    dev.off()
  }
  else{
    message('Skipping combinatoric_group...')
  }
}

message('Interactive html plots')
gp <- ggplotly(gp)
path_name <- file.path(plotsDir_summary, 'PCAplot_All_RlogNormalized.html')
htmlwidgets::saveWidget(gp, file = path_name, selfcontained = TRUE)
delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
unlink(x = delDir, recursive = TRUE, force = TRUE)

gpDNC <- ggplotly(gpDNC)
path_name <- file.path(plotsDir_summary, 'PCAplot_All_DepthNormalized.html')
htmlwidgets::saveWidget(gpDNC, file = path_name, selfcontained = TRUE)
delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
unlink(x = delDir, recursive = TRUE, force = TRUE)

gpRC <- ggplotly(gpRC)
path_name <- file.path(plotsDir_summary, 'PCAplot_All_Raw.html')
htmlwidgets::saveWidget(gpRC, file = path_name, selfcontained = TRUE)
delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
unlink(x = delDir, recursive = TRUE, force = TRUE)

z <- 1
for (i in names(sampleColData[1:length(sampleColData)])) {
  if (i != 'combinatoric_group'){
    p <- plotPCA(rld, intgroup = i) #get PCA components
    Group <- factor(unlist(sampleColData[z]))
    Replicates <- factor(unlist(replicateColData[z]))
    g <- ggplot(p$data, aes(x = PC1, y = PC2, color = Replicates, shape = Group)) + xlab(p$labels[2]) + ylab(p$labels[1]) + scale_shape_manual(values=1:nlevels(Group)) + geom_point(size=2) + ggtitle(label = as.character(i)) + theme(plot.title = element_text(hjust = 0.5)) + theme_classic()
    g <- ggplotly(g)
    path_name <- file.path(plotsDir_comparison, paste(i,'PCAplot.html', sep = '/'))
    message(paste0('saving: [',path_name,']'))
    htmlwidgets::saveWidget(g,file = path_name)
    delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
    unlink(x = delDir, recursive = TRUE, force = TRUE)
    z <- z + 1
  }
  else{
    message('Skipping combinatoric_group...')
  }
}

cat('building plots\n')
cat('\tdispersions\n')
#plot dispersions
pdf(file = paste(plotsDir_summary,'Dispersion.pdf', sep = '/'), onefile = FALSE)
disp <- plotDispEsts(dds)
dev.off()

cat('\tsample heatmap\n')
#heatmap of normalized data, sample distibution matrix
sampleDists <- dist(t(assay(rld))) #rlog normalized data
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$combinatoric_group, sep='-')
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, 'Blues')))(255)
pdf(file = paste(plotsDir_summary,'Heatmap_Samples.pdf', sep = '/'), onefile = FALSE)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

cat('\ttop variant heatmap\n')
#heatmap with top 500 variant or expressed genes, rlog normalized data
colors <- colorRampPalette(brewer.pal(9, 'Blues'))(255)
select <- order(rowVars(assay(rld)), decreasing=TRUE)[1:500]
df <- data.frame(Group = colData(rld)[,c('combinatoric_group')], row.names = rownames(colData(dds)))
pdf(file = paste(plotsDir_summary,'/Heatmap_TopVar.pdf', sep = '/'), onefile = FALSE, width=10, height=20)
pheatmap(assay(rld)[select,], scale="row",  cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, fontsize = 7, las = 2, fontsize_row = 7, color = colors, main = '500 Top Variably Expressed Genes Heatmap')
dev.off()

select <- order(rowMeans(assay(rld)), decreasing=TRUE)[1:500]
df <- data.frame(Group = colData(rld)[,c('combinatoric_group')], row.names = rownames(colData(dds)))
pdf(file = paste(plotsDir_summary,'Heatmap_TopExp.pdf', sep = '/'), onefile = FALSE, width=10, height=20)
pheatmap(assay(rld)[select,], scale="row",  cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, fontsize = 7, las = 2, fontsize_row = 7, color = colors, main = '500 Top Expressed Genes Heatmap')
dev.off()

cat('\tbox plots\n')
#boxplot of non-normalized and normalized data
pdf(file = paste(plotsDir_summary,'BoxPlot.pdf', sep = '/'), onefile = TRUE)
rawCountsDf <- as.data.frame(rawCounts)
df <- melt(log2(rawCountsDf), variable.name = 'Samples', value.name = 'count') # reshape the matrix
df$Condition <- colData$combinatoric_group[match(df$Samples,colData$sample_name)]
ggplot(df, aes(x = df$Samples, y = count, fill = Condition)) + geom_boxplot(notch = TRUE, outlier.shape = NA) + ggtitle('Non-normalized Counts') + xlab('') + ylab('Log2 counts') + theme_classic() + theme(axis.text.x  = element_text(angle=90, vjust=0.5))

normCountsDf <- as.data.frame(normCounts)
dfn <- melt(log2(normCountsDf), variable.name = 'Samples', value.name = 'count') # reshape the matrix
dfn$Condition <- colData$combinatoric_group[match(dfn$Samples,colData$sample_name)]
ggplot(dfn, aes(x = dfn$Samples, y = count, fill = Condition)) + geom_boxplot(notch = TRUE, outlier.shape = NA) + ggtitle('Depth-normalized Counts') + xlab('') + ylab('Log2 depth-normalized counts') + theme_classic() + theme(axis.text.x  = element_text(angle=90, vjust=0.5))

rldDf <- as.data.frame(assay(rld))
dfr <- melt(rldDf, variable.name = 'Samples', value.name = 'count') # reshape the matrix
dfr$Condition <- colData$combinatoric_group[match(dfr$Samples,colData$sample_name)]
ggplot(dfr, aes(x = dfr$Samples, y = count, fill = Condition)) + geom_boxplot(notch = TRUE, outlier.shape = NA) + ggtitle('Rlog-normalized Counts') + xlab('') + ylab('Regularized-Log2 normalized counts') + theme_classic() + theme(axis.text.x  = element_text(angle=90, vjust=0.5))
dev.off()

#interactive boxplots
bp <- ggplot(df, aes(x = df$Samples, y = count, fill = Condition)) + geom_boxplot(notch = TRUE, outlier.shape = NA) + ggtitle('Non-normalized Counts') + xlab('') + ylab('Log2 counts') + theme_classic() + theme(axis.text.x  = element_text(angle=90, vjust=0.5))
bp <- ggplotly(bp)
path_name <- file.path(plotsDir_summary, 'BoxPlot_RawCounts.html')
htmlwidgets::saveWidget(bp, file = path_name)
delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
unlink(x = delDir, recursive = TRUE, force = TRUE)

bpn <- ggplot(dfn, aes(x = dfn$Samples, y = count, fill = Condition)) + geom_boxplot(notch = TRUE, outlier.shape = NA) + ggtitle('Depth-normalized Counts') + xlab('') + ylab('Log2 depth-normalized counts') + theme_classic() + theme(axis.text.x  = element_text(angle=90, vjust=0.5))
bpn <- ggplotly(bpn)
path_name <- file.path(plotsDir_summary, 'BoxPlot_DepthNormalizedCounts.html')
htmlwidgets::saveWidget(bpn, file = path_name)
delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
unlink(x = delDir, recursive = TRUE, force = TRUE)

bpr <- ggplot(dfr, aes(x = dfr$Samples, y = count, fill = Condition)) + geom_boxplot(notch = TRUE, outlier.shape = NA) + ggtitle('Rlog-normalized Counts') + xlab('') + ylab('Regularized-Log2 normalized counts') + theme_classic() + theme(axis.text.x  = element_text(angle=90, vjust=0.5))
bpr <- ggplotly(bpr)
path_name <- file.path(plotsDir_summary, 'BoxPlot_RlogNormalizedCounts.html')
htmlwidgets::saveWidget(bpr, file = path_name)
delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
unlink(x = delDir, recursive = TRUE, force = TRUE)

cat('\tdensity plots\n')
#raw and normalized count density, removing rows with 0 values
pdf(file = paste(plotsDir_summary,'DensityPlot.pdf', sep = '/'), onefile = TRUE)
df <- as.data.frame(log2(rawCounts[idx.nz,])) # raw counts, removed 0s
df <- melt(df, variable.name = 'Samples', value.name = 'count') # reshape the matrix
ggplot(df, aes(x = count, colour = Samples)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.5, size = 0.25)  +
  theme(legend.position = 'right') + ylab('Density') + xlab('Log2 counts') + ggtitle('Non-normalized Counts') + theme_classic()

dfn <- as.data.frame(log2(normCounts[idx.nz,])) #normalized counts (counts/size factors)
dfn <- melt(dfn, variable.name = 'Samples', value.name = 'count') # reshape the matrix
ggplot(dfn, aes(x = count, colour = Samples)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.5, size = 0.25)  +
  theme(legend.position = 'right') + ylab('Density') + xlab('Log2 depth-normalized counts') + ggtitle('Depth-normalized Counts') + theme_classic()

dfr <- as.data.frame(rldDf[idx.nz,]) #normalized counts (counts/size factors)
dfr <- melt(dfr, variable.name = 'Samples', value.name = 'count') # reshape the matrix
ggplot(dfr, aes(x = count, colour = Samples)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.5, size = 0.25)  +
  theme(legend.position = 'right') + ylab('Density') + xlab('Regularized-Log2 normalized counts') + ggtitle('Rlog-normalized Counts') + theme_classic()
dev.off()

#interactive density plots
dp <- ggplot(df, aes(x = count, colour = Samples)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.5, size = 0.25)  +
  theme(legend.position = 'right') + ylab('Density') + xlab('Log2 counts') + ggtitle('Non-normalized Counts') + theme_classic()
dp <- ggplotly(dp)
path_name <- file.path(plotsDir_summary, 'DensityPlot_RawCounts.html')
htmlwidgets::saveWidget(dp, file = path_name)
delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
unlink(x = delDir, recursive = TRUE, force = TRUE)

dpn <- ggplot(dfn, aes(x = count, colour = Samples)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.5, size = 0.25)  +
  theme(legend.position = 'right') + ylab('Density') + xlab('Log2 depth-normalized counts') + ggtitle('Depth-normalized Counts') + theme_classic()
dpn <- ggplotly(dpn)
path_name <- file.path(plotsDir_summary, 'DensityPlot_DepthNormalizedCounts.html')
htmlwidgets::saveWidget(dpn, file = path_name)
delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
unlink(x = delDir, recursive = TRUE, force = TRUE)

dpr <- ggplot(dfr, aes(x = count, colour = Samples)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.5, size = 0.25)  +
  theme(legend.position = 'right') + ylab('Density') + xlab('Regularized-Log2 normalized counts') + ggtitle('Rlog-normalized Counts') + theme_classic()
dpr <- ggplotly(dpr)
path_name <- file.path(plotsDir_summary,'DensityPlot_RlogNormalizedCounts.html')
htmlwidgets::saveWidget(dpr, file = path_name)
delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
unlink(x = delDir, recursive = TRUE, force = TRUE)

cat('\tcorrelation plots\n')
#correlation plot between samples
uniqGroups <- unique(colData$combinatoric_group) # identify unique groups
pdf(file = paste(plotsDir_summary,'CorrelMatrix.pdf', sep = '/'), onefile = TRUE, compress = TRUE)
for (i in 1:length(uniqGroups)){
  g <- colData[colData$combinatoric_group %in% uniqGroups[i], 'sample_name'] # collect sample names in those groups
  if (length(g) < 2) next
  gc <- rldDf[, colnames(rldDf) %in% g] # pull out columns for those groups
  gcp <- ggpairs(data = gc[1:ncol(gc)], upper = list(continuous = wrap(ggally_cor, use = 'pairwise.complete.obs', method = 'spearman')), title = as.character(uniqGroups[i])) + theme_bw()
  print(gcp)
}
ggpairs(data = rldDf[1:ncol(rldDf)], upper = list(continuous = wrap(ggally_cor, use = 'pairwise.complete.obs', method = 'spearman')), title = 'All samples') + theme_bw()
dev.off()

#interactive correlation plots
for (i in 1:length(uniqGroups)){
  g <- colData[colData$combinatoric_group %in% uniqGroups[i], 'sample_name'] # collect sample names in those groups
  if (length(g) < 2) next
  gc <- rldDf[, colnames(rldDf) %in% g] # pull out columns for those groups
  gcp <- ggpairs(data = gc[1:ncol(gc)], upper = list(continuous = wrap(ggally_cor, use = 'pairwise.complete.obs', method = 'spearman')), title = as.character(uniqGroups[i])) + theme_bw()
  gcp <- ggplotly(gcp)
  path_name <- file.path(paste0(plotsDir_summary,'/CorrelMatrix_',uniqGroups[i],'.html'))
  htmlwidgets::saveWidget(gcp, file = path_name)
  delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
  unlink(x = delDir, recursive = TRUE, force = TRUE)
}

if(length(1:ncol(rldDf)) < 10){
   cat('\tcorrelation plots:all v all\n')
   acp <- ggpairs(data = rldDf[1:ncol(rldDf)],upper = list(continuous = wrap(ggally_cor, use = 'pairwise.complete.obs', method = 'spearman')), title = 'All samples') + theme_bw()
   acp <- ggplotly(acp)
   path_name <- file.path(plotsDir_comparison,'CorrelMatrix_All.html')
   htmlwidgets::saveWidget(acp, file = path_name)
   delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
   unlink(x = delDir, recursive = TRUE, force = TRUE)
 } else{
   message('Too many samples for all vs all CorrelMatrix. Skipping...')
 }

cat('\twriting counts files\n')
setDT(rawCountsDf, keep.rownames = TRUE)[] #set rownames to valid column
setnames(rawCountsDf, 'rn', 'id')
write.table(x = rawCountsDf, file=paste(countsDir,'raw_counts.txt', sep = '/'), append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE) #write file

setDT(normCountsDf, keep.rownames = TRUE)[] #set rownames to valid column
setnames(normCountsDf, 'rn', 'id')
write.table(x = normCountsDf, file=paste(countsDir,'depth_normalized_counts.txt', sep = '/'), append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE) #write file

setDT(rldDf, keep.rownames = TRUE)[] #set rownames to valid column
setnames(rldDf, 'rn', 'id')
write.table(x = rldDf, file=paste(countsDir,'rlog_normalized_counts.txt', sep = '/'), append = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE) #write file

####
# Diffex analysis in loop
####

for (i in 1:nrow(contrastData)){
  #dir assignment
  dir_diffex <- paste(diffexDir, contrastData$factor[i], sep = '/') 
  dir_plots <- paste(plotsDir_comparison, contrastData$factor[i], sep = '/') 
  
  #collect references, etc
  referenceName <- contrastData$reference_level[i]
  testName <- contrastData$test_level[i]
  factorName <- contrastData$factor[i]
  
  cat('\tcontrasting\n')
  #diffex analysis use generic name, subset DEGs
  testSamples <- subset(x = colData, subset = colData[,factorName] == testName & !duplicated(colData$combinatoric_group), select = combinatoric_group) #collect test samples in a list, from groups column
  testSamples$combinatoric_group <- sub('^', 'combinatoric_group', testSamples$combinatoric_group) # put 'combinatoric_group' at beginning of name
  testSamples$combinatoric_group <- gsub('^', '.', testSamples$combinatoric_group, fixed = TRUE) # replace '^' with '.' to match name in deseq obj
  testSamples$combinatoric_group <- gsub('-', '.', testSamples$combinatoric_group, fixed = TRUE) # replace '-' with '.' to match name in deseq obj
  referenceSamples <- subset(x = colData, subset = colData[,factorName] == referenceName & !duplicated(colData$combinatoric_group), select = combinatoric_group) #collect reference samples in a list
  referenceSamples$combinatoric_group <- sub('^', 'combinatoric_group', referenceSamples$combinatoric_group) # put 'combinatoric_group' at beginning of name
  referenceSamples$combinatoric_group <- gsub('^', '.', referenceSamples$combinatoric_group, fixed = TRUE) # replace '^' with '.' to match name in deseq obj
  referenceSamples$combinatoric_group <- gsub('-', '.', referenceSamples$combinatoric_group, fixed = TRUE)# replace '-' with '.' to match name in deseq obj
  res <- results(dds, parallel = TRUE, contrast = list(c(unlist(testSamples)), c(unlist(referenceSamples))), listValues = c(1/nrow(testSamples), -1/nrow(referenceSamples)))
  res <- res[order(res$padj),] #order by adjusted p-value
  diffexData <- as.data.frame(res)
  setDT(diffexData, keep.rownames = TRUE)[] #set rownames to valid column
  colnames(diffexData) <- c('id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj') #rename first column
  diffexData$Condition <- testName
  diffexData$Control <- referenceName
  
  cat('\tMA plot\n')
  #### MA plot and volcano plot
  #assign an calculate values based upon logFC and padj
  df <- diffexData
  df$dot <- rep(3, nrow(df))
  df$dot[which(df$padj <= pval & df$log2FoldChange < 0 & abs(df$log2FoldChange) >= log2(fc))] = 2
  df$dot[which(df$padj <= pval & df$log2FoldChange > 0 & abs(df$log2FoldChange) >= log2(fc))] = 1
  df$sig <- df$dot
  
  #take top 10 up, down, then combine, assign label
  top <- rbind(head(subset(df, df$dot == 1), 10),head(subset(df, df$dot == 2), 10))
  top$label <- top$id
  df <- merge(x = df, y = top[,c('id','label')], by = "id", all.x = TRUE)
  
  #count the number of significan up and down genes, assign value for legend
  df$dot <- factor(df$dot,levels = c(1,2,3), labels = c(paste0('Up: ', sum(df$dot == 1)),paste0('Down: ', sum(df$dot == 2)),'NS'))
  
  #MA plot
  pdf(file = paste0(dir_plots,'/MAplot_',contrastData$base_file_name[i],'.pdf'), onefile = FALSE)
  p <- ggplot(df, aes(x = log2(baseMean+1), y = log2FoldChange)) + geom_point(aes(color = df$dot), size = 1) + theme_classic() + xlab('Log2 mean normalized expression') + ylab('Log2 fold-change')
  p <- p + scale_color_manual(name = '', values=c('#B31B21', '#1465AC', 'darkgray'))
  p <- p + scale_x_continuous(breaks=seq(0, max(log2(df$baseMean+1)), 2)) + geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2), color = c('black', 'black', 'black'))
  if (sum(is.na(df$label)) < nrow(df)) {
    p <- p + geom_label_repel(label = df$label, force = 3, segment.alpha = 0.4) + ggtitle(as.character(contrastData$base_file_name[i]))
  } else {
    p <- p + ggtitle(as.character(contrastData$base_file_name[i]))
  }
  print(p)
  dev.off()
  
  #interactive MA plot
  p <- ggplot(df, aes(x = log2(baseMean+1), y = log2FoldChange, colour = df$dot, label = id)) + geom_point(size = 1) + theme_classic() + xlab('Log2 mean normalized expression') + ylab('Log2 fold-change')
  p <- p + scale_color_manual(name = '', values=c('#B31B21', '#1465AC', 'darkgray'))
  p <- p + scale_x_continuous(breaks=seq(0, max(log2(df$baseMean+1)), 2)) + geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2), color = c('black', 'black', 'black'))
  p <- p + ggtitle(as.character(contrastData$base_file_name[i]))
  mp <- ggplotly(p)
  path_name <- file.path(paste0(dir_plots,'/MAplot_',contrastData$base_file_name[i],'.html'))
  htmlwidgets::saveWidget(mp, file = path_name)
  delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
  unlink(x = delDir, recursive = TRUE, force = TRUE)
  
  cat('\tvolcano plot\n')
  #Volcano plot
  pdf(file = paste0(dir_plots,'/VolcanoPlot_',contrastData$base_file_name[i],'.pdf'), onefile = FALSE)
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) + geom_point(aes(color = df$dot), size = 1) + theme_classic() + xlab('Log2 fold-change') + ylab('-Log10 adjusted p-value')
  p <- p + scale_color_manual(name = '', values=c('#B31B21', '#1465AC', 'darkgray'))
  p <- p + geom_vline(xintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2), color = c('black', 'black', 'black')) + geom_hline(yintercept = -log10(pval), linetype = 2, color = 'black')
  if (sum(is.na(df$label)) < nrow(df)) {
    p <- p + geom_label_repel(label = df$label, force = 3, segment.alpha = 0.4) + ggtitle(as.character(contrastData$base_file_name[i]))
  } else {
    p <- p + ggtitle(as.character(contrastData$base_file_name[i]))
  }
  print(p)
  dev.off()
  
  #interactive Volcano plot
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), colour = df$dot, label = id)) + geom_point(size = 1) + theme_classic() + xlab('Log2 fold-change') + ylab('-Log10 adjusted p-value')
  p <- p + scale_color_manual(name = '', values=c('#B31B21', '#1465AC', 'darkgray'))
  p <- p + geom_vline(xintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2), color = c('black', 'black', 'black')) + geom_hline(yintercept = -log10(pval), linetype = 2, color = 'black')
  p <- p + ggtitle(as.character(contrastData$base_file_name[i]))
  vp <- ggplotly(p)
  path_name <- file.path(paste0(dir_plots,'/VolcanoPlot_',contrastData$base_file_name[i],'.html'))
  htmlwidgets::saveWidget(vp, file = path_name)
  delDir <- gsub(pattern = '.html', replacement = '_files', x = path_name)
  unlink(x = delDir, recursive = TRUE, force = TRUE)
  
  cat('\twriting diffex genes\n')
  #make DEG calls and select DEGs
  diffexData$Call <- rep('NO', nrow(df))
  diffexData$Call[which(diffexData$padj <= pval & abs(diffexData$log2FoldChange) >= log2(fc))] = 'YES'
  diffexData <- diffexData[order(-rank(diffexData$Call), diffexData$pvalue), ]
  
  #write to individual tab-delimited txt files
  write.table(x = diffexData, file=paste0(dir_diffex,'/',contrastData$base_file_name[i], '.txt'), append = FALSE, sep = '\t', na = 'NA', row.names = FALSE, quote = FALSE)
}

cat('deseq2_diffex.R done.\n')

