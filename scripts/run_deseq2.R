args <- commandArgs(T)
USAGE <- "Rscript <script> outDir=<dir> htseqDir=<dir> sampleConditionsFileName=<file> comparisonsFileName=<file>"
for (i in args) {
  arg = strsplit(i,"=",fixed=TRUE)
  assign(arg[[1]][1],arg[[1]][2])	
}

args_to_validate <- c("outDir", "htseqDir", "sampleConditionsFileName", "comparisonsFileName")

for (arg in args_to_validate) {
  if (!exists(arg)) {
    stop(paste0("Missing ", arg, " option.\nUsage: ", USAGE, "\n")) 
  } else {
    value <- tryCatch({
      normalizePath(get(arg))
    }, warning = function(w) {
      stop(paste0(arg, " not an existing file or dir"))
    })
    assign(arg, value)
  }
  message(paste0(arg, " = ", get(arg)))
}

message("loading libraries...")
suppressMessages(library(DESeq2))

sampleConditionFileName <- "sample_conditions.txt"
comparisonsFileName <- 'comparisons.txt'
htseqDir <- "./07-htseq"
sampleTable <- read.table(file=sampleConditionFileName, sep ="\t", header = TRUE)
countFileName <- paste(sampleTable$Sample_ID, "_counts.txt", sep ="")
sampleTable <- cbind(countFileName, sampleTable)
sampleTable <- sampleTable[,c(2,1,3)]

comparisonsTable <- read.table(file=comparisonsFileName, sep="\t", header=TRUE, stringsAsFactors=FALSE)
dds <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=htseqDir, design=~Condition)
dds <- DESeq(dds)

runContrast <- function (deseqDataset, condition_1, condition_2, outDir) {
  message(paste('Processing', condition_1, 'v', condition_2))
  results <- results(deseqDataset,
                     contrast=c("Condition", condition_1, condition_2))
  results <- results[order(results$padj),]
  df <-data.frame(results)
  geneSymbol <- rownames(df)
  df <- cbind(geneSymbol, df)
  fileName <- paste0(outDir, "/", condition_1, "_", condition_2, ".txt")
  write.table(df, row.names=FALSE, file=fileName, sep="\t", quote=FALSE)
}

for(i in 1:nrow(comparisonsTable)) {
  comparison <- comparisonsTable[i,]
  runContrast(dds, comparison$condition_1, comparison$condition_2, outDir)
}

