### Manjusha Pande, mpande@umich.edu, September 16, 2014
########################################################
### Usage: Rscript <path to the mergeCuffDiffCountFiles.R file> DEDir=<dir> 
#####################################################################################################
### Example:
## $ Rscript /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/mergeCuffDiffCountFiles.R DEDir="/ccmb/BioinfCore/Projects/Spence_959_spencejr_mceachin_RS4_all_compare/cuffdiff" 
#####################################################################################################

args <- commandArgs(T)

for (i in args) {
	arg = strsplit(i,"=",fixed=TRUE)
	assign(arg[[1]][1],arg[[1]][2])	
}

if (!exists("DEDir")) {
	cat (paste ("Incorrect input.\n Usage: Rscript mergeCuffDiffCountFiles.R  DEDir=<dir> \n")) 
} else {
	#cat (paste("DEDir is ", DEDir, "\n"))
}

DESubdirs = list.dirs(DEDir, recursive = FALSE)

countFiles = c()
for (dir in 1:length(DESubdirs)) {
	files = list.files(DESubdirs[dir])
	countFile = files[grep(glob2rx("*_repRawCounts.txt"), files)]
	#cat (paste ("Dir is", DESubdirs[dir], ", countFile is", countFile, "\n"))
	#cat (paste ("Length of countFile is", length(countFile), "\n"))
	if (length(countFile) != 0) {countFiles = c(countFiles, paste(DESubdirs[dir], countFile, sep="/"))}
}

#cat("Creating the raw counts matrix\n")
rawCountMatrix = read.table(countFiles[1], header=T, sep="\t")
if (length(countFiles) > 1 ){
	for (f in 2:length(countFiles)) {
		df = (read.table(countFiles[f], header=T, sep="\t"))
		commonCol = c()
		for (c1 in 2:ncol(df)) { 
			for (c2 in 2:ncol(rawCountMatrix)) {
				if (colnames(df)[c1] == colnames(rawCountMatrix)[c2]) {
					commonCol = c(commonCol, c1)
				}
			}
		}
		df = df[, -commonCol]
		rawCountMatrix = merge(rawCountMatrix, df, by=1)
	}
}
sampleCols = colnames(rawCountMatrix)[3:ncol(rawCountMatrix)]
rawCountMatrix = rawCountMatrix[c(colnames(rawCountMatrix[1:2]), sort(sampleCols))]
write.table (rawCountMatrix, file=paste(DEDir, "CuffDiff_rawCountsMatrix.txt", sep="/"), sep="\t",quote = FALSE, row.names=F)


