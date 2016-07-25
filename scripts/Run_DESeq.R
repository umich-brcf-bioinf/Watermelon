####  Pande, mpande@umich.edu, July 17, 2014
########################################################
### Usage: Rscript <path to the Run_DESeq2.R file> baseDir=<baseDir> grpRepFile=<grpRepFile> countsFile=<countsFile>
#####################################################################################################
### Example:
## $ Rscript /ccmb/BioinfCore/Projects/RNA-seq_pipeline/RSQ_scripts/Run_DESeq2.R baseDir=/ccmb/BioinfCore/Projects/Rsq_test_v3/SEnoNovelTranscripts/DESeq2/WT_v_mutant grpRepFile=/ccmb/BioinfCore/Projects/Rsq_test_v3/SEnoNovelTranscripts/runInfo/Run_2014-05-15_16-49-20/Group_replicate_list.txt countsFile=/ccmb/BioinfCore/Projects/Rsq_test_v3/SEnoNovelTranscripts/HTSeq/HTSeq_counts.txt
#####################################################################################################
# baseDir="/ccmb/BioinfCore/Projects/Rsq_test_v3/SEnoNovelTranscripts/DESeq2/WT_v_mutant_v_E8.5_mutant"
# grpRepFile="/ccmb/BioinfCore/Projects/Rsq_test_v3/SEnoNovelTranscripts/runInfo/Run_2014-05-15_16-49-20/Group_replicate_list.txt" 
# countsFile="/ccmb/BioinfCore/Projects/Rsq_test_v3/SEnoNovelTranscripts/HTSeq/HTSeq_counts.txt"
library(DESeq)
library(DESeq2)
args <- commandArgs(T)

for (i in args) {
	arg = strsplit(i,"=",fixed=TRUE)
	assign(arg[[1]][1],arg[[1]][2])	
}

dataSetPath = strsplit(baseDir, "/")[[1]]
dataSet = dataSetPath[length(dataSetPath)]
#cat(paste("Counts file is", countsFile, "\n", sep=" "))
Grps = unlist(strsplit(dataSet, "_", fixed=TRUE))    #abhasi : "_" instead of "_v_"
grpReps = read.table(grpRepFile, sep="\t")
reps = c()
for (i in 1:length(grpReps[,1])) { 
	tmp = c()
	tmp = gsub ("replicate", grpReps[i,]$V1, grpReps[i,]$V2)
	tmp = unlist (strsplit(tmp, ", ",fixed=T))
	reps = c(reps, tmp)
}
countTable1 = read.table(countsFile, header=T, row.names=1)
rownames = row.names(countTable1)
countTable = sapply(countTable1[,], as.integer) ### convert data frame to integer matrix
row.names(countTable) = rownames
col = colnames(countTable)

## Remove samples that are not included in the analysis
exclude = c() ## samples to be excluded
for (i in 1:length(col)) {  ## for each sample 
	sample_grps = c() ## groups this sample belongs to 
	for (j in 1:length(reps)) { 
		tmp = unlist(strsplit(reps[j], ": ", fixed=T)) 
		if (tmp[2] == col[i]) { 
			grp = sub("_\\d+$", "",tmp[1]) ## sample is a replicate of this group
			if (grp %in% Grps) { ## this group is being compared
				sample_grps = c (sample_grps, grp)
			}
		} 
	}
	if (length(sample_grps) == 0) { #cat (paste (col[i], "does not belong to any of the groups being compared. This sample will be excluded from DESeq2 analysis.\n"))
	exclude = c(exclude, col[i])
 	} else if (length(sample_grps) > 1) {cat (paste (col[i], "belongs to multiple groups: "))
		cat(sample_grps, sep=", "); cat("\n") 
		stop("Overlapping groups cannot be normalized together.\n")
	}
}
if (length(exclude) > 0) {
	cat ("Following sample/s do not belong to any of the groups being compared and will be excluded from DESeq2 analysis.\n")
	for (s in 1:length(exclude)) {
		cat (paste(exclude[s], "\n"))
	}
}
## remove samples to be excluded from countTable

index = which(colnames(countTable)%in%exclude)
if (length(index)) {
	countTable = countTable[, -index]  ### try droplevels??
}
col = colnames(countTable)
for (i in 1:length(col)) {  ## for each sample 
	for (j in 1:length(reps)) { 
		tmp = unlist(strsplit(reps[j], ": ", fixed=T)) 
		grp = sub("_\\d+$", "",tmp[1]) ## sample is a replicate of this group
		if (grp %in% Grps && col[i] == tmp[2]) { 
			col[i] = tmp[1] ## assign replicate number 
		}
	}
}

countTable = countTable[, order(col)] ## rearrange column by groups
col = col[order(col)]
conditions = c()
for (i in 1:length(col)) { 
	cond = sub("_\\d+$", "",col[i]) 
	conditions = c(conditions, cond)
}

load_DESeq2 = require("DESeq2") ## returns FALSE if package is not installed
if(load_DESeq2 == FALSE) {
	stop("DESeq2 package could not be loaded. Please check if the package is installed.")
}

countDataSet = newCountDataSet(countTable, conditions)
countDataSet = estimateSizeFactors(countDataSet) ### estimate size factors from the count data
if (length(unique(conditions))==length(col)) { ## no conditions are replicated
	countDataSet = estimateDispersions(countDataSet, method="blind", sharingMode = "fit-only")
} else {
	countDataSet = estimateDispersions(countDataSet) ### estimate biological variation, default method=pooled, sharingMode=maximum
}

### DESeq22
#load_DESeq22 = require("DESeq22") ## returns FALSE if package is not installed
#if(load_DESeq22 == FALSE) {
#	stop("DESeq22 package could not be loaded. Please check if the package is installed.")
#}
#sampleGrps = data.frame(colnames(countTable), conditions, row.names = 1)
#DESDataSet = DESeq2DataSetFromMatrix(countData = countTable, colData = sampleGrps, design = ~conditions)
#DESDataSet = DESeq2(DESDataSet)

## perform pair-wise DE comparison
for (g1 in 1:(length(Grps)-1)) { 
	for (g2 in (g1+1):length(Grps)) {
		grp1 = Grps[g1]
		grp2 = Grps[g2]
		grp1_samples = which(conditions==grp1)
		grp2_samples = which(conditions==grp2)
		cat (paste("Group 1 is ", grp1, "Group 2 is ", grp2, "\n")) 
		if (!length(grp1_samples)) {cat (paste("No samples found in ", grp1, "\n")); next
		} else if (!length(grp2_samples)) {cat (paste("No samples found in ", grp2, "\n")); next}

		DE = nbinomTest(countDataSet, grp1, grp2) ### pair-wise DE
		DE = DE[order(DE$pval),]
		#outDir = paste("DESeq2", dataSet, sep="/")
		#outDir = ""
		dir.create(file.path(baseDir), recursive=T, showWarnings = FALSE)
		DEFile = paste("DESeq2_" ,Grps[g1], "_",  Grps[g2], "_DE.txt", sep="")  #abhasi : "_" instead of "_v_"
		write.table(DE, file=paste(baseDir, DEFile , sep="/"), sep="\t", row.names=F)
		DEsig = DE[DE$padj < 0.05 & (DE$foldChange>=1.5 | DE$foldChange<=1/1.5),]
		DEsig = na.omit(DEsig)
		DEsig = DEsig[order(DEsig$pval),]
		DESigFile = paste("DESeq2_" ,Grps[g1], "_",  Grps[g2], "_DESig.txt", sep="")  #abhasi : "_" instead of "_v_"
		write.table(DEsig, file=paste(baseDir, DESigFile , sep="/"), sep="\t", row.names=F)
	}	
}
cat("All done!\n")



