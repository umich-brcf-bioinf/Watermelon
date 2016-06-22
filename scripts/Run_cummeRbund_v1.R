
### Manjusha Pande, mpande@umich.edu, Jan 23, 2013, modified Apr 16, 2013
########################################################
### Usage: Rscript <path to the Run_cummeRbund.R file> baseDir=<dir> gtfFile=<file.gtf> genome=<genome> [grpRepFile=<grpRepFile>]
#####################################################################################################
### Example:
## $ Rscript /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/Run_cummeRbund_v1.R baseDir="/ccmb/home/mpande/CummeRbund/WT_v_SCA3KI" gtfFile="/ccmb/home/mpande/mm9_cuffcmp.combined.gtf" genome="mm9"
#####################################################################################################

args <- commandArgs(T)

for (i in args) {
	arg = strsplit(i,"=",fixed=TRUE)
	assign(arg[[1]][1],arg[[1]][2])	
}

if (!exists("baseDir") | !exists("gtfFile") | !exists("genome") | !exists("cuffDiffDir")) {
	cat (paste ("Incorrect input.\n Usage: Rscript Run_cummeRbund.R baseDir=<dir> gtfFile=<file.gtf> genome=<genome> \n")) 
} else {
	#cat (paste("baseDir is ", baseDir, "\ngtfFile is ", gtfFile, "\ngenome is ", genome, "\n"))
}
if (!exists("grpRepFile")) { cat (paste ("Group replicate file not specified.\n")) }

if (baseDir == ".") {
	baseDir = getwd()
}

setwd(baseDir) 
outDir = "Plots"
dir.create(file.path(getwd(), baseDir, outDir),showWarnings = FALSE)

dataSetPath = strsplit(baseDir, "/")[[1]]
dataSetName = dataSetPath[length(dataSetPath)]

#library(cummeRbund) ## throws error if package is not installed
load_cummeRbund = require("cummeRbund") ## returns FALSE if package is not installed
if(load_cummeRbund == FALSE) {
	stop("cummeRbund package could not be loaded - is it installed?")
}

#tryCatch ({cuff = readCufflinks(genome=genome, gftFile=gtfFile, rebuild=T)
#	}, error = function(e) {cat (paste((e$message), ". Could not read cuffdiff output using the provided gtf file\n"))}) ## incorrectly formatted gtf file 

#if (!exists("cuff")) {
#	cat ("Reading cuffdiff output without using the gtf file\n")
#	cuff = readCufflinks()  # does not work well with rebuild=T 
#}

cuffdiff_dir = cuffDiffDir
cuff = readCufflinks(dir=cuffdiff_dir, genome=genome, gftFile=gtfFile, rebuild=T)
#cuff = readCufflinks(dir=cuffdiff_dir, genome=genome, gftFile=gtfFile)

#cuff = readCufflinks()
#cuff 

samples = samples(cuff)$sample_name  ## sample groups
#nSamples = length(samples)
#reps = replicates(cuff)$replicate
genes.features = annotation (genes(cuff)) 
features = subset(genes.features, select = c(gene_id, gene_short_name))

## per gene raw and normalized counts
genes.readgroup = read.table(paste(cuffdiff_dir, "genes.read_group_tracking", sep="/"), sep="\t", header=T)
reps = unique(genes.readgroup[,c("condition", "replicate")])
repNames = c()
for (i in 1:length(reps[,1])) { repNames = c(repNames, paste(reps[i,"condition"],reps[i,"replicate"], sep="_"))}

genes.rawFrags = genes.readgroup[genes.readgroup[,"condition"]== reps[1,"condition"] & genes.readgroup[,"replicate"]== reps[1,"replicate"], c("tracking_id","raw_frags")]
colnames(genes.rawFrags)[2] <- paste("RawFrags",repNames[1], sep=" ")

for(i in 2:length(repNames)) {
    genes.rawFrags = merge(genes.rawFrags, genes.readgroup[genes.readgroup[,"condition"]== reps[i,"condition"] & genes.readgroup[,"replicate"]== reps[i,"replicate"], c("tracking_id","raw_frags")],
                          by.x=1, by.y=1, all.x=TRUE, all.y=TRUE)
    colnames(genes.rawFrags)[ncol(genes.rawFrags)] = paste("RawFrags", repNames[i], sep=" ")
}
repRawCountMatrix = merge(features, genes.rawFrags, by.x = "gene_id", by.y = "tracking_id")

sampleCols = colnames(repRawCountMatrix)[3:ncol(repRawCountMatrix)]
repRawCountMatrix = repRawCountMatrix[c(colnames(repRawCountMatrix[1:2]), sort(sampleCols))]
colnames(repRawCountMatrix) = sub("gene_id", "Gene Id", colnames(repRawCountMatrix))
colnames(repRawCountMatrix) = sub("gene_short_name", "Gene Name", colnames(repRawCountMatrix))

genes.repCount.matrix = repCountMatrix(genes(cuff)) ## scaled counts per gene
gene_id = rownames(genes.repCount.matrix)
genes.repCount.matrix = cbind(gene_id, genes.repCount.matrix)
repCountMatrix = merge(features, genes.repCount.matrix, by.x = "gene_id", by.y = "gene_id") ## gene_id, gene_short_name, per sample count
nCountsFileName = paste(dataSetName, "repScaledCounts.txt", sep = "_")
write.table (repCountMatrix, file=paste(getwd(), nCountsFileName,  sep ="/"), sep="\t", quote = FALSE, row.names=F)

# isoforms.repCount.matrix = repCountMatrix(isoforms(cuff)) ## counts per isoform

## per gene FPKMs
genes.repFpkm.matrix = repFpkmMatrix(genes(cuff)) 
gene_id = rownames(genes.repFpkm.matrix)
genes.repFpkm.matrix = cbind(gene_id, genes.repFpkm.matrix)
repFpkmMatrix = merge(features, genes.repFpkm.matrix, by.x = "gene_id", by.y = "gene_id")
fpkmFileName = paste(dataSetName, "repFpkms.txt", sep = "_")
write.table (repFpkmMatrix, file=paste(getwd(), fpkmFileName,  sep ="/"), sep="\t", quote = FALSE, row.names=F)

## get FPKMs for DE genes
DEGFile = paste(dataSetName, "DEG_names.txt", sep = "_") ## file with gene names and FC
if (file.exists(DEGFile)) {
	tryCatch ({DEGs = read.table(paste(getwd(), "../11-gene_annotation", DEGFile, sep="/")) ## data frame with gene names and FC
		DEGs = DEGs[order(-DEGs$V2), ] ## descending order by FC
		DEGsAll = DEGs$V1
	
		### select top <= 100 and bottom <= 100 genes for heatmap
		DEGsUp = DEGs[DEGs$V2>=1, ] ## up-regulated genes
		DEGsDown = DEGs[DEGs$V2<1, ] ## down-regulated genes
		DEGsDown = DEGsDown[order(DEGsDown$V2),] ## ascending order by FC for down-regulated genes 

		if (length(DEGsUp$V1)>100) {  ## select top 100 up-regulated
		DEGsUp = DEGsUp[1:100,]
		}
		if (length(DEGsDown$V1)>100) { ## select top 100 down-regulated
			DEGsDown = DEGsDown[1:100,]
		}
		nUp = length(DEGsUp$V1)
		nDown = length(DEGsDown$V1)
		geneNames = rbind(DEGsUp, DEGsDown)
		geneNames = as.vector(geneNames$V1) ## vector of gene names
		myFpkms = repFpkmMatrix[repFpkmMatrix$gene_short_name%in%geneNames,] ### repFpkmsMatrix (data frame) for genes of interest
		### Get a matrix from the dataframe
		rnames = myFpkms[,1] # assign labels in column 1 to "rnames"
		mat = data.matrix(myFpkms[,3:ncol(myFpkms)])  # transform column 2-5 into a matrix
		rownames(mat) = rnames                  # assign row names 
		transform = function(x) log10(x+1)
		mat_transform =  apply(mat, 2, transform)

		### Generate heatmaps with dendrograms for top <= 200 DE genes 
		if (!require("gplots")) {
  			install.packages("gplots", dependencies = TRUE)
  			library(gplots)
   		}
		if (!require("RColorBrewer")) {
   			install.packages("RColorBrewer", dependencies = TRUE)
   			library(RColorBrewer)
  		 }
		file = paste(dataSetName, paste("DEGs_", paste("up", nUp, sep = ""), paste("down" ,nDown, sep = ""), "_fpkm_heatmap.png", sep = ""), sep = "_")
		#col_breaks = c(seq(min(mat_transform),1,length=1000), seq(1.01,2,length=1000), seq(2.01,max(mat_transform),length=1000)) ## non-overlapping breaks
		png(paste(outDir,file,sep="/"),  height=10+2/3, width=10+2/3, units="in", res=1200, pointsize = 6)
		#heatmap.2(mat_transform, trace="none", margins =c(12,9), col=rev(heat.colors(2999)), breaks=col_breaks, dendrogram="both", keysize=0.5, density.info="none") ##
		#heatmap.2(mat_transform, trace="none", margins =c(12,9), col=bluered(2999), breaks=col_breaks, dendrogram="both", keysize=0.5, density.info="none")
		heatmap.2(mat_transform, trace="none", margins =c(12,9), col=bluered(2999), dendrogram="both", keysize=0.5, density.info="none")
		#heatmap.2(mat_transform, trace="none", margins =c(12,9), col=bluered(2999), keysize=0.5, density.info="none")
		legend("bottomleft", "value=log10(FPKM+1)")
		dev.off() 

		### use all DEGs for heatmap
		DEGNames = as.vector(DEGsAll) ## vector of DE gene names
		myDEGFpkms = repFpkmMatrix[repFpkmMatrix$gene_short_name%in%DEGNames,] ### repFpkmsMatrix (data frame) for all DE genes
		#DEGFpkmFileName = paste(dataSetName, "DEGFpkms.txt", sep = "_")
		#write.table (myDEGFpkms, file=paste(baseDir, DEGFpkmFileName,  sep ="/"), sep="\t", quote = FALSE, row.names=F)
		rnames_DEG = myDEGFpkms[,1] # assign labels in column 1 to "rnames"
		mat_DEG = data.matrix(myDEGFpkms[,3:ncol(myDEGFpkms)])  # transform column 2-5 into a matrix
		rownames(mat_DEG) = rnames_DEG                  # assign row names 
		mat_DEG_transform =  apply(mat_DEG, 2, transform)
	
		file = paste(dataSetName, "DEGsAll_fpkm_heatmap.png", sep = "_")
		#col_breaks = c(seq(min(mat_DEG_transform),1,length=1000), seq(1.01,2,length=1000), seq(2.01,max(mat_DEG_transform),length=1000)) ## non-overlapping breaks
		png(paste(outDir,file,sep="/"),  height=10+2/3, width=10+2/3, units="in", res=1200, pointsize = 4)
		#heatmap.2(mat_DEG_transform, trace="none", margins =c(12,9), col=rev(heat.colors(2999)), breaks=col_breaks, dendrogram="both", keysize=0.5, density.info="none") ##
		#heatmap.2(mat_DEG_transform, trace="none", margins =c(12,9), col=bluered(2999), breaks=col_breaks, dendrogram="both", keysize=0.5, density.info="none")
		heatmap.2(mat_DEG_transform, trace="none", margins =c(12,9), col=bluered(2999), dendrogram="both", keysize=0.5, density.info="none")
		legend("bottomleft", "value=log10(FPKM+1)")
		dev.off()

		### stretch the image to make gene names legible
		if (length(DEGNames) > 2000) {
			file = paste(dataSetName, "DEGsAll_fpkm_heatmap_stretched.png", sep = "_")
			#col_breaks = c(seq(min(mat_DEG_transform),1,length=1000), seq(1.01,2,length=1000), seq(2.01,max(mat_DEG_transform),length=1000)) ## non-overlapping breaks
			png(paste(outDir,file,sep="/"),  height=30+2/3, width=10+2/3, units="in", res=1200, pointsize = 1)
			#heatmap.2(mat_DEG_transform, trace="none", margins =c(12,9), col=rev(heat.colors(2999)), breaks=col_breaks, dendrogram="both", keysize=0.5, density.info="none") ##
			#heatmap.2(mat_DEG_transform, trace="none", margins =c(12,9), col=bluered(2999), breaks=col_breaks, dendrogram="both", keysize=0.5, density.info="none")
			heatmap.2(mat_DEG_transform, trace="none", margins =c(12,9), col=bluered(2999), dendrogram="both", keysize=0.5, density.info="none")
			legend("bottomleft", "value=log10(FPKM+1)")
			dev.off() 
		}
	}, error = function(e) {cat (paste ("Error occurred in read.table(DEGFile).", e, "\n"))})
}

## Plot correlation matrix 
library(lattice)
rgb.palette <- colorRampPalette(c("blue", "yellow"), space = "rgb")

count_cor = cor(data.matrix(repRawCountMatrix[,3:ncol(repRawCountMatrix)])) ## pearson correlation

if (!any(is.na(count_cor))) {
    count_corFileName = paste(dataSetName, "rawCounts_correlation.txt", sep = "_")
    write(paste(c("", colnames(count_cor)),collapse="\t"), count_corFileName)
    write.table(count_cor, paste(getwd(), count_corFileName, sep="/"), col.names=F, quote = FALSE, sep= "\t", append=T)

    file = paste(dataSetName, "rawCounts_correlation.pdf", sep = "_")
    print(paste(getwd(), outDir, file, sep="/"))
    pdf(paste(getwd(), outDir, file, sep="/"))
    levelplot(count_cor, main="Raw Counts Correlation Matrix", xlab="", ylab="", scales=list(x=list(rot=45)), col.regions=rgb.palette(120), cuts=100, at=seq(min(count_cor), 1, (max(count_cor) - min (count_cor))/10))
    dev.off()
}

fpkm_cor = cor(repFpkmMatrix(genes(cuff)))

if (!any(is.na(fpkm_cor))) {
    fpkm_corFileName = paste(dataSetName, "fpkm_correlation.txt", sep = "_")
    write(paste(c("", colnames(fpkm_cor)),collapse="\t"), fpkm_corFileName)
    write.table(fpkm_cor, paste(getwd(), fpkm_corFileName, sep="/"), col.names=F, quote = FALSE, sep= "\t", append=T)

    #fpkm_cor = cor(repFpkmMatrix(genes(cuff)))
    file = paste(dataSetName, "fpkm_correlation.pdf", sep = "_")
    pdf(paste(getwd(), outDir, file, sep="/"))
    levelplot(fpkm_cor, main="FPKM Correlation Matrix", xlab="", ylab="", scales=list(x=list(rot=45)), col.regions=rgb.palette(120), cuts=100, at=seq(min(fpkm_cor), 1, (max(fpkm_cor) - min (fpkm_cor))/10))
    dev.off()
}

###################################################################################################
## replace colnames in raw counts matrix with sample Ids
if (exists("grpRepFile") && file.exists(grpRepFile)) {
	tryCatch ({grpRepSamps = read.table(grpRepFile, sep="\t")
		repSamps = c()
		for (i in 1:length(grpRepSamps[,1])) { 
			tmp = c()
			tmp = gsub ("replicate", grpRepSamps[i,]$V1, grpRepSamps[i,]$V2)
			tmp = unlist (strsplit(tmp, ", ",fixed=T))
			repSamps = c(repSamps, tmp)
		}

		for(i in 1:length(colnames(repRawCountMatrix))) {
			for(j in 1:length(repSamps)) {
				tmpgrpRep = sub("RawFrags ", "",  colnames(repRawCountMatrix)[i])
				match = grep(tmpgrpRep, repSamps[j])
				tmp = "" 
				if (length(match) != 0) {tmp = gsub (paste(tmpgrpRep, ": ", sep=""), "", repSamps[j])
				 #colnames(repRawCountMatrix)[i] = paste("RawFrags ", tmp, sep = "")
					colnames(repRawCountMatrix)[i] = tmp 
				}
			}
		}
	}, error = function(e) {cat (paste ("Error occurred in read.table(grpRepFile).", e, "\n"))})
}
###################################################################################################
rawCountsFileName = paste(dataSetName, "repRawCounts.txt", sep = "_")
write.table (repRawCountMatrix, file=paste(getwd(), rawCountsFileName,  sep ="/"), sep="\t", quote = FALSE, row.names=F)

### Generate cummeRbund plots 
#cat (paste ("Generating Plots.\n")) 

### Estimated overdispersion for each sample group
cat (paste ("Dispersion Plots.\n")) 
tryCatch ({
	disp = dispersionPlot(genes(cuff)) 	
	file = paste(dataSetName, "disp.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/")) 
	print (disp)	
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in Dispersion Plots.", e, "\n"))})
graphics.off()

### Squared coefficient of variation 
cat (paste ("Fpkm SCV Plots.\n")) 
tryCatch ({
	genes.scv = fpkmSCVPlot(genes(cuff))
	file = paste(dataSetName, "genes_scv.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (genes.scv)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in genes SCV Plots.", e, "\n"))})
graphics.off()

tryCatch ({
	isoforms.scv = fpkmSCVPlot(isoforms(cuff))
	file = paste(dataSetName, "isoform_scv.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (isoforms.scv) 
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in isforms SCV Plots.", e, "\n"))})
graphics.off()

### Distribution of FPKM scores across samples
cat (paste ("Density plot.\n")) 
tryCatch ({
	dens = csDensity(genes(cuff))
	file = paste(dataSetName, "dens.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print(dens)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in Density Plots.", e, "\n"))})
graphics.off()

tryCatch ({
	densRep = csDensity(genes(cuff), replicates=T)
	file = paste(dataSetName, "densRep.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (densRep)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in replicates Density Plots.", e, "\n"))})
graphics.off()

### boxplots of FPKM distribution
cat (paste ("Box Plots.\n")) 
tryCatch ({
	b = csBoxplot(genes(cuff))
	file = paste(dataSetName, "boxplot.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (b)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in Boxplots.", e, "\n"))})
graphics.off()

tryCatch ({
	brep = csBoxplot(genes(cuff), replicates=T)
	file = paste(dataSetName, "boxplotRep.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (brep)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in replicates Boxplots.", e, "\n"))})
graphics.off()

### pair-wise scatterplots (FPKMs)
cat (paste ("Scatter Plots.\n"))
tryCatch ({
	s = csScatterMatrix(genes(cuff)) 
	file = paste(dataSetName, "scatterplots.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file, sep="/"))
	print (s)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in Scatterplots.", e, "\n"))})
graphics.off()

### dendrogram based on the JS (Jensen-Shannon divergence) distance
cat (paste ("Dendrograms.\n"))
tryCatch ({
	dend = csDendro(genes(cuff))
	file = paste(dataSetName, "dend.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"), width=10, height=12) ## stretch vertically to accomodate long labels
	plot(dend)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in Dendrograms.", e, "\n"))})
graphics.off()

tryCatch ({
	dend.rep = csDendro(genes(cuff), replicates=T)
	file = paste(dataSetName, "dendRep.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"), width=10, height=12) ## stretch vertically to accomodate long labels
	plot(dend.rep)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in replicates Dendrograms.", e, "\n"))})
graphics.off()


### MA plot: log intensity ratio vs. average log intensity 
### MA plot (FPKM)
cat (paste ("MA Plots.\n"))
tryCatch ({
	m = MAplot(genes(cuff), samples[1], samples[2])
	file = paste(dataSetName, "MA.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (m)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in MA plots.", e, "\n"))})
graphics.off()

### MA plot (normalized counts)
tryCatch ({
	mCount = MAplot(genes(cuff), samples[1], samples[2], useCount=T)
	file = paste(dataSetName, "MA_count.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (mCount)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in normalized counts MA plots.", e, "\n"))})
graphics.off()

### Vocano plots fold change vs. significance
cat (paste ("Volcano Plots.\n"))
tryCatch ({
	v = csVolcanoMatrix(genes(cuff))
	#abline (v = c(-1,1))
	file = paste(dataSetName, "volcano.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (v)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in Volcano plots.", e, "\n"))})
graphics.off()

## Principal Component Analysis (PCA)
cat (paste ("PCA Plots.\n"))
tryCatch ({
	PCA = PCAplot(genes(cuff),"PC1","PC2")
	file = paste(dataSetName, "PCA.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (PCA)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in PCA plots.", e, "\n"))})
graphics.off()

tryCatch ({
	PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
	file = paste(dataSetName, "PCARep.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (PCA.rep)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in replicates PCA plots.", e, "\n"))})
graphics.off()

## Multi-dimensional Scaling (MDS)
cat (paste ("MDS Plots.\n"))
tryCatch ({
	MDS = MDSplot(genes(cuff))
	file = paste(dataSetName, "MDS.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (MDS)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in MDS plots.", e, "\n"))})
graphics.off()

tryCatch ({
	MDS.rep = MDSplot(genes(cuff),replicates=T)
	file = paste(dataSetName, "MDSRep.pdf", sep = "_")
	pdf(paste(getwd(),outDir,file,sep="/"))
	print (MDS.rep)
	dev.off()
}, error = function(e) {cat (paste ("Error occurred in replicates MDS plots.", e, "\n"))})
graphics.off()


