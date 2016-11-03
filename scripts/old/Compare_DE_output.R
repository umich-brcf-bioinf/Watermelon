### Manjusha Pande, mpande@umich.edu, Feb 9, 2015
########################################################
### Usage: Rscript <path to the Compare_DE_output.R file> baseDir=<dir> inFile1=<cuffdiff_output.txt> inFile2=<deseq_output.txt>
#####################################################################################################
### Example:
## $ Rscript /ccmb/BioinfCore/Projects/Manjusha/RNA-seq_pipeline/RSQ_scripts/Compare_DE_output.R baseDir="/ccmb/BioinfCore/Projects/Myers_1176_mgmyers_mpande_RS15/DE_output_comparison/LepR_WT_v_Non_WT" inFile1="/ccmb/BioinfCore/Projects/Myers_1176_mgmyers_mpande_RS15/cuffdiff/LepR_WT_v_Non_WT/LepR_WT_v_Non_WT_gene.foldchange.1.5.txt" inFile2="/ccmb/BioinfCore/Projects/Myers_1176_mgmyers_mpande_RS15/DESeq/DESeq_LepR_WT_v_Non_WT_DE.txt"
#####################################################################################################

args <- commandArgs(T)

for (i in args) {
	arg = strsplit(i,"=",fixed=TRUE)
	assign(arg[[1]][1],arg[[1]][2])	
}

if (!exists("baseDir") | !exists("gtfFile") | !exists("genome")) {
	cat (paste ("Incorrect input.\n Usage: Rscript Run_cummeRbund.R baseDir=<dir> gtfFile=<file.gtf> genome=<genome> \n")) 
}else {
	#cat (paste("baseDir is ", baseDir, "\ngtfFile is ", gtfFile, "\ngenome is ", genome, "\n"))
}

if (baseDir == ".") {
	baseDir = getwd()
}

setwd(baseDir) 
