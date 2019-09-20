setwd("/nfs/med-bfx-activeprojects/trsaari/sandbox/20190910_gtf_manipulation/")

library(plyr)
library(rtracklayer)
library(biomaRt)

gtf <- import("hg19_noRibo.gtf")

#Uses the Danio mart from ENSEMBL
mart = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org") #Change host here to use older version
humanMart<- useDataset("hsapiens_gene_ensembl",mart=mart)

#Keeping track of current version used 
currentMarts <- listMarts()
currentMarts$version[currentMarts$biomart == "ENSEMBL_MART_ENSEMBL"]
#Ensembl Genes 97

id_list = unique(gtf$gene_id)

#Gets attributes using transcript id from biomaRt
attr<- getBM(filters="external_gene_name", attributes=c("external_gene_name", "entrezgene_id", "description"), values=id_list, mart=humanMart)

#Since gene symbols are used as the gene_id in this gtf, the standard gene_id column will be identical to the external_gene_name column
attr$gene_id <- attr$external_gene_name

#For queries with multiple entrez gene IDs, return these IDs in a single cell, as a string with commas separating the individual IDs
collapsed.entrezIDs <- ddply(attr, .(gene_id), summarize, entrezgene_id=paste(entrezgene_id, collapse=","))

#Sanity check. Gene symbols should be consistent, without conflicting results
test1 <- ddply(attr, .(gene_id), summarize, external_gene_name=paste(external_gene_name, collapse=","))
View(test1[grepl(",", test1$external_gene_name),])
#Another way of testing automagically
test1.split <- strsplit(test1$external_gene_name[grepl(",", test1$external_gene_name)], ",")
unique(sapply(test1.split, function(x) {length(unique(x))})) == c(1)

#Sanity check. See if descriptions are consistent without conflicting results
test2 <- ddply(attr, .(gene_id), summarize, description=paste(description, collapse="@"))
test2.split <- strsplit(test2$description[grepl("@", test2$description)], "@")
unique(sapply(test2.split, function(x) {length(unique(x))})) == c(1)
#Nope - investigate
conflicting.desc.mask <- sapply(test2.split, function(x) {length(unique(x))}) != 1
test2.split[conflicting.desc.mask]
#Looks like only the sources differ - Will just remove this source info and carry on
attr$description <- sub(" \\[Source.*\\]$", "", attr$description)
#And now check
test2 <- ddply(attr, .(gene_id), summarize, description=paste(description, collapse="@"))
test2.split <- strsplit(test2$description[grepl("@", test2$description)], "@")
unique(sapply(test2.split, function(x) {length(unique(x))})) == c(1)

#Since entrez gene IDs are the only ones with conflicting entries, can merge our collapsed entrezIDs from before with the attr table

mapping_table <- merge(collapsed.entrezIDs, attr, on="gene_id", all.x=T)

write.table(mapping_table, "hg19_noRibo_annotation.tsv", quote=F, sep="\t", row.names=F)
