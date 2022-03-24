#!/usr/bin/env Rscript

#This script uses bioMart to retrieve all records for a given dataset - e.g. hsapiens_gene_ensembl
#And outputs an annotation file mapping ENSEMBL ids to entrez IDs, gene symbols, and descriptions

lib.vector <- c("optparse", "plyr", "rtracklayer", "biomaRt")
libs.loaded <- suppressMessages(lapply(lib.vector, library, character.only=T, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-o", "--outfile"), action="store", default=NA, type='character',
              help="Name of output file"),
  make_option(c("-d", "--dataset"), action="store", default=NA, type='character',
              help="ENSEMBL dataset e.g. hsapiens_gene_ensembl"),
  make_option(c("-m", "--mart"), action="store", default="ENSEMBL_MART_ENSEMBL",
              help="Biomart to use [default %default]"),
  make_option(c("-H", "--host"), action="store", default="www.ensembl.org",
              help="Biomart host to use [default %default]"),
  make_option(c("-a", "--attributes"), action="store", default="ensembl_gene_id,entrezgene_id,external_gene_name,description",
              help="Comma-separated list of attributes to gather from biomaRt - must be equivalent to defaults - [default %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

mart = useMart(opt$mart, host=opt$host) #Change host here to use older version
species.mart<- useDataset(opt$dataset, mart=mart)

#Keeping track of current version used
currentMarts <- listMarts(host=opt$host)
thisMart <- currentMarts$version[currentMarts$biomart == opt$mart]
message(paste0("Using ", thisMart))
#Ensembl Genes 105

#Gets attributes using ensembl gene id from biomaRt
attr.cols <- trimws(unlist(strsplit(opt$attributes, ",")))
attr<- getBM(attributes=attr.cols, mart=species.mart)

rename_col_by_name <- function (df, old, new) {
  #https://stackoverflow.com/a/16490387
  names(df)[names(df) == old] <- new
  return(df)
}

#Standardize column name of unique ID to 'gene_id' - This corresponds to the unique ID used in GTF
attr <- rename_col_by_name(attr, attr.cols[1], 'gene_id')
#Standardize other column names
attr <- rename_col_by_name(attr, attr.cols[2], 'entrezgene_id')
attr <- rename_col_by_name(attr, attr.cols[3], 'external_gene_name')
attr <- rename_col_by_name(attr, attr.cols[4], 'description')

#For queries with multiple entrez gene IDs, return these IDs in a single cell, as a string with commas separating the individual IDs
collapsed.entrezIDs <- ddply(attr, .(gene_id), summarize, entrezgene_id=paste(entrezgene_id, collapse=","))

#Sanity check. Gene symbols should be consistent, without conflicting results
test1 <- ddply(attr, .(gene_id), summarize, external_gene_name=paste(external_gene_name, collapse=","))
test1.split <- strsplit(test1$external_gene_name[grepl(",", test1$external_gene_name)], ",")
test1.result <- unique(sapply(test1.split, function(x) {length(unique(x))})) == c(1)
if(test1.result != TRUE && length(test1.split) != 0) {
  stop("Gene symbols have conflicting entries for at least one ensembl_gene_id")
}

#Sanity check. Descriptions are also consistent, no conflicting results
test2 <- ddply(attr, .(gene_id), summarize, description=paste(description, collapse="@"))
test2.split <- strsplit(test2$description[grepl("@", test2$description)], "@")
test2.result <- unique(sapply(test2.split, function(x) {length(unique(x))})) == c(1)
if(test2.result != TRUE && length(test2.split) != 0) {
  stop("Descriptions have conflicting entries for at least one description")
}

#Since entrez gene IDs are the only ones now with conflicting entries, can merge our collapsed entrezIDs from before with a deduplicated attr table
#Rows of attr with duplicated gene_ids should be identical besides entrezgeneID cols. 
#So just grab the other cols and the first instance of each gene_id in attr
non.entrezID.cols <- c("gene_id", "external_gene_name", "description")
attr <- attr[,non.entrezID.cols]
attr <- attr[!duplicated(attr),]
#Verify that gene_ids are unique within attr table
if(nrow(attr) != length(unique(attr$gene_id))) {
  stop("Could not collapse attribute table due to conflicting entries")
}

#And merge them
mapping_table <- merge(collapsed.entrezIDs, attr, by="gene_id", all.x=T)

write.table(mapping_table, opt$outfile, quote=F, sep="\t", row.names=F)

