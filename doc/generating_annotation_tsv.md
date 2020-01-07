# Example - Generating an annotation_tsv

This Rscript was uses BiomaRt to create an annotation_tsv for S. scrofa (pig)


    library(plyr)
    library(biomaRt)

    mart = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org") #Change host here to use older version
    pigMart<- useDataset("sscrofa_gene_ensembl",mart=mart)

    #Keeping track of current version used
    currentMarts <- listMarts()
    currentMarts$version[currentMarts$biomart == "ENSEMBL_MART_ENSEMBL"]
    #Ensembl Genes 97

    #Gets attributes using transcript id from biomaRt
    attr<- getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "external_gene_name", "description"), mart=pigMart)

    rename_col_by_name <- function (df, old, new) {
      #https://stackoverflow.com/a/16490387
      names(df)[names(df) == old] <- new
      return(df)
    }

    #Standardize column name of unique ID to 'gene_id' - This corresponds to the unique ID used in GTF
    attr <- rename_col_by_name(attr, 'ensembl_gene_id', 'gene_id')

    #For queries with multiple entrez gene IDs, return these IDs in a single cell, as a string with commas separating the individual IDs
    collapsed.entrezIDs <- ddply(attr, .(gene_id), summarize, entrezgene_id=paste(entrezgene_id, collapse=","))

    #Sanity check. Gene symbols should be consistent, without conflicting results
    test1 <- ddply(attr, .(gene_id), summarize, external_gene_name=paste(external_gene_name, collapse=","))
    View(test1[grepl(",", test1$entrezgene_id),])
    #Another way of testing automagically
    test1.split <- strsplit(test1$external_gene_name[grepl(",", test1$external_gene_name)], ",")
    unique(sapply(test1.split, function(x) {length(unique(x))})) == c(1)

    #Sanity check. Descriptions are also consistent, no conflicting results
    test2 <- ddply(attr, .(gene_id), summarize, description=paste(description, collapse="@"))
    test2.split <- strsplit(test2$description[grepl("@", test2$description)], "@")
    unique(sapply(test2.split, function(x) {length(unique(x))})) == c(1)

    #Since entrez gene IDs are the only ones with conflicting entries, can merge our collapsed entrezIDs from before with the attr table

    mapping_table <- merge(collapsed.entrezIDs, attr, on="gene_id", all.x=T)

    write.table(mapping_table, "Sscrofa11.1_annotation.tsv", quote=F, sep="\t", row.names=F)
