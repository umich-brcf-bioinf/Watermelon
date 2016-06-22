####################################################################################
##
## Processes either gene or isoform based cuffdiff2 output by combining data from
## files _exp.diff, .fpkm_tracking and .read_group_tracking into a single table, 
## including gene annotations
##
## shell> Rscript process_cuffdiff2.r target=[genes|isoforms] significant=[TRUE|FALSE]
##        logFC=[int] fpkm=[int] path=[path to cuffdiff output]
##
## 'target' and 'significant' arguments are required
## other arguments have defaults of 0, 0 and ., respectively
##
## Example: shell> Rscript process_cuffdiff2.r target=genes significant=TRUE
##
## Note: to use this function you will have to install the cummeRbund package, 
##       version 1.99.1 or later:
## R> source("http://www.bioconductor.or/biocLite.R")
## R> biocLite()
## R> biocLite("cummeRbund")
## R> biocLite("BiocUpgrade") # if the above grabs an older version of cummeRbund
## R> q()
##
## The first time you run this, cummeRbund will read all info from cuffdiff
## output directory into a database (cuffData.db). This can take sometime. 
## Subsequent runs for the same dataset should be much quicker. 
## When a new version of cummeRbund is installed you'll probably have to delete 
## the cuffData.db file so that a compatible database can be built. 
##
## Ana / Oct 1, 2012, edited Dec 13 to add RawFrags

getArgs <- function(set) {

  ## [ Function shared by Barry ]
  ## Parse command line options (for Rscript IO handeling)
  ## If there is an equals sign (=) in the input args then it
  ## is parsed and the first part (before the = sign) treated
  ## as the variable name, with the second (after the = sign )
  ## taken as the value of the argument.
  ## If there is no = sign then the arg value is assigned
  ## logical TRUE
  ##
  ## Note that all values will be characters (i.e. strings).
  ## so YOU must sort them out afterwords
  ##
  ## E.g.
  ##
  ##  #!/usr/bin/env Rscript
  ##  got <- getArgs( c("dir","id","nfile") )
  ##
  ##  if(!got["dir"])
  ##    stop("Command line argument 'dir=?' is required")
  ##  if(!got["id"]) {
  ##    warning("Assuming 'id' is the same as basename(dir)")
  ##    id <- basename(dir) }
  ##  if(!got["nfile"])
  ##    nfile <- NULL
  ##  ...


  passed.args <- commandArgs(TRUE)
#  cat( paste("   Called with", length(passed.args),"arguments:\n\t-  ",
#             paste(passed.args, collapse=",\n\t-   ")),"\n\n")

  for (e in passed.args) {
    ta = strsplit(e,"=",fixed=TRUE)
    if(ta[[1]][1] %in% set) {
      if(! is.na(ta[[1]][2])) {
        assign(ta[[1]][1], ta[[1]][2], env = .GlobalEnv)
        cat("   Assigned ",ta[[1]][1]," the value of:",ta[[1]][2],"\n")
      } else {
        assign(ta[[1]][1],TRUE, env = .GlobalEnv)
        cat("   Assigned ",ta[[1]][1]," the value of: TRUE\n")
      }
    } else {
      warning(paste(" Command line argument:", ta[[1]][1],
                    "is undefined and will be ignored"))
    }
  }
  sapply(set, exists)
}

process.cuffdiff <- function(path, attribute) {

  ## process.cuffdiff, by Ana Oct 1, 2012
  ## function to process cuffdiff2 output using R package cummeRbund
  ## can process either gene or isoform based output
  ## must suply the path to the cuffdiff output directory (this can be .)
  ## combines data from files _exp.diff, .fpkm_tracking and .read_group_tracking
  ## returns a single combined table, including gene annotations

  load.cummeRbund <- require("cummeRbund")
  if(load.cummeRbund==FALSE) {
    stop("cummeRbund package could not be loaded - is it installed?")
  }

  # first time this is run it reads in all cuffdiff input into a database
  cuff <- readCufflinks(path)

  # retrieve  data for the chosen attribute (genes or isoforms)
  all             <- get(attribute)(cuff) # ie, genes(cuff) or isoforms(cuff)
  
  # retrieve annotations
  all.features    <- annotation(all)
  
  # retrieve and reformat differential expression info
  all.diffData    <- diffData(all) 
  fpkm.cols <- grep("value_", colnames(all.diffData))
  colnames(all.diffData)[fpkm.cols] <- paste("FPKM", colnames(all.diffData)[fpkm.cols],
                                             sep="_")
  
  # retrieve and reformat fpkm per replicate, rather than sample
  all.repFpkmMat  <- repFpkmMatrix(all)
  colnames(all.repFpkmMat) <- paste("FPKM", colnames(all.repFpkmMat), sep="_")

  # retrieve and reformat external counts (ie normalized) per replicate, rather than sample
  all.readgroup   <- read.table(paste(attribute, "read_group_tracking", sep="."),
                                sep="\t", header=TRUE) # counts per replicate
  reps <- unique(all.readgroup[,c("condition", "replicate")])
  reps.names <- apply(reps, 1, paste, collapse="_")
  all.externalScaledFrags <- all.readgroup[all.readgroup[,"condition"]==
                                           reps[1,"condition"] &
                                           all.readgroup[,"replicate"]==
                                           reps[1,"replicate"],
                                           c("tracking_id","external_scaled_frags")]
  colnames(all.externalScaledFrags)[2] <- paste("ExternalScaledFrags",reps.names[1],
                                                sep="_")
  
  for(i in 2:length(reps.names)) {
    all.externalScaledFrags <- merge(all.externalScaledFrags,
                                     all.readgroup[all.readgroup[,"condition"]==
                                                   reps[i,"condition"] &
                                                   all.readgroup[,"replicate"]==
                                                   reps[i,"replicate"],
                                                   c("tracking_id","external_scaled_frags")],
                                     by.x=1, by.y=1, all.x=TRUE, all.y=TRUE)
    colnames(all.externalScaledFrags)[ncol(all.externalScaledFrags)] <-
      paste("ExternalScaledFrags", reps.names[i], sep="_")
  }

  # retrieve and reformat raw counts per replicate, rather than sample
  all.rawFrags <- all.readgroup[all.readgroup[,"condition"]==
                                reps[1,"condition"] &
                                all.readgroup[,"replicate"]==
                                reps[1,"replicate"],
                                c("tracking_id","raw_frags")]
  colnames(all.rawFrags)[2] <- paste("RawFrags",reps.names[1],
                                     sep="_")
  
  for(i in 2:length(reps.names)) {
    all.rawFrags <- merge(all.rawFrags,
                          all.readgroup[all.readgroup[,"condition"]==
                                        reps[i,"condition"] &
                                        all.readgroup[,"replicate"]==
                                        reps[i,"replicate"],
                                        c("tracking_id","raw_frags")],
                          by.x=1, by.y=1, all.x=TRUE, all.y=TRUE)
    colnames(all.rawFrags)[ncol(all.rawFrags)] <-
      paste("RawFrags", reps.names[i], sep="_")
  }


  
  # combine and return differential expression, replicate fpkm, replicate counts and annotations
  all.data1 <- merge(all.diffData, all.repFpkmMat, by.x=1, by.y="row.names",
                     all.x=TRUE, all.y=TRUE)
  all.data2 <- merge(all.data1, all.externalScaledFrags, by.x=1, by.y=1,
                     all.x=TRUE, all.y=TRUE)
  all.data3 <- merge(all.data2, all.rawFrags, by.x=1, by.y=1,
                     all.x=TRUE, all.y=TRUE)
  all.data.annot <- merge(all.data3, all.features, by.x=1, by.y=1,
                          all.x=TRUE, all.y=FALSE)
  return(all.data.annot)
}

#### read and error check command line 
args <- getArgs( c("target","significant","logFC", "fpkm", "path") )
if(sum(args)<2) {
  stop(paste("Command line arguments 'target' and 'significant' are required\n",
             "shell> Rscript process_cuffdiff2.r target=[genes|isoforms] ",
             "significant=[TRUE|FALSE] logFC=[int] fpkm=[int] path=[path]", sep=""))
}
if(!args["target"]) {
  stop(paste("Command line argument 'target=' is required. Should be 'genes' or 'isoforms'\n",
             "shell> Rscript process_cuffdiff2.r target=[genes|isoforms] ",
             " significant=[TRUE|FALSE] logFC=[int] fpkm=[int] path=[path] ", sep=""))
}
if(!args["significant"]) {
  stop(paste("Command line argument 'significant=' is required. Should be 'TRUE' or 'FALSE'\n",
             "shell> Rscript process_cuffdiff2.r target=[genes|isoforms] ",
             " significant=[TRUE|FALSE] logFC=[int] fpkm=[int] path=[path]", sep=""))
}
if(!args["logFC"]) {
  logFC <- 0
  warning(paste("Assuming 'logFC' cutoff is ", logFC, ", ie no cutoff"), sep="")
}
if(!args["fpkm"]) {
  fpkm <- 0
  warning(paste("Assuming 'fpkm' cutoff is ", fpkm, ", ie no cutoff"), sep="")
}
if(!args["path"]) {
  path <- getwd()
  warning(paste("Assuming 'path' is current working directory: ", path), sep="")
}

#### process data according to command line arguments

## 1. check target is either genes or isoforms and process
if(target=="genes" || target=="isoforms") {
  my.cuffdiff <- process.cuffdiff(path, target)
} else {
  stop("Command line argument 'target=' should be 'genes' or 'isoforms'")
}

## 2. check significant is TRUE or FALSE, grab significant data if required
if(significant == TRUE) {
  my.cuffdiff <- my.cuffdiff[my.cuffdiff[,"significant"]=="yes",] 
} else {
  if(significant !=FALSE) {
    stop("Command line argument 'significant=' should be 'TRUE' or 'FALSE'")
  }
}

## 3. check logFC cutoff is numeric and positive or 0, grab data that passes the cutoff 
logFC <- as.numeric(logFC)
if(!is.na(logFC) && logFC>=0) {
  my.cuffdiff <- my.cuffdiff[abs(my.cuffdiff[,"log2_fold_change"])>=logFC,]
} else {
  stop("'logFC' should be numeric and positive") 
}

## 4. check fpkm cutoff is numeric and positive or 0, grab data that passes the cutoff
fpkm <- as.numeric(fpkm)
if(!is.na(fpkm) && fpkm>=0) {
  my.cuffdiff <- my.cuffdiff[my.cuffdiff[,"FPKM_value_1"]>=fpkm |
                             my.cuffdiff[,"FPKM_value_2"]>=fpkm,]
} else {
  stop("'fpkm' should be numeric and positive")
}

## 5. create output file name using command line arguments
file.name <- paste(paste(target,
                         paste("significant", significant, sep=""),
                         paste("minlogFC", logFC, sep=""),
                         paste("minFPKM", fpkm, sep=""), sep="_"),
                   "txt", sep=".")

## 6. write requested output to output file
write.table(my.cuffdiff, file=file.name, col.names=TRUE, row.names=FALSE,
            sep="\t", quote=FALSE)

