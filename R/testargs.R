
#clear all libraries to prevent conflicts
#from mjaniec
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
detachAllPackages()

library(doParallel)
library(readr)
library(reshape2)
library(homologene)
library(ROCR)
library(optparse)
library(dplyr)


option_list = list(
  make_option(c("-a", "--ageSource"), type="character", default=NULL, help="AW, Blood or Mistry", metavar="character"),
  make_option(c("-s", "--cellSource"), type="character", default=NULL, help="Tasic or Zeisel", metavar="character"),
  make_option(c("-l", "--cellLevel"), type="character", default=NULL, help="TranscriptomicName (lowest level) or mainClass (highest level)", metavar="character"),
  make_option(c("-i", "--iterations"), type="integer", default=NULL, help=""),
  make_option(c("-c", "--cores"), type="integer", default=NULL, help="number of cores when multithreading")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (interactive()) { #set the variables manually if in Rstudio, for testing
  agingGeneSource <- "AW"
  #agingGeneSource <- "Blood"
  #agingGeneSource <- "Mistry"
  
  cellTypeSource <- "Tasic"
  #cellTypeSource <- "Zeisel"
  
  #level <- "TranscriptomicName"
  level <- "mainClass"
  
  cores <- 4
  iterations <- 4
  
} else if (!is.null(opt$ageSource) & !is.null(opt$cellSource) & !is.null(opt$iterations)) {
  agingGeneSource <- opt$ageSource
  cellTypeSource <- opt$cellSource
  iterations <- opt$iterations
  cores<-opt$cores
  level <- opt$cellLevel
} else {
  print_help(opt_parser)
  stop()
}

registerDoParallel(cores=cores)

#other params to add to options, currently hardcoded
filterForCoreTranscriptomeTypes <- T #tasic param
enrichmentThreshold <- 2

name <- paste(agingGeneSource, cellTypeSource,"emperical", "iter", iterations,level,sep=".")
print(name)