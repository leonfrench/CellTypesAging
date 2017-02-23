rm(list = ls())

#clear all libraries to prevent conflicts
#from mjaniec
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
detachAllPackages()

library(readr)
library(reshape2)
library(homologene)
library(dplyr)
library(org.Hs.eg.db)
library(annotate)
library(GO.db)
library(gplots) #for heatmap2

level <- "TranscriptomicName"
filterForCoreTranscriptomeTypes <- T #tasic param
targetCellType <- "Sst Cdk6"
goGroupName <- "GO:0007268" #synaptic transmission

#targetCellType <- "Vip Mybpc1"
#goGroupName <-"GO:0007267" #cell-cell signalling

source("./R/TasicProcessingCode.R")

#to know cell type-specific genes
load("/Users/lfrench/Desktop/results/CellTypesAging/server results/AW.Tasic.emperical.iter.10000.TranscriptomicName.1479930034/rpkmByType2xHuman.RData", v=T)

go_object <- as.list(org.Hs.egGO2ALLEGS)

symbolsInGO <- getSYMBOL(unlist(go_object[goGroupName]), data='org.Hs.eg')

humanTargetGenesTypeAndGO <- filter(rpkmByType2xHuman, CellClassID == targetCellType, humanGene %in% unlist(symbolsInGO))
targetMouseGenes <- humanTargetGenesTypeAndGO$geneName

rpkmDataCore <- filter(rpkmDataCore, geneName %in% targetMouseGenes)
#average across the cell types - needs loaded data
rpkmDataCore <- group_by(rpkmDataCore, CellClassID, geneName) #slow
#average across cell types for each gene
rpkmByType <- summarise(rpkmDataCore, log1Expression = mean(log1Expression))

#scale/z-score - do z score after averaging in a cell type
rpkmByType <- group_by(rpkmByType, geneName)
rpkmByType <- mutate(rpkmByType, log1ExpressionZ = (log1Expression - mean(log1Expression)) / sd(log1Expression))

targetByType <- rpkmByType
targetTypeMatrix <- dcast(targetByType, geneName ~ CellClassID, value.var = "log1ExpressionZ")
rownames(targetTypeMatrix) <- targetTypeMatrix$geneName
targetTypeMatrix <- targetTypeMatrix[,-1]
heatmap(as.matrix(targetTypeMatrix),scale="none")
#heatmap.2(as.matrix(targetTypeMatrix),scale="none")
