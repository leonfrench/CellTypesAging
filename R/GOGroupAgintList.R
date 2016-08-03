load(file="/Users/lfrench/Desktop/results/CellTypesAging/aging_genes_aw_fisher_result.cleaned.v2.RData",verbose=T)

library(ggplot2)
library(org.Hs.eg.db)
library(annotate)
library(GO.db)
library(dplyr)

#filter for all genes with any GO group!
symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data='org.Hs.eg')
aw_result <- filter(aw_result, gene_symbol %in% symbolsInGO)

#taken from https://github.com/sarbal/EGAD/blob/master/EGAD/R/auroc_analytic.R
# by Sara Ballouz
source("./R/AUCFunction.R") #load AUC function

go_object <- as.list(org.Hs.egGO2ALLEGS)

geneGroupResults <- data.frame( id=character(), name=character(), geneCount= numeric(), p=numeric(), stringsAsFactors=F)
resultposition <- 0
Sys.time()
count <- 0

awRanking <- rank(aw_result$agingPValuesWithDirection)

for(goGroupName in names(go_object)) {
  
  goGroup <- go_object[goGroupName]
  geneIDs <- unique(unlist(goGroup, use.names=F))  #discard evidence codes
  genesymbols <- unique(getSYMBOL(geneIDs, data='org.Hs.eg'))
  genesymbols <- intersect(aw_result$gene_symbol, genesymbols) #work within the cell type enriched list

  geneCount <- length(genesymbols)
  #if (geneCount < 40 | geneCount > 310) next();
  #if (geneCount > 60 & geneCount < 140) next();
  #if (geneCount > 160 & geneCount < 290) next();
  if (!(length(genesymbols) > 10 & length(genesymbols) < 200)) next();
  resultposition <- resultposition + 1
  #print(paste("GO group", goGroupName, "position:", resultposition))
  if (resultposition %% 100 == 0) {
    print(resultposition)
  }
  
  labelMask <- aw_result$gene_symbol %in% genesymbols
  
  auroc <- auroc_analytic(awRanking, labelMask)
  
  wilcoxP <- wilcox.test(aw_result[!labelMask,]$agingPValuesWithDirection, aw_result[labelMask, ]$agingPValuesWithDirection )$p.value

  geneGroupResults[resultposition, "geneCount"] <- length(genesymbols)
  geneGroupResults[resultposition, "id"] <- goGroupName
  geneGroupResults[resultposition, "name"] <- Term(goGroupName)
  geneGroupResults[resultposition, "AUROC"] <- auroc
  geneGroupResults[resultposition, "p"] <- wilcoxP

}

Sys.time()

geneGroupResults$p.adjusted <- p.adjust(geneGroupResults$p)
geneGroupResults <- geneGroupResults[order(geneGroupResults$p),]
head(geneGroupResults, n=20)
