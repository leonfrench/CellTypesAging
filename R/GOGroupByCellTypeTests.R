library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(annotate)
library(GO.db)
library(hash)
library(gplots) #for heatmap2
library(homologene)

#library(reshape2)


#taken from https://github.com/sarbal/EGAD/blob/master/EGAD/R/auroc_analytic.R
# by Sara Ballouz
source("./R/AUCFunction.R") #load AUC function

#show version information
org.Hs.eg.db
go_object <- as.list(org.Hs.egGO2ALLEGS)
length(go_object)
length(go_object$`GO:0007268`)

#filter for all genes with any GO group!
symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data='org.Hs.eg')

#create a hash table that links GO group to gene symbol and name
goTable <- hash()
goNames <- hash()
system.time(
  for(goGroupName in names(go_object)) {
    name <- Term(goGroupName)
    goGroup <- go_object[goGroupName]
    geneIDs <- unique(unlist(goGroup, use.names=F))  #discard evidence codes
    genesymbols <- unique(getSYMBOL(geneIDs, data='org.Hs.eg'))
    goNames[[goGroupName]] <- name
    goTable[[goGroupName]] <- genesymbols
  }
)
print("Done creating GO Hash")
start <- Sys.time()
#prune the GO table
for(goGroupName in keys(goTable)) {
  genesymbols <- goTable[[goGroupName]]
  genesymbols <- intersect(rpkmByType2xHuman$humanGene, genesymbols) 
  if (length(genesymbols) < 11) {
    delete(goGroupName, goTable)
  } else {
    goTable[[goGroupName]] <- genesymbols #put back trimmed list
  }
}
Sys.time() - start

print("Starting gene ontology runs")
start <- Sys.time()
#multiproc
geneGroupResults = foreach(targetCellType=unique(rpkmByType2xHuman$CellClassID), .combine=rbind) %dopar% {  

  thisResult <- data.frame( id=character(), name=character(), geneCount= numeric(), p=numeric(), stringsAsFactors=F)
  resultposition <- 0
  #iterate cell type lists here
  print(targetCellType)
  resultsOfInterest <- subset(rpkmByType2xHuman, CellClassID == targetCellType)
  
  #filter for genes with GO group annotations
  resultsOfInterest <- filter(resultsOfInterest, humanGene %in% symbolsInGO)
  
  for(goGroupName in keys(goTable)) {
    genesymbols <- goTable[[goGroupName]]
    genesymbols <- intersect(resultsOfInterest$humanGene, genesymbols) #work within the cell type enriched list
    #if (length(genesymbols) == 0) next()
    resultsOfInterest <- dplyr::mutate(resultsOfInterest, inGOGroup = humanGene %in% genesymbols)
    
    geneCount <- length(genesymbols)
    #if (geneCount < 40 | geneCount > 310) next();
    #if (geneCount > 60 & geneCount < 140) next();
    #if (geneCount > 160 & geneCount < 290) next();
    if (!(length(genesymbols) > 10 & length(genesymbols) < 200)) next();
    if (all(resultsOfInterest$inGOGroup)) next();
    #print(paste("GO group", goGroupName, "position:", resultposition))
    
    auroc <- auroc_analytic(rank(resultsOfInterest$agingPValuesWithDirection), resultsOfInterest$inGOGroup)
    wilcoxP <- tryCatch({
      wilcox.test(agingPValuesWithDirection ~ inGOGroup, resultsOfInterest)$p.value  #, alternative = "less"? force same direction as whole cell type?
    }, error = function(e) {
      print("ERROR on wilcox")
      print(resultsOfInterest$inGOGroup)
      print(resultsOfInterest)
      print(nrow(resultsOfInterest))
    })
    resultposition <- resultposition + 1
    thisResult[resultposition, "geneCount"] <- length(genesymbols)
    thisResult[resultposition, "id"] <- goGroupName
    thisResult[resultposition, "name"] <- goNames[[goGroupName]]
    thisResult[resultposition, "AUROC"] <- auroc
    thisResult[resultposition, "p"] <- wilcoxP
    
    thisResult[resultposition, "targetCellType"] <- targetCellType
    
    auroc <- auroc_analytic(rank(resultsOfInterest$log1ExpressionZ), resultsOfInterest$inGOGroup)
    wilcoxP <- wilcox.test(log1ExpressionZ ~ inGOGroup, resultsOfInterest)$p.value
    thisResult[resultposition, "spec.AUROC"] <- auroc
    thisResult[resultposition, "spec.p"] <- wilcoxP
  }
  thisResult
}
print("Done runs for GO")
Sys.time() - start

geneGroupResults$uniqueStats <- paste(geneGroupResults$geneCount, geneGroupResults$AUROC, geneGroupResults$spec.AUROC, geneGroupResults$targetCellType)
geneGroupResults <- tbl_df(geneGroupResults)
geneGroupResults <- dplyr::group_by(geneGroupResults, uniqueStats)
geneGroupResults <- dplyr::summarize(geneGroupResults,  geneCount = first(geneCount), p=first(p), AUROC=first(AUROC), targetCellType=first(targetCellType), spec.AUROC= first(spec.AUROC), spec.p=first(spec.p),
                                           synonymIDs=toString(id), synonyms=toString(name), combinedGroups = n(), name=first(name),id=first(id))

geneGroupResults <- geneGroupResults %>% mutate(synonymIDs = if_else(combinedGroups ==1, "", synonymIDs)) %>% mutate(synonyms = if_else(combinedGroups ==1, "", synonyms)) 


geneGroupResults$p.adjusted <- p.adjust(geneGroupResults$p)
geneGroupResults <- geneGroupResults[order(geneGroupResults$p),]
geneGroupResults <- dplyr::select(geneGroupResults, id, name, geneCount, AUROC, p, p.adjusted, targetCellType, spec.AUROC, spec.p, synonyms, synonymIDs, combinedGroups)
print(geneGroupResults, n=50)
head(subset(geneGroupResults,targetCellType!="CelltypeNonSpecific"), n=20)

#output folder from previous run (mainCode.R)
#save(file=paste0(outputFolder,"rpkmByType2xHuman.RData"),rpkmByType2xHuman)
#save(file=paste0(outputFolder,"aging_genes_aw_fisher_data.RData"),aw_result_mouseFilter)
write.csv(geneGroupResults, file=paste0(outputFolder,"CellbyGOEnrichment.csv"),row.names=F)

allCellByGO <- geneGroupResults
allCellByGO <- allCellByGO %>% group_by(targetCellType) %>% arrange(p) %>% summarize(name.1st=first(name), AUROC.1st = first(AUROC), p.adjusted.1st = first(p.adjusted),name.2nd=nth(name,2), AUROC.2nd = nth(AUROC,2), p.adjusted.2nd = nth(p.adjusted,2) ,name.3rd=nth(name,3), AUROC.3rd = nth(AUROC,3), p.adjusted.3rd = nth(p.adjusted,3)) %>% arrange(p.adjusted.1st)
allCellByGO <- allCellByGO %>% mutate_each(funs(signif(.,2)), starts_with("AUROC")) %>% mutate_each(funs(signif(.,2)), starts_with("p.adjusted")) 
write.csv(allCellByGO, paste0(outputFolder, "CellbyGOEnrichmentSummary.csv"),row.names=F)


if (interactive()) { #set the variables manually if in Rstudio, for testing
  
  printAUC <- function(target, goName) {
    resultsOfInterestCellType <- subset(rpkmByType2xHuman, CellClassID == target)
    resultsOfInterestCellType <- dplyr::filter(resultsOfInterestCellType, humanGene %in% symbolsInGO)
    #filter for genes with GO group annotations
    genesymbols <- goTable[[goName]]
    resultsOfInterestCellTypeAndGO <- dplyr::filter(resultsOfInterestCellType, humanGene %in% genesymbols)
    
    (auroc <- auroc_analytic(rank(resultsOfInterestCellType$agingPValuesWithDirection), resultsOfInterestCellType$humanGene %in% resultsOfInterestCellTypeAndGO$humanGene ))
    (wilcoxP <- wilcox.test(resultsOfInterestCellType$agingPValuesWithDirection[resultsOfInterestCellType$humanGene %in% resultsOfInterestCellTypeAndGO$humanGene],resultsOfInterestCellType$agingPValuesWithDirection[!(resultsOfInterestCellType$humanGene %in% resultsOfInterestCellTypeAndGO$humanGene)])$p.value)  #, alternative = "less"? force same direction as whole cell type?
    print(goNames[[goName]])
    print(target)
    print(paste("genes:", nrow(resultsOfInterestCellTypeAndGO)))
    print(paste("AUROC:", signif(auroc,2)))
    print(paste("p:", signif(wilcoxP,3)))
  }

  #load rpkm data - tasic
  load("/Users/lfrench/Desktop/results/CellTypesAging/server results/AW.Tasic.emperical.iter.10000.TranscriptomicName.1479930034/rpkmByType2xHuman.RData", v=T)
  
  #load rpkm data - tasic+Mistry
  load("/Users/lfrench/Desktop/results/CellTypesAging/server results/Mistry.Tasic.emperical.iter.10000.TranscriptomicName.1479930072/rpkmByType2xHuman.RData", v=T)

  #load rpkm data - tasic+Blood
  load("/Users/lfrench/Desktop/results/CellTypesAging/server results/Blood.Tasic.emperical.iter.10000.TranscriptomicName.1479930030/rpkmByType2xHuman.RData", v=T)
  
  
  #get genes for a specific group and cell type
  targetCellType <- "Sst Cdk6"
  goGroupName <- "GO:0007268" #synaptic transmission
  printAUC(targetCellType, goGroupName)
  goGroupName <-"GO:0007267" #cell-cell signalling
  printAUC(targetCellType, goGroupName)
  
  targetCellType <- "Vip Mybpc1"
  goGroupName <-"GO:0007267" #cell-cell signalling
  printAUC(targetCellType, goGroupName)

  ##############Zeisel ####################################  
  load("/Users/lfrench/Desktop/results/CellTypesAging/results/AW.Zeisel.emperical.iter.1.TranscriptomicName.1479938336/rpkmByType2xHuman.RData")
  load("/Users/lfrench/Desktop/results/CellTypesAging/results/Mistry.Zeisel.emperical.iter.1.TranscriptomicName.1479938797/rpkmByType2xHuman.RData")
  load("/Users/lfrench/Desktop/results/CellTypesAging/results/Blood.Zeisel.emperical.iter.1.TranscriptomicName.1479939252/rpkmByType2xHuman.RData")
  
  targetCellType <- "Int2"
  goGroupName <- "GO:0007268" #synaptic transmission
  printAUC(targetCellType, goGroupName)
  
  targetCellType <- "Int2"
  goGroupName <-"GO:0007267" #cell-cell signalling
  printAUC(targetCellType, goGroupName)

  targetCellType <- "Int10"
  goGroupName <-"GO:0007267" #cell-cell signalling
  printAUC(targetCellType, goGroupName)
  
  print(goNames[[goGroupName]])
  resultsOfInterestCellType <- subset(rpkmByType2xHuman, CellClassID == targetCellType)
  
  #filter for genes with GO group annotations
  genesymbols <- goTable[[goGroupName]]
  resultsOfInterestCellTypeAndGO <- dplyr::filter(resultsOfInterestCellType, humanGene %in% genesymbols)
  
  (auroc <- auroc_analytic(rank(resultsOfInterestCellType$agingPValuesWithDirection), resultsOfInterestCellType$humanGene %in% resultsOfInterestCellTypeAndGO$humanGene ))
  (wilcoxP <- wilcox.test(resultsOfInterestCellType$agingPValuesWithDirection[resultsOfInterestCellType$humanGene %in% resultsOfInterestCellTypeAndGO$humanGene],resultsOfInterestCellType$agingPValuesWithDirection[!(resultsOfInterestCellType$humanGene %in% resultsOfInterestCellTypeAndGO$humanGene)])$p.value)  #, alternative = "less"? force same direction as whole cell type?
  
  geneSubsetTable <- resultsOfInterestCellTypeAndGO %>% ungroup() %>% dplyr::select(humanGene, log1ExpressionZ, agingPValuesWithDirection) %>% arrange(agingPValuesWithDirection)
  write.table(geneSubsetTable, file=paste0(outputFolder, "specGenes.", targetCellType,".", goGroupName,".tsv"),sep="\t")
  #########################
  #write heatmap - old code - see HeatMapForCombination.R
  mouseGenes <- resultsOfInterestCellTypeAndGO$geneName

  #call processing code here
  targetMatrix <- filter(rpkmData, geneName %in% mouseGenes)
  targetMatrix <- as.data.frame(targetMatrix)
  rownames(targetMatrix) <- targetMatrix$geneName
  targetMatrix <- as.matrix(log(1+targetMatrix[,-1]))
  targetMatrix <- targetMatrix[,cellTable$long_name]
  colSideColors <- (dplyr::select(cellTable, long_name, major_class) %>% arrange(colnames(targetMatrix)))$major_class
  
  #colSideColors <- (dplyr::select(cellTable, long_name, TranscriptomicName) %>% arrange(colnames(targetMatrix)))$TranscriptomicName
  #colSideColors <- as.character(colSideColors == "Sst Cdk6")
  library(RColorBrewer)
  myColors <- brewer.pal(length(unique(colSideColors)),"Set1")
  names(myColors) <- levels(as.factor(colSideColors))
  colSideColors <- myColors[colSideColors]
  
  heatmap(targetMatrix, ColSideColors = colSideColors)
  heatmap.2(targetMatrix, trace = "none", ColSideColors = colSideColors, scale="row")


  ######
  dim(geneGroupResults)
  #how many have identical p-values, AUROC and geneCount - 14%
  
  
  #test a group on the whole aw ranking - no cell type specificity
  goGroupName <- "GO:0007268"
  genesymbols <- goTable[[goGroupName]]
  resultsOfInterest <- aw_result_mouseFilter$gene_symbol
  (auroc <- auroc_analytic(rank(aw_result_mouseFilter$agingPValuesWithDirection), aw_result_mouseFilter$gene_symbol %in% genesymbols ))
  (wilcoxP <- wilcox.test(aw_result_mouseFilter$agingPValuesWithDirection[aw_result_mouseFilter$gene_symbol %in% genesymbols],aw_result_mouseFilter$agingPValuesWithDirection[!(aw_result_mouseFilter$gene_symbol %in% genesymbols)])$p.value)  #, alternative = "less"? force same direction as whole cell type?
}
