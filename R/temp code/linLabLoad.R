library(homologene)
library(dplyr)
library(reshape2)
library(readr)

linLabMatrix <- read_tsv("/Users/lfrench/Downloads/expression_mRNA_17-Aug-2014.tsv", col_names=F)

cellTable <- tbl_df(t(linLabMatrix[1:10,2:ncol(linLabMatrix)]))
colnames(cellTable) <- cellTable[1,]
(cellTable <- cellTable[-1,])

colnames(linLabMatrix) <- linLabMatrix[8,]
linLabMatrix[1:15,1:15]
linLabMatrix <- linLabMatrix[12:nrow(linLabMatrix),]
linLabMatrix <- dplyr::select(linLabMatrix, -cell_id)
colnames(linLabMatrix)[1] <- "geneName"

#now it can be melted
linLabMelted <- reshape2::melt(linLabMatrix,factorsAsStrings = TRUE, id.vars=c("geneName"), variable.name = "cell_id", value.name = "moleculeCount")
linLabMelted <- tbl_df(linLabMelted)
linLabMelted$moleculeCount <- as.numeric(linLabMelted$moleculeCount)

linLabMelted <- mutate(linLabMelted, log1Expression = log(1+moleculeCount))

#summarized across cell types
linLabMelted <- inner_join(linLabMelted, dplyr::select(cellTable,level2class,cell_id), by="cell_id")

linLabMelted <- dplyr::rename(linLabMelted, CellClassID= level2class)

linLabMelted <- group_by(linLabMelted, CellClassID, geneName)
#average across cell types
rpkmByType <- summarise(linLabMelted, log1Expression = mean(log1Expression))

#match other source by creating z-score
rpkmByType <- group_by(rpkmByType, geneName)

#scale/z-score - do z score after averaging in a cell type
rpkmByType <- mutate(rpkmByType, log1ExpressionZ = (log1Expression - mean(log1Expression)) / sd(log1Expression))

backgroundHumanGenes <- tbl_df(mouse2human(rpkmByType$geneName))

#needs work here
#rpkmByType
#geneCellResults <- left_join(geneCellResults, rpkmDataCore %>% group_by(CellClassID) %>% summarise(CellCount = n_distinct(long_name)), by="CellClassID")

#cellCountsTable <- 

#create empirical p-values
#############################################################################################
start.time <- Sys.time()
cellTable$CellClassID <- cellTable$level2class
rpkmDataCore <- linLabMelted
rownames(geneCellResults) <- geneCellResults$CellClassID
for (iteration in 1:iterations) {
  
  print(paste(iteration, ", time: " , Sys.time() - start.time))
  cellTableShuffled <- cellTable
  cellTableShuffled$CellClassID <- sample(cellTableShuffled$CellClassID)

  rpkmDataCore <- rpkmDataCore %>% ungroup() %>% select(-CellClassID)
  rpkmDataCore <- inner_join(rpkmDataCore, dplyr::select(cellTableShuffled,CellClassID,cell_id), by="cell_id")
  
  rpkmDataCore <- group_by(rpkmDataCore, CellClassID, geneName) #slow
  #average across cell types
  rpkmByType <- summarise(rpkmDataCore, log1Expression = mean(log1Expression)) %>% group_by(geneName)
  
  #scale/z-score - do z score after averaging in a cell type
  rpkmByType <- mutate(rpkmByType, log1ExpressionZ = (log1Expression - mean(log1Expression)) / sd(log1Expression))
  
  #cut by enrichment - leaving how many genes per cell type?
  rpkmByType2x <- filter(rpkmByType, log1ExpressionZ > enrichmentThreshold) #enrichmentThreshold defined at top
  
  rpkmByType <- group_by(rpkmByType, geneName)
  stableGenes <- summarise(rpkmByType, log1ExpressionZ = max(log1ExpressionZ)) %>% filter(log1ExpressionZ < enrichmentThreshold)
  
  #add stable genes to rpkmByType2x
  stableGenes$CellClassID <- "CelltypeNonSpecific"
  
  rpkmByType2x <- bind_rows(dplyr::select(rpkmByType2x,CellClassID,log1ExpressionZ, geneName), stableGenes)
  
  #convert to human genes - remove those with no ortholog
  rpkmByType2xHuman <- inner_join(rpkmByType2x, backgroundHumanGenes, by=c("geneName" = "mouseGene"))
  
  #filter out aw_result human genes that cannot be reached from the mouse genes
  aw_result_mouseFilter  <- filter(aw_result, gene_symbol %in% backgroundHumanGenes$humanGene)
  
  #remove genes not in the aging dataset
  rpkmByType2xHuman <- filter(rpkmByType2xHuman, humanGene %in% aw_result_mouseFilter$gene_symbol)
  
  #now the gene lists are ready
  rpkmByType2xHuman <- dplyr::inner_join(rpkmByType2xHuman, dplyr::select(aw_result_mouseFilter, gene_symbol, agingPValuesWithDirection), by=c("humanGene"="gene_symbol"))
  
  #iterate cell type lists here
  
  for(targetCellType in unique(rpkmByType2xHuman$CellClassID)) {
    resultsOfInterest <- subset(rpkmByType2xHuman, CellClassID == targetCellType)
    genesOfInterest <- resultsOfInterest$humanGene
    labelMask <- aw_result_mouseFilter$gene_symbol %in% genesOfInterest
    auroc <- auroc_analytic(rank(aw_result_mouseFilter$agingPValuesWithDirection), labelMask)
    #wilcoxP <- wilcox.test(aw_result_mouseFilter[!labelMask,]$agingPValuesWithDirection, aw_result_mouseFilter[labelMask, ]$agingPValuesWithDirection)$p.value
    
    #geneCellResults[targetCellType, paste0("geneCount",iteration)] <- length(genesOfInterest)
    #geneCellResults[targetCellType, paste0("CellClassID",iteration)] <- targetCellType
    geneCellResults[targetCellType, paste0("AUROC",iteration)] <- auroc
    #geneCellResults[targetCellType, paste0("p",iteration)] <- wilcoxP
  }
}
Sys.time() - start.time