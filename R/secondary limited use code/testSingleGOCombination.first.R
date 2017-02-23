#get genes for a specific group and cell type
targetCellType <- "Sst Cdk6"
goGroupName <- "GO:0007268"
#goGroupName <-"GO:0007267" #cell-cell signalling
resultsOfInterest <- subset(rpkmByType2xHuman, CellClassID == targetCellType)



#filter for genes with GO group annotations
genesymbols <- goTable[[goGroupName]]
resultsOfInterest <- dplyr::filter(resultsOfInterest, humanGene %in% genesymbols)

cellCellSignal <- resultsOfInterest$geneName
synapticTrans <- resultsOfInterest$geneName

geneSubsetTable <- resultsOfInterest %>% ungroup() %>% dplyr::select(humanGene, log1ExpressionZ, agingPValuesWithDirection) %>% arrange(agingPValuesWithDirection)
write.table(geneSubsetTable, file=paste0(outputFolder, "specGenes.", targetCellType,".", goGroupName,".tsv"),sep="\t")
#########################
#write heatmap
mouseGenes <- c("Cdk6", human2mouse(geneSubsetTable$humanGene)$mouseGene)
targetMatrix <- filter(rpkmDataCore, geneName %in% mouseGenes)
targetMatrix <- group_by(targetMatrix, geneName, CellClassID)
targetMatrix <- summarise(targetMatrix, expression = mean(log1Expression))
targetMatrix <- dcast(targetMatrix, geneName ~ CellClassID)

targetMatrix <- as.data.frame(targetMatrix)
rownames(targetMatrix) <- targetMatrix$geneName
targetMatrix <- as.matrix(targetMatrix[,-1])
heatmap(targetMatrix)

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


print.data.frame(filter(rpkmDataCore, geneName == "Adcy7" & log1Expression > 0))
hist(filter(rpkmDataCore, geneName == "Adcy7" & log1Expression > 0)$log1Expression)

cdk6 <- filter(rpkmDataCore, geneName == "Cdk6") %>% ungroup() 
cdk6 <- group_by(cdk6, CellClassID) 
cdk6 <- dplyr::summarise(cdk6, log1Expression = mean(log1Expression)) 
cdk6 <- mutate(cdk6, exp = (log1Expression - mean(log1Expression))/ sd(log1Expression))
print(cdk6,n=50)
filter(cdk6, exp > 2)
