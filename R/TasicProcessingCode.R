
############## 

cellTable <- read_csv("./data/Allen data/case study/cell_metadata.csv")
cellClass <- read_csv("./data/Allen data/case study/cell_classification.csv")

#clusterNames1 <- tbl_df(read.xlsx("./data/Allen data/case study/Cluster Names.xlsx",1)) #xlsx reader not available on cluster
clusterNames <- tbl_df(read.csv("./data/Allen data/case study/Cluster Names.csv"))
clusterNames$X <- as.character(clusterNames$X)
clusterNames$X[is.na(clusterNames$X)] <- ""

clusterNames$TranscriptomicName <- paste(clusterNames$Transcriptomic.type,clusterNames$X)
clusterNames$TranscriptomicName <- sub("\\s+$", "", clusterNames$TranscriptomicName) #remove trailing whitespace

#add in the name "TranscriptomicName"
cellClass <- inner_join(cellClass, clusterNames, by= c("primary" = "Final.Cluster.ID"))
colnames(cellClass)[1] <-  "short_name"

if (filterForCoreTranscriptomeTypes) {
  cellsToUseWithLongNames <- filter(cellClass, coretype == "Core")[,1]
} else {
  cellsToUseWithLongNames <- cellClass[,1]
}


cellTable <- inner_join(cellTable, dplyr::select(cellClass, short_name, TranscriptomicName, Transcriptomic.type), by="short_name")

###########################################################################################################################
#choose how we select cell types

#major class
#unique(cellTable$CellClassID <- paste0(cellTable$major_class))
#enrichmentThreshold <- 1.4

#sort(unique(cellTable$CellClassID <- cellTable$sub_class)) 
#unique(cellTable$CellClassID <- paste0(cellTable$major_class,".", cellTable$sub_class)) #~37 types

if (level == "TranscriptomicName") {
  cellTable$CellClassID <- cellTable$TranscriptomicName
} else if (level == "mainClass") {
  cellTable$CellClassID <- cellTable$major_class
  enrichmentThreshold <- 1.4
} else if (level == "Vip.Sst.Pvalb.Excite") {
  cellTable$CellClassID <- NA
  cellTable[!is.na(cellTable$sub_class) & cellTable$sub_class == "Vip", "CellClassID"] <- "Vip"
  cellTable[!is.na(cellTable$sub_class) & cellTable$sub_class == "Pvalb", "CellClassID"] <- "Pvalb"
  cellTable[!is.na(cellTable$sub_class) & grepl("Sst", cellTable$sub_class), "CellClassID"] <- "Sst"
  cellTable[!is.na(cellTable$major_class) & cellTable$major_class == "Excitatory", "CellClassID"] <- "Excitatory"
} else {
  print("unkown level to group cells - set level variable")
  stop()
}

#remove unknown or NA cell type names
#remove NA and unkown?
cellsToUseWithLongNames <- setdiff(cellsToUseWithLongNames, filter(cellTable, is.na(CellClassID) | CellClassID == "Unknown" | major_class == "Unknown") %>% dplyr::select(short_name))
print(paste("Number of cells:",nrow(cellsToUseWithLongNames)))

#extra filtering
#cellsToUseWithLongNames <- setdiff(cellsToUseWithLongNames, filter(cellTable, is.na(CellClassID) | major_class == "Glia") %>% select(short_name))
#cellsToUseWithLongNames <- setdiff(cellsToUseWithLongNames, filter(cellTable, is.na(CellClassID) | major_class == "Endothelial") %>% select(short_name))
#cellsToUseWithLongNames <- setdiff(cellsToUseWithLongNames, filter(cellTable, is.na(CellClassID) | major_class == "Excitatory") %>% select(short_name))
#cellsToUseWithLongNames <- setdiff(cellsToUseWithLongNames, filter(cellTable, is.na(CellClassID) | major_class == "Inhibitory") %>% select(short_name))

#cells with disagreeing labels
#subset(cellTable, short_name %in% cellsToUseWithLongNames$short_name & TranscriptomicName == "SMC Myl9")
#subset(cellTable, short_name %in% cellsToUseWithLongNames$short_name & TranscriptomicName == "L2 Ngb")
#subset(cellTable, short_name %in% cellsToUseWithLongNames$short_name & TranscriptomicName == "L5a Pde1c")
#subset(cellTable, short_name %in% cellsToUseWithLongNames$short_name & TranscriptomicName == "Sncg")
#subset(cellTable, short_name %in% cellsToUseWithLongNames$short_name & TranscriptomicName == "Pvalb Gpx3")
#missclassifiyedAsGlia <- c("L2 Ngb", "L5a Pde1c", "Sncg" ,"Pvalb Gpx3")
##too braod of a removal - fix if used again! -removes whole cell classes
#cellsToUseWithLongNames <- setdiff(cellsToUseWithLongNames, filter(cellTable,TranscriptomicName %in% missclassifiyedAsGlia) %>% select(short_name))

###########################################################################################################################

#add long cell names
coreCellNames <- inner_join(cellsToUseWithLongNames, cellTable[,c("long_name","short_name")], by= "short_name")
#filter cell table
cellTable <- filter(cellTable, long_name %in% coreCellNames$long_name)  #use just the core cell types

rpkmData <- read_csv("./data/Allen data/case study/genes_rpkm.csv")
colnames(rpkmData)[1] <- "geneName"


#find genes with zero counts for all cells
countData <- read_csv("./data/Allen data/case study/genes_counts.csv")
colnames(countData)[1] <- "geneName"
countDataCore <- melt(countData[,c("geneName", coreCellNames$long_name)],factorsAsStrings = TRUE)
countDataCore$variable <- as.character(countDataCore$variable)
countDataCore <- group_by(countDataCore, geneName)
#z-score across the columns (cells)
zeroCountGenes <- summarise(countDataCore, max = max(value, na.rm = TRUE)) %>% filter(max <= 0) %>% dplyr::select(geneName)
zeroCountGenes <- zeroCountGenes$geneName

rm(countData, countDataCore)

#filter cells for those with clear cluster (core cells) and melt
rpkmDataCore <- melt(rpkmData[,c("geneName", coreCellNames$long_name)],factorsAsStrings = TRUE)
rpkmDataCore$variable <- as.character(rpkmDataCore$variable)
rpkmDataCore <- group_by(rpkmDataCore, geneName)
#filter zero count genes
rpkmDataCore <- filter(rpkmDataCore, !(geneName %in% zeroCountGenes))

#convert to log(+1)
rpkmDataCore <- mutate(rpkmDataCore, log1Expression = log(1+value))
#rpkmDataCore <- dplyr::select(rpkmDataCore, -value) #keep original value

colnames(rpkmDataCore)[2] <- "long_name"

#update gene name here


#join with mapping file for new symbol
newSymbols <- read_csv("./data/NewGeneSymbols.csv")
rpkmDataCore <- left_join(rpkmDataCore, newSymbols, by=c("geneName" = "oldSymbol"))
rpkmDataCore <- ungroup(rpkmDataCore)
rpkmDataCore <- mutate(rpkmDataCore, geneName = ifelse(is.na(newSymbol), geneName, newSymbol))
rpkmDataCore <- dplyr::select(rpkmDataCore, -newSymbol)

#add cell type groups
rpkmDataCore <- inner_join(rpkmDataCore, dplyr::select(cellTable,CellClassID,long_name), by="long_name")

#print cell counts
rpkmDataCore <- group_by(rpkmDataCore, CellClassID)
summarise(rpkmDataCore, CellCount = n_distinct(long_name), genes = n_distinct(geneName))
range(summarise(rpkmDataCore, CellCount = n_distinct(long_name), genes = n_distinct(geneName))$CellCount)
