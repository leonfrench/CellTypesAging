linLabMatrix <- read_tsv("./data/Zeisel/expression_mRNA_17-Aug-2014.tsv", col_names=F)

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


#add cell type information

if (level == "TranscriptomicName") {
  cellTable <- dplyr::rename(cellTable, CellClassID= level2class)
} else if (level == "mainClass") {
  cellTable <- dplyr::rename(cellTable, CellClassID= level1class)
} else {
  print("unkown level to group cells - set level variable")
  stop()
}
linLabMelted <- inner_join(linLabMelted, dplyr::select(cellTable,CellClassID,cell_id), by="cell_id")

linLabMelted <- dplyr::rename(linLabMelted, long_name=cell_id)
linLabMelted <- dplyr::select(linLabMelted, -moleculeCount)
rpkmDataCore <- linLabMelted #rpkmDataCore is done
cellTable <- dplyr::rename(cellTable, long_name=cell_id)
rm(linLabMatrix)