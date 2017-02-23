pastList <- read.table("/Users/lfrench/Desktop/results/CellTypesAging/Summary_SingleCellEnrichment.order.txt",sep="|")
pastList$oldOrder <- rownames(pastList)
newList <- read.csv("/Users/lfrench/Desktop/results/CellTypesAging/results/AW.emperical.TranscriptomicName.49.1468342974/CellEnrichmentEmpericalP.csv")
joined <- dplyr::inner_join(newList, pastList,by=c("CellClassID"= "V1"))

cor.test(as.numeric(joined$oldOrder), as.numeric(joined$empericalP),m='s')
