library(readr)
library(dplyr)
#combine results
resultsBA1147 <- read.csv("/Users/lfrench/Desktop/results/CellTypesAging/server results/AW.Tasic.emperical.iter.10000.TranscriptomicName.1479930034/CellEnrichmentEmpericalPMinCols.csv", row.names=1)
cor.test(resultsBA1147$empericalP, resultsBA1147$CellCount, m='s')
plot(resultsBA1147$empericalP, resultsBA1147$CellCount, m='s')
plot(resultsBA1147$AUROC, resultsBA1147$CellCount, m='s')
cor.test(resultsBA1147$AUROC, resultsBA1147$CellCount, m='s')


plot(resultsBA1147$AUROC, resultsBA1147$geneCount, m='s')
cor.test(resultsBA1147$AUROC, resultsBA1147$geneCount, m='s')

plot(resultsBA1147$AUROC, resultsBA1147$CellCount, m='s')
cor.test(resultsBA1147$AUROC, resultsBA1147$CellCount, m='s')


plot(resultsBA1147$empericalP, resultsBA1147$CellCount, m='s')
cor.test(resultsBA1147$empericalP, resultsBA1147$CellCount, m='s')

plot(resultsBA1147$empericalP, resultsBA1147$geneCount, m='s')
cor.test(resultsBA1147$empericalP, resultsBA1147$geneCount, m='s')


resultsBA1147 <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/server results/AW.Tasic.emperical.iter.10000.TranscriptomicName.1479930034/CellEnrichmentEmpericalPMinCols.csv")
resultsMistry <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/server results/Mistry.Tasic.emperical.iter.10000.TranscriptomicName.1479930072/CellEnrichmentEmpericalPMinCols.csv")
resultsBlood <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/server results/Blood.Tasic.emperical.iter.10000.TranscriptomicName.1479930030/CellEnrichmentEmpericalPMinCols.csv")

resultsMerged <- inner_join(resultsBA1147, resultsMistry, by="CellClassID",suffix=c(".BA1147", ".Mistry"))
resultsBlood <- resultsBlood %>% setNames(paste0(names(.), ".Blood"))
resultsMerged <- inner_join(resultsMerged, resultsBlood, by=c("CellClassID" = "CellClassID.Blood"))

write.csv(resultsMerged, "/Users/lfrench/Desktop/results/CellTypesAging/server results/Tasic.TranscriptomicName.combined.csv")
resultsMergedSlim <- select(resultsMerged, CellClassID, AUROC.BA1147, p.BA1147 = empericalP.BA1147, q.BA1147 = empericalFDRP.BA1147, AUROC.Mistry, p.Mistry= empericalP.Mistry, AUROC.Blood, p.Blood = empericalP.Blood)

resultsMergedSlim <- resultsMergedSlim %>% mutate_each(funs(signif(.,2)), -CellClassID) 

write.csv(resultsMergedSlim, "/Users/lfrench/Desktop/results/CellTypesAging/server results/Tasic.TranscriptomicName.combined.slim.csv",row.names=F)


#combine Zeisel results
zresultsBA1147 <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/server results/AW.Zeisel.emperical.iter.10000.TranscriptomicName.1480347201/CellEnrichmentEmpericalPMinCols.csv")
zresultsMistry <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/server results/Mistry.Zeisel.emperical.iter.10000.TranscriptomicName.1480347194/CellEnrichmentEmpericalPMinCols.csv")
zresultsBlood <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/server results/Blood.Zeisel.emperical.iter.10000.TranscriptomicName.1480347179/CellEnrichmentEmpericalPMinCols.csv")
resultsBA1147 <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/server results/AW.Tasic.emperical.iter.10000.TranscriptomicName.1479930034/CellEnrichmentEmpericalPMinCols.csv")
resultsBA1147$directionTasicBA1147 <- sign(resultsBA1147$AUROC - resultsBA1147$meanEmpAUC)

zresultsMerged <- inner_join(zresultsBA1147, zresultsMistry, by="CellClassID",suffix=c(".BA1147", ".Mistry"))
zresultsBlood <- zresultsBlood %>% setNames(paste0(names(.), ".Blood"))
zresultsMerged <- inner_join(zresultsMerged, zresultsBlood, by=c("CellClassID" = "CellClassID.Blood"))

zresultsMerged

#load mapping file
ztoTasicMap <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/data/Top6TypesMappedToZeisel.csv")

zresultsMergedTop <- inner_join(zresultsMerged, ztoTasicMap, by = c("CellClassID" ="Zeisel.CellClassID"))
#add in Tasic results
zresultsMergedTop <- inner_join(zresultsMergedTop,select(resultsBA1147, CellClassID, empericalFDRP,directionTasicBA1147, Tasic.BA1147.AUROC = AUROC), by=c("Tasic.CellClassID" = "CellClassID"))

zresultsMergedTop$direction.BA1147 <- sign(zresultsMergedTop$AUROC.BA1147 - zresultsMergedTop$meanEmpAUC.BA1147 )
zresultsMergedTop$direction.Mistry <- sign(zresultsMergedTop$AUROC.Mistry-zresultsMergedTop$meanEmpAUC.Mistry)
#plus direction in Tasic results
#test if direction agrees
zresultsMergedTop <- mutate(zresultsMergedTop, directionAgrees=abs(directionTasicBA1147 + direction.BA1147 + direction.Mistry) == 3)
zresultsMergedTop$directionAgrees

print.data.frame(zresultsMergedTop)

write.csv(zresultsMergedTop, "/Users/lfrench/Desktop/results/CellTypesAging/server results/Zeisel.top.TranscriptomicName.combined.csv",row.names=F)

zresultsMergedTopSlim <- zresultsMergedTop %>% mutate_each(funs(signif(.,2)), -CellClassID, -Tasic.CellClassID, -Tasic.label) 
zresultsMergedTopSlim <- select(zresultsMergedTopSlim, Zeisel.label = CellClassID, Cell.Count = CellCount.BA1147, Tasic.label, Tasic.BA1147.AUROC = Tasic.BA1147.AUROC, Tasic.BA1147.q = empericalFDRP, AUROC.BA1147, p.BA1147 = empericalP.BA1147, AUROC.Mistry,p.Mistry = empericalP.Mistry, empericalP.Blood) %>% arrange(Tasic.BA1147.q,Zeisel.label)
colnames(zresultsMergedTopSlim) <- gsub("Mistry", "Prefrontal", colnames(zresultsMergedTopSlim))

write.csv(zresultsMergedTopSlim, "/Users/lfrench/Desktop/results/CellTypesAging/server results/Zeisel.top.TranscriptomicName.combined.final.csv",row.names=F)
