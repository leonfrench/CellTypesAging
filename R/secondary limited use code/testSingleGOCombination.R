library(dplyr)
library(hash)

#overlaps with GOGroupByCellTypeTests - interactive setup at bottom
#tasic ranking
load("/Users/lfrench/Desktop/results/CellTypesAging/server results/AW.Tasic.emperical.iter.10000.TranscriptomicName.1479930034/aging_genes_aw_fisher_data.RData", v=T)
load("/Users/lfrench/Desktop/results/CellTypesAging/server results/AW.Tasic.emperical.iter.10000.TranscriptomicName.1479930034/rpkmByType2xHuman.RData", v=T)

#mistry ranking
Mistry.Tasic.emperical.iter.10000.TranscriptomicName.1479930072

#blood ranking


targetCellType <- "Pvalb Tacr3"
targetCellType <- "Vip Gpc3"
targetCellType <- "Vip Sncg"
targetCellType <- "Vip Mybpc1"
targetCellType <- "Sst Cdk6"

goGroupName <-"GO:0007267" #cell-cell signalling
goGroupName <- "GO:0007268" #synaptic transmission

genesymbolsInGOGroup <- goTable[[goGroupName]]

#make supplement table
(resultsForCellType <- subset(rpkmByType2xHuman, CellClassID == targetCellType))
mean(resultsForCellType$log1ExpressionZ) #checking if it's more specific

aw_result_mouseFilterCellSpec <- filter(aw_result_mouseFilter, gene_symbol %in% resultsForCellType$humanGene)
aw_result_mouseFilterCellSpecAndGO <- filter(aw_result_mouseFilter, gene_symbol %in% resultsForCellType$humanGene, gene_symbol %in% genesymbolsInGOGroup)

supplementTable <- arrange(aw_result_mouseFilterCellSpecAndGO, agingPValuesWithDirection) %>% mutate(direction = if_else(agingPValuesWithDirection<0, "Down-regulated", "Up-regulated")) %>% dplyr::select(gene_symbol, p = p_value_aw, q= AW_qvalue, direction )
supplementTable <- inner_join(supplementTable, dplyr::select(ungroup(resultsForCellType), humanGene, SstCdk6EnrichmentZ = log1ExpressionZ), by=c("gene_symbol"= "humanGene"))
###########################
## supplement written
###########################
#write.csv(supplementTable, "/Users/lfrench/Desktop/results/CellTypesAging/Tables/TableS2.csv",row.names=F)


sum(aw_result_mouseFilterCellSpec$gene_symbol %in% genesymbolsInGOGroup)
intersect(genesymbolsInGOGroup, rpkmByType2xHuman$humanGene)

resultsForCellTypeGO <- filter(resultsForCellType, humanGene %in% genesymbolsInGOGroup)
mean(resultsForCellTypeGO$log1ExpressionZ) #more specific than any other gene

(auroc <- auroc_analytic(rank(aw_result_mouseFilterCellSpec$agingPValuesWithDirection), aw_result_mouseFilterCellSpec$gene_symbol %in% genesymbolsInGOGroup ))
(wilcoxP <- wilcox.test(aw_result_mouseFilterCellSpec$agingPValuesWithDirection[aw_result_mouseFilterCellSpec$gene_symbol %in% genesymbolsInGOGroup],aw_result_mouseFilterCellSpec$agingPValuesWithDirection[!(aw_result_mouseFilterCellSpec$gene_symbol %in% genesymbolsInGOGroup)])$p.value)  #, alternative = "less"? force same direction as whole cell type?

#whole list
(auroc <- auroc_analytic(rank(aw_result_mouseFilter$agingPValuesWithDirection), aw_result_mouseFilter$gene_symbol %in% genesymbolsInGOGroup ))
(wilcoxP <- wilcox.test(aw_result_mouseFilter$agingPValuesWithDirection[aw_result_mouseFilter$gene_symbol %in% genesymbolsInGOGroup],aw_result_mouseFilter$agingPValuesWithDirection[!(aw_result_mouseFilter$gene_symbol %in% genesymbolsInGOGroup)])$p.value)  #, alternative = "less"? force same direction as whole cell type?
