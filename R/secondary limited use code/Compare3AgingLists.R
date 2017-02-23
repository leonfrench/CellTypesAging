library(scatterD3)
#temp code
joinedResults <- inner_join(geneCellResultsOrig, geneCellResultsBlood, by="CellClassID", suffix=c(".brain",".blood"))
joinedResults <- inner_join(geneCellResults, joinedResults, by="CellClassID") #for mistry
joinedResults <- dplyr::rename(joinedResults, AUROC.mistry = AUROC)

scatterD3(x = joinedResults$AUROC.mistry, y = joinedResults$AUROC.brain, lab = joinedResults$CellClassID)
scatterD3(x = joinedResults$AUROC.blood, y = joinedResults$AUROC.brain, lab = joinedResults$CellClassID)
cor.test(joinedResults$AUROC.mistry, joinedResults$AUROC.brain)
cor.test(joinedResults$AUROC.mistry, joinedResults$AUROC.blood)
cor.test(joinedResults$AUROC.brain, joinedResults$AUROC.blood)
joinedResults$AUCROC.meanBrain <- (joinedResults$AUROC.mistry +  joinedResults$AUROC.brain)/2

scatterD3(x = joinedResults$AUROC.blood, y = joinedResults$AUCROC.meanBrain, lab = joinedResults$CellClassID)

lmResult <- lm(AUCROC.meanBrain ~ AUROC.blood, data= joinedResults)
joinedResults$AUROC.bloodResiduals <- lmResult$residuals
joinedResults[order(joinedResults$AUROC.bloodResiduals),c("CellClassID", "AUCROC.meanBrain", "geneCount", "AUROC.bloodResiduals")]

cor.test(joinedResults$AUROC.bloodResiduals,  joinedResults$AUROC.blood,m='s')

scatterD3(x = joinedResults$AUROC.bloodResiduals, y = joinedResults$AUCROC.meanBrain, lab = joinedResults$CellClassID)


bloodAging <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/Blood aging list/GOGroupResults.vMajor.reducedCols.csv")
bloodAging <- mutate(bloodAging, agingPValuesWithDirection = ifelse(Direction == "+", 1-P, -1+P))
bloodAging <- dplyr::rename(bloodAging, gene_symbol = `NEW-Gene-ID`) %>% dplyr::rename(p_value_aw = P)


##############
mistryAgeDown <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/Meeta Age Tables/Supplementary Table 10.ageDown.csv")
mistryAgeDown$direction <- -1
mistryAgeUp <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/Meeta Age Tables/Supplementary Table 10.ageUp.csv")
mistryAgeUp$direction <- 1
mistryAge <- bind_rows(mistryAgeUp, mistryAgeDown)

#take direction with best pvalue
mistryAge <- mistryAge %>% group_by(GeneSymbol) %>% dplyr::slice(which.min(Pvalue)) %>% arrange(Pvalue)
mistryAge <- mutate(mistryAge, agingPValuesWithDirection = ifelse(direction > 0, 1-Pvalue, -1+Pvalue))
mistryAge <- dplyr::rename(mistryAge, gene_symbol = GeneSymbol) %>% dplyr::rename(p_value_aw = Pvalue)

joined <- inner_join(dplyr::select(bloodAging, gene_symbol, agingPValuesWithDirection), dplyr::select(mistryAge, gene_symbol, agingPValuesWithDirection), by="gene_symbol",suffix=c(".Blood",".Mistry"))

#load("/Users/lfrench/Desktop/results/CellTypesAging/aging_RLM_results.RData",verbose = T) #for direction
#load("/Users/lfrench/Desktop/results/CellTypesAging/aging_corrected_gene_symbol.RData",verbose=T)
#load("/Users/lfrench/Desktop/results/CellTypesAging/aging_genes_aw_fisher_result.RData",verbose=T)
#aw_result$p_value_aw <- as.numeric(aw_result$p_value_aw)
#Gene_symbol <- Corrected_Gene_symbol
#aw_result$gene_symbol <- Gene_symbol
#subset(aw_result, gene_symbol=="SEPT14")
#aw_result$RLM.coef.1 <- summary_rlm_BA11$RLM.coef
#aw_result$RLM.coef.2 <- summary_rlm_BA47$RLM.coef
#aw_result$RLM.coef.meta <- (aw_result$study_weight.1 * aw_result$RLM.coef.1 +  aw_result$study_weight.2 * aw_result$RLM.coef.2)/aw_result$total_weight
#length(unique(Gene_symbol))

load("/Users/lfrench/Desktop/results/CellTypesAging/AgingFullResults_UniqueGenes.RData",verbose=T)
all(rownames(p.dat.unique.genes) == rownames(es.dat.unique.genes))
all(rownames(p.dat.unique.genes) == rownames(aw.out))
aw.out <- data.frame(aw.out)
aw.out$RLM.coef.BA11 <- es.dat.unique.genes[,"BA11"]
aw.out$RLM.coef.BA47 <- es.dat.unique.genes[,"BA47"]
aw.out$gene_symbol <- rownames(aw.out)
aw.out$RLM.coef.meta <- (aw.out$weight_BA11 * aw.out$RLM.coef.BA11 +  aw.out$weight_BA47 * aw.out$RLM.coef.BA47)/(aw.out$weight_BA11 + aw.out$weight_BA47)
aw_result <- tbl_df(aw.out)


#aw_result <- filter(aw_result, gene_symbol != "---")
#show duplicated genes
#group_by(aw.out, gene_symbol) %>% summarise(n= n()) %>% arrange(n) %>% tail(n=20)

#pick the best p-value per duplicated gene
#aw_result <- aw_result %>% group_by(gene_symbol) %>% arrange(p_value_aw) %>% filter(row_number()==1)
aw_result <- dplyr::rename(aw_result, p_value_aw = AW_pvalue)

#convert p-value to -1 to 1 scale 
aw_result[aw_result$RLM.coef.meta >= 0,"agingPValuesWithDirection"] <- 1-aw_result[aw_result$RLM.coef.meta >= 0, "p_value_aw"]
aw_result[aw_result$RLM.coef.meta < 0,"agingPValuesWithDirection"] <- -1+aw_result[aw_result$RLM.coef.meta <0, "p_value_aw"]


joined <- inner_join(joined, dplyr::select(aw_result, gene_symbol, agingPValuesWithDirection), by="gene_symbol")
