library(readr)
library(dplyr)

#duplicated code
bloodAging <- read_csv("./data/Blood aging list/GOGroupResults.vMajor.reducedCols.csv")
bloodAging <- dplyr::rename(bloodAging, gene_symbol = `NEW-Gene-ID`) %>% dplyr::rename(p_value_aw = P)
bloodAging$agingPValuesWithDirection <- rank(-1*bloodAging$p_value_aw)
bloodAging[bloodAging$Direction == "-","agingPValuesWithDirection"] <- -1*bloodAging[bloodAging$Direction == "-", "agingPValuesWithDirection"]

mistryAgeDown <- read_csv("./data/Mistry Age Tables/Supplementary Table 10.ageDown.csv")
mistryAgeDown$direction <- -1
mistryAgeUp <- read_csv("./data/Mistry Age Tables/Supplementary Table 10.ageUp.csv")
mistryAgeUp$direction <- 1
mistryAge <- bind_rows(mistryAgeUp, mistryAgeDown)

#take direction with best pvalue
mistryAge <- mistryAge %>% group_by(GeneSymbol) %>% dplyr::slice(which.min(Pvalue)) %>% arrange(Pvalue)
mistryAge <- mutate(mistryAge, agingPValuesWithDirection = ifelse(direction > 0, 1-Pvalue, -1+Pvalue))
mistryAge <- dplyr::rename(mistryAge, gene_symbol = GeneSymbol) %>% dplyr::rename(p_value_aw = Pvalue)


load("./data/AWLists/AgingFullResults_UniqueGenes.RData",verbose=T)
all(rownames(p.dat.unique.genes) == rownames(es.dat.unique.genes))
all(rownames(p.dat.unique.genes) == rownames(aw.out))
aw.out <- data.frame(aw.out)
aw.out$RLM.coef.BA11 <- es.dat.unique.genes[,"BA11"]
aw.out$RLM.coef.BA47 <- es.dat.unique.genes[,"BA47"]
aw.out$gene_symbol <- rownames(aw.out)

#aw.out <- subset(aw.out, RLM.coef.BA11 * RLM.coef.BA47 > 0) #filter out different directions

aw.out$RLM.coef.meta <- (aw.out$weight_BA11 * aw.out$RLM.coef.BA11 +  aw.out$weight_BA47 * aw.out$RLM.coef.BA47)/(aw.out$weight_BA11 + aw.out$weight_BA47)
aw_result <- tbl_df(aw.out)

#pick the best p-value per duplicated gene
aw_result <- dplyr::rename(aw_result, p_value_aw = AW_pvalue)

#convert p-value to -1 to 1 scale - not used - changed due to ties at 1/-1
#aw_result[aw_result$RLM.coef.meta >= 0,"agingPValuesWithDirection"] <- 1-aw_result[aw_result$RLM.coef.meta >= 0, "p_value_aw"]
#aw_result[aw_result$RLM.coef.meta < 0,"agingPValuesWithDirection"] <- -1+aw_result[aw_result$RLM.coef.meta <0, "p_value_aw"]

aw_result$agingPValuesWithDirection <- rank(-1*aw_result$p_value_aw)
aw_result[aw_result$RLM.coef.meta < 0,"agingPValuesWithDirection"] <- -1*aw_result[aw_result$RLM.coef.meta <0, "agingPValuesWithDirection"]



joined <- inner_join(bloodAging, aw_result, by="gene_symbol", suffix=c(".blood", ".BA1147"))
joined <- inner_join(joined, mistryAge, by="gene_symbol", suffix=c("", "Mistry"))
joined <- rename(joined, agingPValuesWithDirection.Mistry = agingPValuesWithDirection)
cor.test(joined$agingPValuesWithDirection.Mistry, joined$agingPValuesWithDirection.BA1147,m='s')
plot(joined$agingPValuesWithDirection.Mistry, joined$agingPValuesWithDirection.BA1147)
cor.test(joined$agingPValuesWithDirection.blood, joined$agingPValuesWithDirection.BA1147,m='s')
cor.test(joined$agingPValuesWithDirection.blood, joined$agingPValuesWithDirection.Mistry,m='s')

#check duplicates
length(unique(joined$agingPValuesWithDirection.Mistry))
length(joined$agingPValuesWithDirection.blood)
length(unique(joined$agingPValuesWithDirection.blood))


