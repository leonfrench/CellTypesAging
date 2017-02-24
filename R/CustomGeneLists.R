library(tidyr)
library(readr)
library(homologene)
source("./R/AUCFunction.R") #load AUC function

#Other cell type lists
prefixes <- c("Zeisel", "Cahoy", "Darmanis", "Custom")
prefixes <- c("Darmanis")
#prefixes <- c("Custom")
#prefixes <- c("MouseSource")
resultposition <- 0
geneCellResultsOtherLists <- data.frame( CellClassID=character(), geneCount= numeric(), p=numeric(), stringsAsFactors=F)



for(agingGeneSource in c("Blood", "Mistry", "AW")) {
  #warning - duplicated code
  if (agingGeneSource == "Blood") {
    bloodAging <- read_csv("./data/Blood aging list/GOGroupResults.vMajor.reducedCols.csv")
    #bloodAging <- mutate(bloodAging, agingPValuesWithDirection = ifelse(Direction == "+", 1-P, -1+P)) #not used due to precision trouble
    bloodAging <- dplyr::rename(bloodAging, gene_symbol = `NEW-Gene-ID`) %>% dplyr::rename(p_value_aw = P)
    
    #use ranks for ordering due to precision trouble
    bloodAging$agingPValuesWithDirection <- rank(-1*bloodAging$p_value_aw)
    bloodAging[bloodAging$Direction == "-","agingPValuesWithDirection"] <- -1*bloodAging[bloodAging$Direction == "-", "agingPValuesWithDirection"]
    
    aw_result <- bloodAging
  } else if (agingGeneSource == "Mistry") {
    mistryAgeDown <- read_csv("./data/Mistry Age Tables/Supplementary Table 10.ageDown.csv")
    mistryAgeDown$direction <- -1
    mistryAgeUp <- read_csv("./data/Mistry Age Tables/Supplementary Table 10.ageUp.csv")
    mistryAgeUp$direction <- 1
    mistryAge <- bind_rows(mistryAgeUp, mistryAgeDown)
    
    #take direction with best pvalue
    mistryAge <- mistryAge %>% group_by(GeneSymbol) %>% dplyr::slice(which.min(Pvalue)) %>% arrange(Pvalue)
    mistryAge <- mutate(mistryAge, agingPValuesWithDirection = ifelse(direction > 0, 1-Pvalue, -1+Pvalue))
    mistryAge <- dplyr::rename(mistryAge, gene_symbol = GeneSymbol) %>% dplyr::rename(p_value_aw = Pvalue)
    aw_result <- mistryAge
  } else if (agingGeneSource == "AW") {
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
    
  }
  
  for(prefix in prefixes) {
    #IF using Darminasis human genes, switch to aw_result! as it's not filtered for mouse reachable genes!
    if (prefix == "Darmanis" | prefix == "Custom") {
      print("using unfiltered aw result")
      agingData <- aw_result
    } else {
      reachableHumanGenes <- mouse2human(human2mouse(aw_result$gene_symbol)$mouseGene)$humanGene
      agingData <- filter(aw_result, gene_symbol %in% reachableHumanGenes)
    }
    
    for(filename in list.files("./data/GeneLists/", pattern = paste0(prefix,".*txt"), full.names = T)) {
      print(filename)
      #skip fetal cell types
      if (grepl("Darmanis.QuiescentNewlyBornNeurons.txt", filename)  || grepl( "ReplicatingNeuronalProgenitors",filename)) next()
      resultposition <- resultposition + 1
      genesOfInterest <- read.csv(filename,header=F,stringsAsFactors = F)[,"V1"]
      #IF USING human genes, switch to aw_result!!
      labelMask <- agingData$gene_symbol %in% genesOfInterest
      auroc <- auroc_analytic(rank(agingData$agingPValuesWithDirection), labelMask)
      wilcoxP <- wilcox.test(agingData[!labelMask,]$agingPValuesWithDirection, agingData[labelMask, ]$agingPValuesWithDirection)$p.value
      
      geneCellResultsOtherLists[resultposition, "ranking"] <- agingGeneSource
      
      geneCellResultsOtherLists[resultposition, "geneCount"] <- sum(labelMask)
      geneCellResultsOtherLists[resultposition, "CellClassID"] <- gsub(paste0(".*/",prefix,"."),"", filename)
      geneCellResultsOtherLists[resultposition, "CellClassID"] <- gsub(paste0(".txt"),"", geneCellResultsOtherLists[resultposition, "CellClassID"])
      geneCellResultsOtherLists[resultposition, "AUROC"] <- auroc
      geneCellResultsOtherLists[resultposition, "p"] <- wilcoxP
      
    }
    
    write.csv(geneCellResultsOtherLists, file=paste0(outputFolder,prefix,".CellTypeLists.csv"),row.names=F)
  }
}

#correct for the number of cell types but not the rankings
geneCellResultsOtherLists$p.adjusted <- geneCellResultsOtherLists$p * length(unique(geneCellResultsOtherLists$CellClassID))
(geneCellResultsOtherLists <- geneCellResultsOtherLists[order(geneCellResultsOtherLists$p),])
geneCellResultsOtherLists


addStar <- function(x) {
  #print(x)
  pAdjust <- as.numeric(x["p.adjusted"])
  AUROC <- x["AUROC"]
  print(pAdjust)
  print(AUROC)
  result <- signif(as.numeric(AUROC),2)
  if (pAdjust < 0.0005) {
    result <- paste0(result, " ***")
  } else if (pAdjust < 0.005) {
    result <- paste0(result, " **")
  } else if(pAdjust < 0.05) {
    result <- paste0(result, " *")
  }
  print(result)
  as.character(result)
}
geneCellResultsOtherLists$AUROCStar <- apply(geneCellResultsOtherLists, 1, addStar)
AUCTable <- spread(dplyr::select(geneCellResultsOtherLists, -AUROC, -p, -geneCount, -p.adjusted), ranking, AUROCStar)
(AUCTable <- dplyr::select(tbl_df(AUCTable), Name = CellClassID, BA1147=AW, Mistry, Blood))

write.csv(AUCTable, file=paste0(outputFolder,prefix,".OtherLists.ROC.csv"),row.names=F)

#spread(dplyr::select(geneCellResultsOtherLists, -AUROC, -p, -geneCount), ranking, p.adjusted)
#spread(dplyr::select(geneCellResultsOtherLists, -p.adjusted, -p, -geneCount), ranking, AUROC)
#spread(dplyr::select(geneCellResultsOtherLists, -p.adjusted, -p, -AUROC), ranking, geneCount)


write.csv(file = "/Users/lfrench/Desktop/results/CellTypesAging/data/GeneLists/Custom.housekeeping.intersect.SynTranGO.txt", intersect(goTable[[goGroupName]], read.csv("/Users/lfrench/Desktop/results/CellTypesAging/data/GeneLists/Custom.housekeeping.txt",header=F,stringsAsFactors = F)[,"V1"]),row.names=F, quote=F)
