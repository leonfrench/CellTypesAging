
#Other cell type lists
prefixes <- c("Zeisel", "Cahoy", "Darmanis", "Custom")
for(prefix in prefixes) {
  #IF using Darminasis human genes, switch to aw_result! as it's not filtered for mouse reachable genes!
  if (prefix == "Darmanis" | prefix == "Custom") {
    print("using unfiltered aw result")
    agingData <- aw_result
  } else {
    agingData <- aw_result_mouseFilter
  }
  geneCellResultsOtherLists <- data.frame( CellClassID=character(), geneCount= numeric(), p=numeric(), stringsAsFactors=F)
  resultposition <- 0
  for(filename in list.files("/Users/lfrench/Desktop/data/Gene Lists/Human/", pattern = paste0(prefix,".*txt"), full.names = T)) {
    print(filename)
    #skip fetal cell types
    if (grepl("Darmanis.QuiescentNewlyBornNeurons.txt", filename)  || grepl( "ReplicatingNeuronalProgenitors",filename)) next()
    resultposition <- resultposition +1
    genesOfInterest <- read.csv(filename,header=F,stringsAsFactors = F)[,"V1"]
    #IF USING human genes, switch to aw_result!!
    labelMask <- agingData$gene_symbol %in% genesOfInterest
    auroc <- auroc_analytic(rank(agingData$agingPValuesWithDirection), labelMask)
    wilcoxP <- wilcox.test(agingData[!labelMask,]$agingPValuesWithDirection, agingData[labelMask, ]$agingPValuesWithDirection)$p.value
    
    geneCellResultsOtherLists[resultposition, "geneCount"] <- length(genesOfInterest)
    geneCellResultsOtherLists[resultposition, "CellClassID"] <- gsub(paste0(".*/",prefix,"."),"", filename)
    geneCellResultsOtherLists[resultposition, "CellClassID"] <- gsub(paste0(".txt"),"", geneCellResultsOtherLists[resultposition, "CellClassID"])
    geneCellResultsOtherLists[resultposition, "AUROC"] <- auroc
    geneCellResultsOtherLists[resultposition, "p"] <- wilcoxP
    
  }
  
  geneCellResultsOtherLists$p.adjusted <- p.adjust(geneCellResultsOtherLists$p)
  (geneCellResultsOtherLists <- geneCellResultsOtherLists[order(geneCellResultsOtherLists$p),])
  
  write.csv(geneCellResultsOtherLists, file=paste0(outputFolder,prefix,".CellTypeLists.csv"),row.names=F)
}

