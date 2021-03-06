library(gplots)
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(annotate)
library(GO.db)
library(tmod)
library(readr)
library(homologene)
library(ggrepel)

targetVariable = "regression"
#targetVariable = "female"
#targetVariable = "male"

metaRegressionResult <- read_csv("./data/metaA/MetaRegression_SexSpecific_2-17-17_MLS.csv")
colnames(metaRegressionResult)[1] <- "SYMBOL"

#fix gene symbols
goodGeneNames <- read_csv("./data/metaA/MDD-metaAR_8cohorts_Final.csv")
goodGeneNames$namesUpper <- toupper(goodGeneNames$SYMBOL)
metaRegressionResult$namesUpper <- toupper(metaRegressionResult$SYMBOL)
metaRegressionResult <- left_join(metaRegressionResult, dplyr::select(goodGeneNames, namesUpper,SYMBOL), by=c("namesUpper"="namesUpper")) %>% dplyr::select(-namesUpper)
metaRegressionResult %<>% mutate(SYMBOL = if_else(!is.na(SYMBOL.y),SYMBOL.y,SYMBOL.x)) %>% dplyr::select(SYMBOL, everything()) %>% dplyr::select(-SYMBOL.x, -SYMBOL.y)


if (targetVariable == "male") {
  metaRegressionResult$medianEffect <- metaRegressionResult$EffectSize_Males #median effect across the studies (split by sexes)
  metaRegressionResult$p_value_aw <- metaRegressionResult$pvalue_Males
} else if (targetVariable == "female") {
  metaRegressionResult$medianEffect <- metaRegressionResult$EffectSize_Females #median effect across the studies (split by sexes)
  metaRegressionResult$p_value_aw <- metaRegressionResult$pvalue_Females
} else if (targetVariable == "regression") {
  metaRegressionResult$medianEffect <- metaRegressionResult$EffectSize_MetaR_sex #median effect across the studies (split by sexes)
  metaRegressionResult$p_value_aw <- metaRegressionResult$pvalue_MetaR_sex
} else {
  metaRegressionResult <- NULL
  stop("no target sex given")
}

metaRegressionResult <- dplyr::select(metaRegressionResult, gene_symbol = SYMBOL, p_value_aw, medianEffect, EffectSize_Females, EffectSize_Males)
metaRegressionResult <- filter(metaRegressionResult, !is.na(p_value_aw))
metaRegressionResult$agingPValuesWithDirection <- rank(-1*metaRegressionResult$p_value_aw)
metaRegressionResult <- mutate(metaRegressionResult, agingPValuesWithDirection = if_else(medianEffect < 0, -1*agingPValuesWithDirection, agingPValuesWithDirection ))
sortedGenes <- arrange(metaRegressionResult, desc(agingPValuesWithDirection))$gene_symbol


if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { #assume it's already loaded - needs a fix to see if the variable is declared
} else {
  go_object <- as.list(org.Hs.egGO2ALLEGS)
  
  symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data='org.Hs.eg')
  
  #build GO sets for tmod -slow
  tmodNames <- data.frame()
  modules2genes <- list()
  goGroupName <- names(go_object)[1]
  showMethods(Term)
  
  goCount <- length(go_object)
  count <- 1
  for(goGroupName in names(go_object)) {
    # if (Ontology(goGroupName) == "CC") next(); #change in order to skip an ontology
    if (count %% 1000 == 0) print(paste(count, "of", goCount))
    count <- count + 1
    
    goGroup <- go_object[goGroupName]
      
    geneIDs <- unique(unlist(goGroup, use.names=F))  #discard evidence codes
    genesymbols <- unique(getSYMBOL(geneIDs, data='org.Hs.eg'))
    
    genesymbols <- intersect(genesymbols, sortedGenes)
    if (!(length(genesymbols) > 10 & length(genesymbols) < 200)) next();
    
    modules2genes[goGroupName] <- list(genesymbols)
    
    tmodNames <- rbind(tmodNames, data.frame(ID=goGroupName, Title = Term(goGroupName)))
  }
  geneSetsGO <- makeTmod(modules = tmodNames, modules2genes = modules2genes)
}
result <- tmodUtest(c(sortedGenes), mset=geneSetsGO, qval = 1, filter = T)
result <- tbl_df(result) %>% dplyr::select(Title, geneCount =N1,AUC,  P.Value, adj.P.Val, ID)
result <- mutate(rowwise(result), aspect = Ontology(ID))
subset(result, AUC < 0.5)
subset(result, AUC > 0.5)
result

#get the average effect sizes
if (targetVariable == "regression") {
  meanEffects <- NULL
  for(group in result$ID) {
    genes <- unlist(geneSetsGO$MODULES2GENES[group])
    meanEffect <- metaRegressionResult %>% filter(gene_symbol %in% genes) %>% ungroup() %>% summarize(ID = group, EffectSize_Females = mean(EffectSize_Females), EffectSize_Males = mean(EffectSize_Males))
    meanEffects <- bind_rows(meanEffects, meanEffect)
  }
  result <- inner_join(result, meanEffects, by="ID")
}
  
#write out
write.csv(result, paste0(gsub(".csv","","./data/metaA/MetaRegression_SexSpecific_2-17-17_MLS.csv"),".GOResults.",targetVariable,".csv"))

#Euler diagrams
source("./R/metaA/MakeEulerDiagram.R")
(plot <- getEulerDiagram(head(result, n=10)))
ggsave(plot = plot, filename=paste0(gsub(".csv","","./data/metaA/MetaRegression_SexSpecific_2-17-17_MLS.csv"),".top10Plot.",targetVariable,".pdf"), width=12, height=12)



###############
###############
###############
###############
###############
#plot a GO group in the median effect data
aw_result <- read_csv("./data/metaA/MDD-metaAR_8cohorts_Final.csv")
aw_result <- filter(aw_result, !is.na(REM_ALL_P))
effectSizeCols <- c("SYMBOL", "MD2_DLPFC_M.1", "MD1_ACC_M.1",  "MD2_ACC_M.1",  "MD1_AMY_M.1",  "MD2_DLPFC_F.1", "MD3_ACC_F.1", "MD2_ACC_F.1", "MD3_AMY_F.1")
aw_result <- aw_result[,effectSizeCols]

forRHeatmap <- as.data.frame(aw_result)
rownames(forRHeatmap) <- as.character(forRHeatmap[,"SYMBOL"])
forRHeatmap <- forRHeatmap[,-1]
termID <- 'GO:0070126'
median(as.matrix(forRHeatmap[,5:8]))
forRHeatmap <- forRHeatmap[unlist(geneSetsGO$MODULES2GENES[termID]),]
heatmap(as.matrix(forRHeatmap),margins=c(12,8), scale = "none", Colv = NA, main = Term(termID))
median(as.matrix(forRHeatmap[,5:8]))

###############
###############
###############
#check darmansis lists
tmodNames <- data.frame()
modules2genes <- list()
#need to set folder here
for(geneListFilename in list.files("/Users/lfrench/Desktop/results/CellTypesAging/data/GeneLists/", pattern = "Darmanis.*txt", full.names = T)) {
  print(geneListFilename)
  genesOfInterest <- read.csv(geneListFilename,header=F,stringsAsFactors = F)
  shortName <- gsub(".txt","",gsub(paste0(".*/"),"", geneListFilename))
  genesOfInterest$term <- shortName
  
  #already a human gene list
  modules2genes[shortName] <- list(genesOfInterest$V1)
  #print(intersect(sortedGenes, unlist(genesOfInterest$V1)))
  
  tmodNames <- rbind(tmodNames, data.frame(ID=shortName, Title = shortName))
}
geneSets <- makeTmod(modules = tmodNames, modules2genes = modules2genes)

result <- tmodUtest(sortedGenes, mset=geneSets, qval = 1.1, filter = F)
result <- tbl_df(result) %>% dplyr::select(Title, geneCount =N1,AUC,  P.Value, adj.P.Val) %>% mutate(Title = gsub("Darmanis.","",Title))
result
write.csv(result, paste0(gsub(".csv","","./data/metaA/MetaRegression_SexSpecific_2-17-17_MLS.csv"),".Darmanis.",targetVariable,".csv"))



##############################################
##############################################
##############################################
##############################################
#heatmap for a Darmanis list
cellType <- "Darmanis.Microglia"
cellType <- "Darmanis.OligoPrecusors"
cellType <- "Darmanis.Astrocytes"
cellType <- "Darmanis.Endothelial"
cellType <- "Darmanis.Oligo"
cellType <- "Darmanis.Microglia"
cellType <- "Darmanis.Neuron"

forRHeatmap <- as.data.frame(aw_result)
rownames(forRHeatmap) <- as.character(forRHeatmap[,"SYMBOL"])
forRHeatmap <- forRHeatmap[,-1]
forRHeatmap <- forRHeatmap[intersect(rownames(forRHeatmap), unlist(geneSets$MODULES2GENES[cellType])),]
heatmap(as.matrix(forRHeatmap),margins=c(12,8), scale = "none", Colv = NA, main = cellType)


######## compare values between old and new p-values
cellType <- "Darmanis.Microglia"
cellType <- "Darmanis.Oligo"
cellType <- "Darmanis.OligoPrecusors"
aw_result <- read_csv("./data/metaA/MDD-metaAR_8cohorts_Final.csv")
metaRegressionResult <- read_csv("./data/metaA/MetaRegression_SexSpecific_2-17-17_MLS.csv")
colnames(metaRegressionResult)[1] <- "SYMBOL"
goodGeneNames <- read_csv("./data/metaA/MDD-metaAR_8cohorts_Final.csv")
goodGeneNames$namesUpper <- toupper(goodGeneNames$SYMBOL)
metaRegressionResult <- inner_join(metaRegressionResult, dplyr::select(goodGeneNames, namesUpper,SYMBOL), by=c("SYMBOL"="namesUpper"))
metaRegressionResult$SYMBOL <- metaRegressionResult$SYMBOL.y

colnames(metaRegressionResult)[1] <- "SYMBOL"

effectSizeCols <- c("MD2_DLPFC_M.1", "MD1_ACC_M.1",  "MD2_ACC_M.1",  "MD1_AMY_M.1")
aw_result$medianEffect <- apply(aw_result[,effectSizeCols ],1,median) #median effect across the studies (split by sexes)

pValues <- left_join(dplyr::select(aw_result, SYMBOL, medianEffect), dplyr::select(metaRegressionResult, SYMBOL,EffectSize_Males, pvalue_Males, pvalue_Females, pvalue_MetaR_sex))
cellTypeOnly <- filter(pValues , SYMBOL %in% unlist(geneSets$MODULES2GENES[cellType])) %>% arrange(pvalue_MetaR_sex)
as.data.frame(cellTypeOnly)
plot(cellTypeOnly$medianEffect, cellTypeOnly$EffectSize_Males)
abline(a=0,b=1)

