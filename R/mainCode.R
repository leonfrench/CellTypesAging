#clear workspace
rm(list = ls())

#clear all libraries to prevent conflicts
#from mjaniec
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}
detachAllPackages()

library(doParallel)
library(readr)
library(reshape2)
library(homologene)
library(ROCR)
library(optparse)
library(dplyr)


option_list = list(
  make_option(c("-a", "--ageSource"), type="character", default=NULL, help="AW, Blood Mistry or MistryPH", metavar="character"),
  make_option(c("-s", "--cellSource"), type="character", default=NULL, help="Tasic or Zeisel", metavar="character"),
  make_option(c("-l", "--cellLevel"), type="character", default=NULL, help="TranscriptomicName (lowest level) or mainClass (highest level)", metavar="character"),
  make_option(c("-i", "--iterations"), type="integer", default=NULL, help=""),
  make_option(c("-c", "--cores"), type="integer", default=NULL, help="number of cores when multithreading")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (interactive()) { #set the variables manually if in Rstudio, for testing
  agingGeneSource <- "AW"
  agingGeneSource <- "metaA"
  #agingGeneSource <- "Blood"
  #agingGeneSource <- "Mistry"

    
  cellTypeSource <- "Tasic"
  #cellTypeSource <- "Zeisel"
  
  level <- "TranscriptomicName"
  #  level <- "mainClass"
  
  cores <- 3
  iterations <- 1000
  
} else if (!is.null(opt$ageSource) & !is.null(opt$cellSource) & !is.null(opt$iterations)) {
  agingGeneSource <- opt$ageSource
  cellTypeSource <- opt$cellSource
  iterations <- opt$iterations
  cores<-opt$cores
  level <- opt$cellLevel
} else {
  print_help(opt_parser)
  stop()
}

registerDoParallel(cores=cores)

#other params to add to options, currently hardcoded
filterForCoreTranscriptomeTypes <- T #tasic param
enrichmentThreshold <- 2

name <- paste(agingGeneSource, cellTypeSource,"emperical", "iter", iterations,level,sep=".")
print(name)

if (grepl("R$", getwd())) { #if in the source folder, move up
  setwd(paste0(getwd(),"/.."))
}
source("./R/AUCFunction.R") #load AUC function
print(getwd())

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
} else if (agingGeneSource == "MistryPH") {
  mistryAgeDown <- read_csv("./data/Mistry Age Tables/Supplementary Table 10.pHDown.csv")
  mistryAgeDown$direction <- -1
  mistryAgeUp <- read_csv("./data/Mistry Age Tables/Supplementary Table 10.pHUp.csv")
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
  
} else if (agingGeneSource == "metaA") { #meta depression score - temporary
  aw_result <- read_csv("./data/metaA/MDD-metaAR_8cohorts_Final.csv")
  effectSizeCols <- c("MD2_DLPFC_M.1", "MD1_ACC_M.1",  "MD2_ACC_M.1",  "MD1_AMY_M.1",  "MD2_DLPFC_F.1", "MD3_ACC_F.1", "MD2_ACC_F.1", "MD3_AMY_F.1")
  
  #male/female only
  effectSizeCols <- c("MD2_DLPFC_M.1", "MD1_ACC_M.1",  "MD2_ACC_M.1",  "MD1_AMY_M.1")
  #effectSizeCols <- c("MD2_DLPFC_F.1", "MD3_ACC_F.1", "MD2_ACC_F.1", "MD3_AMY_F.1")
  write.csv(effectSizeCols, file=paste0(outputFolder,"UsingMales.txt"),row.names=F)
  
  
  aw_result$medianEffect <- apply(aw_result[,effectSizeCols ],1,median) #median effect across the studies (split by sexes)
  
  #two methods for p-value REM_ALL_P or roP_OC_p
  #aw_result$p_value_aw <- aw_result$REM_ALL_P
  #males REM_M_P or females REM_F_P
  #aw_result$p_value_aw <- aw_result$REM_F_P
  aw_result$p_value_aw <- aw_result$REM_M_P
  
  aw_result <- dplyr::select(aw_result, gene_symbol = SYMBOL, p_value_aw, medianEffect)
  aw_result <- filter(aw_result, !is.na(p_value_aw))
  aw_result$agingPValuesWithDirection <- rank(-1*aw_result$p_value_aw)
  aw_result <- mutate(aw_result, agingPValuesWithDirection = if_else(medianEffect < 0, -1*agingPValuesWithDirection, agingPValuesWithDirection ))
  sortedGenes <- arrange(aw_result, desc(agingPValuesWithDirection))$gene_symbol
} else if (agingGeneSource == "metaASeneyMale") { #meta depression score - temporary
  aw_result <- read_csv("./data/metaA/MetaRegression_SexSpecific_2-17-17_MLS.csv")
  colnames(aw_result)[1] <- "SYMBOL"
  #fix symbols
  goodGeneNames <- read_csv("./data/metaA/MDD-metaAR_8cohorts_Final.csv")
  goodGeneNames$namesUpper <- toupper(goodGeneNames$SYMBOL)
  aw_result <- inner_join(aw_result, dplyr::select(goodGeneNames, namesUpper,SYMBOL), by=c("SYMBOL"="namesUpper"))
  aw_result$SYMBOL <- aw_result$SYMBOL.y

  aw_result$medianEffect <- aw_result$EffectSize_Males #median effect across the studies (split by sexes)
  aw_result$p_value_aw <- aw_result$pvalue_Males
  
  aw_result <- dplyr::select(aw_result, gene_symbol = SYMBOL, p_value_aw, medianEffect)
  aw_result <- filter(aw_result, !is.na(p_value_aw))
  aw_result$agingPValuesWithDirection <- rank(-1*aw_result$p_value_aw)
  aw_result <- mutate(aw_result, agingPValuesWithDirection = if_else(medianEffect < 0, -1*agingPValuesWithDirection, agingPValuesWithDirection ))
  sortedGenes <- arrange(aw_result, desc(agingPValuesWithDirection))$gene_symbol
} else if (agingGeneSource == "metaASeneyFemale") { #meta depression score - temporary
  aw_result <- read_csv("./data/metaA/MetaRegression_SexSpecific_2-17-17_MLS.csv")
  colnames(aw_result)[1] <- "SYMBOL"
  #fix symbols
  goodGeneNames <- read_csv("./data/metaA/MDD-metaAR_8cohorts_Final.csv")
  goodGeneNames$namesUpper <- toupper(goodGeneNames$SYMBOL)
  aw_result <- inner_join(aw_result, dplyr::select(goodGeneNames, namesUpper,SYMBOL), by=c("SYMBOL"="namesUpper"))
  aw_result$SYMBOL <- aw_result$SYMBOL.y

  aw_result$medianEffect <- aw_result$EffectSize_Females #median effect across the studies (split by sexes)
  aw_result$p_value_aw <- aw_result$pvalue_Females
  
  #code duplicated!
  aw_result <- dplyr::select(aw_result, gene_symbol = SYMBOL, p_value_aw, medianEffect)
  aw_result <- filter(aw_result, !is.na(p_value_aw))
  aw_result$agingPValuesWithDirection <- rank(-1*aw_result$p_value_aw)
  aw_result <- mutate(aw_result, agingPValuesWithDirection = if_else(medianEffect < 0, -1*agingPValuesWithDirection, agingPValuesWithDirection ))
  sortedGenes <- arrange(aw_result, desc(agingPValuesWithDirection))$gene_symbol
} else if (agingGeneSource == "metaASeneyMetaR") { #meta depression score - temporary
  aw_result <- read_csv("./data/metaA/MetaRegression_SexSpecific_2-17-17_MLS.csv")
  colnames(aw_result)[1] <- "SYMBOL"
  #fix symbols
  goodGeneNames <- read_csv("./data/metaA/MDD-metaAR_8cohorts_Final.csv")
  goodGeneNames$namesUpper <- toupper(goodGeneNames$SYMBOL)
  aw_result <- inner_join(aw_result, dplyr::select(goodGeneNames, namesUpper,SYMBOL), by=c("SYMBOL"="namesUpper"))
  aw_result$SYMBOL <- aw_result$SYMBOL.y

  aw_result$medianEffect <- aw_result$EffectSize_MetaR_sex #median effect across the studies (split by sexes)
  aw_result$p_value_aw <- aw_result$pvalue_MetaR_sex
  
  #code duplicated!
  aw_result <- dplyr::select(aw_result, gene_symbol = SYMBOL, p_value_aw, medianEffect)
  aw_result <- filter(aw_result, !is.na(p_value_aw))
  aw_result$agingPValuesWithDirection <- rank(-1*aw_result$p_value_aw)
  aw_result <- mutate(aw_result, agingPValuesWithDirection = if_else(medianEffect < 0, -1*agingPValuesWithDirection, agingPValuesWithDirection ))
  sortedGenes <- arrange(aw_result, desc(agingPValuesWithDirection))$gene_symbol
} else {
  print("ERROR no aging gene list source")
  stop()
}



#################################################################
#################################################################
#load the cell type data  
if(cellTypeSource == "Tasic") {
  source("./R/TasicProcessingCode.R")
} else if (cellTypeSource == "Zeisel") {
  source("./R/ZeiselProcessingCode.R")
} else {
  print("ERROR no cell type data source")
  stop()
}
#####################################
#the result should be a loaded cellTable (CellClassID) and rpkmDataCore (geneName long_name log1Expression CellClassID)
#####################################
getEnrichedHuman <- function(rpkmDataCore , aw_result_mouseFilter) {
  rpkmDataCore <- group_by(rpkmDataCore, CellClassID, geneName) #slow
  #average across cell types for each gene
  rpkmByType <- summarise(rpkmDataCore, log1Expression = mean(log1Expression))
  
  #scale/z-score - do z score after averaging in a cell type
  rpkmByType <- group_by(rpkmByType, geneName)
  rpkmByType <- mutate(rpkmByType, log1ExpressionZ = (log1Expression - mean(log1Expression)) / sd(log1Expression))
  
  #filter for one to many mappings? or just expand regardless? most mouse genes have only one human homolog so it's probably fine
  #genesWithOneHomolog <- (group_by(backgroundHumanGenes, mouseGene) %>% summarize(n = n_distinct(humanGene)) %>% filter(n == 1))$mouseGene
  
  #cut by enrichment 
  rpkmByType2x <- filter(rpkmByType, log1ExpressionZ > enrichmentThreshold) #enrichmentThreshold defined at top
  rpkmByType <- group_by(rpkmByType, geneName)
  stableGenes <- summarise(rpkmByType, log1ExpressionZ = max(log1ExpressionZ)) %>% filter(log1ExpressionZ < enrichmentThreshold)
  #add stable genes to rpkmByType2x
  stableGenes$CellClassID <- "CelltypeNonSpecific"
  
  range((group_by(rpkmByType2x, CellClassID) %>% summarise(genes = n_distinct(geneName)))$genes)
  
  #################################################################
  #get average amount of non-specitif genes
  averageSize <- median((group_by(rpkmByType2x, CellClassID) %>% summarise(genes = n_distinct(geneName)))$genes)
  averageSize <- min(nrow(stableGenes), averageSize) #in cases when the average size is lower than population size
  stableGenesSample <- dplyr::sample_n(stableGenes,averageSize)
  stableGenesSample$CellClassID <- "CelltypeNonSpecificSample"
  
  rpkmByType2x <- bind_rows(dplyr::select(rpkmByType2x,CellClassID,log1ExpressionZ, geneName), stableGenes)
  #rpkmByType2x <- bind_rows(rpkmByType2x, stableGenesSample)
  #################################################################
  
  #convert to human genes
  backgroundHumanGenes <- tbl_df(mouse2human(rpkmByType2x$geneName))
  backgroundHumanGenes <- filter(backgroundHumanGenes, humanGene %in% aw_result_mouseFilter$gene_symbol)
  
  #convert to human genes - remove those with no ortholog
  rpkmByType2xHuman <- inner_join(rpkmByType2x, backgroundHumanGenes, by=c("geneName" = "mouseGene"))
  
  #now the gene lists are ready
  #join aging list and cell type list
  rpkmByType2xHuman <- dplyr::inner_join(rpkmByType2xHuman, dplyr::select(aw_result_mouseFilter, gene_symbol, agingPValuesWithDirection), by=c("humanGene"="gene_symbol"))
  rpkmByType2xHuman
}

#filter awresult for genes that are reachable from human
getFilteredAWresult <- function(aw_result, mouseGenes) {
  backgroundHumanGenes <- tbl_df(mouse2human(mouseGenes))
  aw_result_mouseFilter  <- filter(aw_result, gene_symbol %in% backgroundHumanGenes$humanGene)
  aw_result_mouseFilter
}

#filter out aw_result human genes that cannot be reached from the mouse genes, only do once
aw_result_mouseFilter <- getFilteredAWresult(aw_result, unique(rpkmDataCore$geneName)) 
print(paste("Number of genes in ranking after filtering:" , nrow(aw_result_mouseFilter)))
#calls above function
rpkmByType2xHuman <- getEnrichedHuman(rpkmDataCore, aw_result_mouseFilter) 

print(left_join(rpkmByType2xHuman %>% group_by(CellClassID) %>% summarise(genes = n_distinct(humanGene)), rpkmDataCore %>% group_by(CellClassID) %>% summarise(CellCount = n_distinct(long_name)), by="CellClassID"),n=51)

#iterate cell type lists here
geneCellResults <- data.frame( CellClassID=character(), geneCount= numeric(), p=numeric(), stringsAsFactors=F)
resultposition <- 0

for(targetCellType in unique(rpkmByType2xHuman$CellClassID)) {
  resultposition <- resultposition + 1
  resultsOfInterest <- subset(rpkmByType2xHuman, CellClassID == targetCellType)
  #print(targetCellType)
  genesOfInterest <- resultsOfInterest$humanGene
  labelMask <- aw_result_mouseFilter$gene_symbol %in% genesOfInterest
  auroc <- auroc_analytic(rank(aw_result_mouseFilter$agingPValuesWithDirection), labelMask)
  wilcoxP <- wilcox.test(aw_result_mouseFilter[!labelMask,]$agingPValuesWithDirection, aw_result_mouseFilter[labelMask, ]$agingPValuesWithDirection)$p.value
  
  geneCellResults[resultposition, "geneCount"] <- length(genesOfInterest)
  geneCellResults[resultposition, "CellClassID"] <- targetCellType
  geneCellResults[resultposition, "AUROC"] <- auroc
  geneCellResults[resultposition, "p"] <- wilcoxP
  
}
geneCellResults <- left_join(geneCellResults, rpkmDataCore %>% group_by(CellClassID) %>% summarise(CellCount = n_distinct(long_name)), by="CellClassID")

geneCellResults$p.adjusted <- p.adjust(geneCellResults$p)
geneCellResults <- geneCellResults[order(geneCellResults$p),]


##############
# write out

outputFolder <- paste0("./results/",name, ".",format(Sys.time(), "%s"), "/")
dir.create(outputFolder)
save(file=paste0(outputFolder,"rpkmByType2xHuman.RData"),rpkmByType2xHuman)
save(file=paste0(outputFolder,"aging_genes_aw_fisher_data.RData"),aw_result_mouseFilter)
write.csv(geneCellResults, file=paste0(outputFolder,"CellEnrichment.csv"),row.names=F)

dir.create(paste0(outputFolder, "filterForCoreTranscriptomeTypes.",filterForCoreTranscriptomeTypes))
dir.create(paste0(outputFolder, "enrichmentThreshold.",enrichmentThreshold))
dir.create(paste0(outputFolder, "aw_result.human gene count.",nrow(aw_result_mouseFilter)))



#############################################################################################
#############################################################################################
#create empirical p-values
#############################################################################################
geneCellResults <- as.data.frame(geneCellResults)
rownames(geneCellResults) <- geneCellResults$CellClassID
start.time <- Sys.time()
cellTypeList <- unique(rpkmByType2xHuman$CellClassID)
print("Starting emperical runs")
empericalAUCs = foreach(i=1:iterations, .combine=data.frame) %dopar% {  
  thisResult <- as.data.frame(geneCellResults[,c()])
  #print(paste(iteration, ", time: " , Sys.time() - start.time))
  if (i %% 100 == 0) print(i)
  #shuffle the core
  cellTableShuffled <- cellTable
  cellTableShuffled$CellClassID <- sample(cellTableShuffled$CellClassID)
  rpkmDataCoreShuffled <- rpkmDataCore %>% ungroup() %>% dplyr::select(-CellClassID)
  rpkmDataCoreShuffled <- inner_join(rpkmDataCoreShuffled, dplyr::select(cellTableShuffled,CellClassID,long_name), by="long_name")
  
  rpkmByType2xHuman <- getEnrichedHuman(rpkmDataCoreShuffled, aw_result_mouseFilter)
  
  #iterate cell type lists here
  for(targetCellType in cellTypeList) {
    resultsOfInterest <- subset(rpkmByType2xHuman, CellClassID == targetCellType)
    genesOfInterest <- resultsOfInterest$humanGene
    labelMask <- aw_result_mouseFilter$gene_symbol %in% genesOfInterest
    auroc <- auroc_analytic(rank(aw_result_mouseFilter$agingPValuesWithDirection), labelMask)
    thisResult[targetCellType, "AUROC"] <- auroc
  }
  thisResult
}
Sys.time() - start.time
print(paste(iterations, "iterations executed"))

write.csv(empericalAUCs, file=paste0(outputFolder,"EmpericalAUCsBackup.csv"),row.names=F)
print("wrote emperical AUCs to file")

empericalAUCs$CellClassID <- rownames(empericalAUCs)
meltedEmperical <- tbl_df(melt(tbl_df(empericalAUCs)))

meltedEmperical <- meltedEmperical[complete.cases(meltedEmperical),]
meltedEmperical <- rename(meltedEmperical, run=variable)

meltedEmperical <- group_by(meltedEmperical, CellClassID)
meltedEmperical <- mutate(meltedEmperical, AUROCMinusMean = value - mean(value))

geneCellResults <- inner_join(geneCellResults, summarise(meltedEmperical, meanEmpAUC = mean(value)))

geneCellResults <- mutate(geneCellResults, AUCDiff = AUROC - meanEmpAUC)

pValues <- inner_join(meltedEmperical, dplyr::select(geneCellResults, CellClassID, AUCDiff)) %>% group_by(CellClassID) %>% dplyr::summarise(countMoreExtreme = sum(abs(AUROCMinusMean) > abs(AUCDiff)), iterations = n(),by = "CellClassID")
pValues <- mutate(pValues, countMoreExtreme= ifelse(countMoreExtreme == 0, 1, countMoreExtreme)) %>% mutate(empericalP = countMoreExtreme/iterations)
geneCellResults <- inner_join(geneCellResults, dplyr::select(pValues, CellClassID, empericalP), by = "CellClassID") %>% dplyr::arrange(empericalP, desc(abs(AUCDiff)))

#do FDR - same as BH..
meltedEmperical <- group_by(meltedEmperical, CellClassID)
meltedEmperical <- mutate(meltedEmperical, pValue = rank(abs(AUROCMinusMean))/n())
minPValues <- summarize(group_by(meltedEmperical, run), minP = min(pValue))
geneCellResults <- group_by(geneCellResults, CellClassID)
geneCellResults <- dplyr::mutate(geneCellResults, empericalFDRP = sum(empericalP >= c(minPValues$minP)) / nrow(minPValues) )
geneCellResults$p.adjustedR <- p.adjust(geneCellResults$empericalP)


write.csv(geneCellResults, file=paste0(outputFolder,"CellEnrichmentEmpericalP.csv"),row.names=F)
write.csv(dplyr::select(geneCellResults,CellClassID, CellCount, geneCount, AUROC, meanEmpAUC, empericalP, empericalFDRP), file=paste0(outputFolder,"CellEnrichmentEmpericalPMinCols.csv"),row.names=F)

print("Resetting data")
#set back to the real data
rpkmByType2xHuman <- getEnrichedHuman(rpkmDataCore, aw_result_mouseFilter)
print("Done execution of cell type stage")

##############################################################################################
startMain <- Sys.time()
source("./R/CustomGeneLists.R")
source("./R/GOGroupByCellTypeTests.R") #cell type by GO group tests
Sys.time() - startMain
print(Sys.time() - startMain)
