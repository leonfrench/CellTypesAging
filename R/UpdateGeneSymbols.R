#this updates mouse symbols in the Tasic cell type data
#this is run once code to output a text file that is loaded by mainCode

library(dplyr)
library(readr)
library(homologene)

load("/Users/lfrench/Desktop/results/CellTypesAging/AgingFullResults_UniqueGenes.RData",verbose=T)
rpkmData <- read_csv("/Users/lfrench/Desktop/results/CellTypesAging/Allen data/case study/genes_rpkm.csv")
colnames(rpkmData)[1] <- "geneName"

detach("package:mygene", unload=TRUE)
detach("package:GenomicFeatures", unload=TRUE)
detach("package:AnnotationDbi", unload=TRUE)

backgroundHumanGenes <- mouse2human(rpkmData$geneName)

#human genes accessable from both datasets
backgroundHumanGeneList <- intersect(backgroundHumanGenes$humanGene, rownames(p.dat.unique.genes))
missedMouseGenes <- setdiff(rpkmData$geneName, backgroundHumanGenes$mouseGene)
length(missedMouseGenes)

library(mygene)
#missed mouse genes - due to bad symbol or not in homologene
foundNewSymbol <- 0
multipleNewSymbols <- 0

newGeneSybols <- data.frame(  newSymbol=character(), stringsAsFactors=F)

for(missedGene in missedMouseGenes) {
  #for(missedGene in head(missedMouseGenes,n=500)) {
  #print(missedGene)
  myGeneResult <- mygene::query(q=missedGene,species="mouse")
  #print(paste("hit:",myGeneResult$hits$symbol))
  
  
  hitSymbols <- unique(myGeneResult$hits$symbol)
  hitSymbols <- setdiff(hitSymbols, missedGene)
  hits <- length(hitSymbols)
  if (hits >1 && "symbol" %in% colnames(myGeneResult$hits)) {
    print(missedGene)
    print(paste("to many hits:", hits))
    print(myGeneResult)
    multipleNewSymbols <- multipleNewSymbols+1
  }
  #if the doesn't match and it's only one remaining
  if (hits == 1 && myGeneResult$hits$symbol[1] != missedGene) {
    
    newGeneSybols[missedGene,"newSymbol"] = hitSymbols[1]

    foundNewSymbol <- foundNewSymbol + 1
  }
}
print(paste("multiple new symbols for", multipleNewSymbols)) #85
print(paste("updated symbols for", foundNewSymbol)) #1090

#unload do to select overloading
detach("package:mygene", unload=TRUE)
detach("package:GenomicFeatures", unload=TRUE)
detach("package:AnnotationDbi", unload=TRUE)

write.csv(newGeneSybols, "/Users/lfrench/Desktop/results/CellTypesAging/NewGeneSymbols.csv")
