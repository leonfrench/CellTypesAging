######## GSEA Code for Leon ####
#> Sys.Date()
#[1] "2016-08-24"

setwd('/Users/lfrench/Desktop/results/CellTypesAging/R/fromCharles/')
#setwd('/Users/matianzhou/Documents/Pitt Spring 2016/Etienne/Singe Cell RNAseq/FromLeon160720')
load('AgingFullResults_UniqueGenes.RData', verbose = T)
## include aw result, individual effect size and  pvalue
load('enrichedGenesPerCellType.rdata')

###### real part########
cell.type <- unique(genesPerCellTypeReal[,1])
enriched_genes <- vector("list",length(cell.type))
for (i in 1:length(cell.type)){
  enriched_genes[[i]] <- genesPerCellTypeReal[which(genesPerCellTypeReal[,1]==cell.type[i]),2]
}
names(enriched_genes) <- cell.type

###### negative control########
cell.type <- unique(genesPerCellTypeRandom[,1]) 
enriched_genes <- vector("list",length(cell.type))
for (i in 1:length(cell.type)){
  enriched_genes[[i]] <- genesPerCellTypeRandom[which(genesPerCellTypeRandom[,1]==cell.type[i]),2]
}
names(enriched_genes) <- cell.type

########## prepare the input data - filters for same direction in both regions
unique.aging.genes <- rownames(es.dat.unique.genes)[which(sign(es.dat.unique.genes[,1])*sign(es.dat.unique.genes[,2])==1)]
unique.up.genes <- unique.aging.genes[which(sign(es.dat.unique.genes[unique.aging.genes,1])==1)]
unique.down.genes <- unique.aging.genes[which(sign(es.dat.unique.genes[unique.aging.genes,1])== -1)]
dat <- data.frame(geneNames = c(unique.up.genes,unique.down.genes),
                  direction=c(rep(1,length(unique.up.genes)),rep(-1,length(unique.down.genes))),
                  agingP = c(-log10(aw.out[unique.up.genes,4]),-log10(aw.out[unique.down.genes,4])), 
                  cellMarker=rep("BG",length(unique.up.genes)+length(unique.down.genes)),stringsAsFactors = F)
rownames(dat) <- dat[,1]
dat$agingPValuesWithDirection <- c(dat[dat$direction >0, "agingP"], -dat[dat$direction <0, "agingP"]) 

#---------------------------- Leon's changes start here
B <- 100000 ## number of permutations
#library(devtools)
#install_github("ctlab/fgsea")
library(data.table)
library(fgsea)
library(ggplot2)

dat <- dat[order(dat$agingPValuesWithDirection, decreasing = T),]
geneStatsForGSEA <- setNames(dat$agingPValuesWithDirection,rownames(dat))
head(geneStatsForGSEA)
tail(geneStatsForGSEA)

fgseaRes <- fgsea(pathways = enriched_genes, 
                  stats = geneStatsForGSEA,
                  minSize=0,
                  maxSize=2000,
                  nperm=B, 
                  gseaParam = 1)

fgseaRes <- fgseaRes[order(fgseaRes$pval),]
fgseaRes$leadingEdgeSize <- sapply(fgseaRes$leadingEdge, function(x) { length(unlist(x))})

head(fgseaRes,n=20) #note the relationship between nMoreExtreme and pval is not the expected 

print(paste("Cell types with adjusted p < 0.05:" ,sum(fgseaRes$padj < 0.05)))


############################ trying the plotting functions
plotEnrichment(enriched_genes[["CelltypeNonSpecific"]],
               geneStatsForGSEA) + labs(title="CelltypeNonSpecific")

topPathwaysUp <- fgseaRes[ES > 0 & padj < 0.05, ][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0 & padj < 0.05, ][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(enriched_genes[topPathways], geneStatsForGSEA, fgseaRes, 
              gseaParam = 1)
