library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)
library(annotate)
library(GO.db)
library(hash)

source("./R/AUCFunction.R") #load AUC function

#show version information
org.Hs.eg.db
go_object <- as.list(org.Hs.egGO2ALLEGS)

#filter for all genes with any GO group!
symbolsInGO <- getSYMBOL(unique(unlist(go_object)), data='org.Hs.eg')

#create a hash table that links GO group to gene symbol and name
goTable <- hash()
goNames <- hash()
system.time(
  for(goGroupName in names(go_object)) {
    name <- Term(goGroupName)
    goGroup <- go_object[goGroupName]
    geneIDs <- unique(unlist(goGroup, use.names=F))  #discard evidence codes
    genesymbols <- unique(getSYMBOL(geneIDs, data='org.Hs.eg'))
    goNames[[goGroupName]] <- name
    goTable[[goGroupName]] <- genesymbols
  }
)
print("Done creating GO Hash")
start <- Sys.time()


#load rpkm data

load("/Users/lfrench/Desktop/results/CellTypesAging/results/AW.Tasic.emperical.iter.1.TranscriptomicName.1479920613/rpkmByType2xHuman.RData",v=T)
load("/Users/lfrench/Desktop/results/CellTypesAging/results/AW.Tasic.emperical.iter.1.TranscriptomicName.1479920613/aging_genes_aw_fisher_data.RData",v=T)
#get genes for a specific group and cell type
targetCellType <- "Sst Tacstd2"
targetCellType <- "Astro Aqp4"
targetCellType <- "Sst Cdk6"
goGroupName <- "GO:0007268" #synaptic transmission
#goGroupName <-"GO:0007267" #cell-cell signalling
####################################################

load("/Users/lfrench/Desktop/results/CellTypesAging/results/Mistry.Zeisel.emperical.iter.1.TranscriptomicName.1479938797/aging_genes_aw_fisher_data.RData",v=T)
load("/Users/lfrench/Desktop/results/CellTypesAging/results/Mistry.Zeisel.emperical.iter.1.TranscriptomicName.1479938797/rpkmByType2xHuman.RData",v=T)
targetCellType <- "Int2"
goGroupName <- "GO:0007268" #synaptic transmission

####################################################

#filter for genes in GO
aw_result_mouseFilter <- aw_result_mouseFilter %>% dplyr::filter(gene_symbol %in% symbolsInGO)
aw_result_mouseFilter$rank <- rank(-1*aw_result_mouseFilter$agingPValuesWithDirection, ties.method = "random")



#need actual rank!
rpkmByType2xHuman <- inner_join(rpkmByType2xHuman, dplyr::select(aw_result_mouseFilter, gene_symbol, rank), by=c("humanGene" = "gene_symbol"))

goFullName <- paste(toupper(substring(goNames[[goGroupName]], 1,1)), substring(goNames[[goGroupName]], 2),sep="", collapse=" ")

resultsOfInterestGOGroup <- dplyr::filter(rpkmByType2xHuman, humanGene %in% goTable[[goGroupName]])
resultsOfInterestGOGroup$name <- paste(goFullName)
resultsOfInterestCellType <- subset(rpkmByType2xHuman, CellClassID == targetCellType) %>% dplyr::filter(humanGene %in% symbolsInGO)
resultsOfInterestCellType$name <- paste(targetCellType)
resultsOfInterestCellTypeAndGO <- dplyr::filter(resultsOfInterestCellType, humanGene %in% goTable[[goGroupName]])
resultsOfInterestCellTypeAndGO$name <- paste(goFullName, "and", targetCellType)

forPlot <- rbind(resultsOfInterestCellType, resultsOfInterestGOGroup,resultsOfInterestCellTypeAndGO)

resultsOfInterestGOGroup


(rasterPlot <- ggplot(forPlot, aes(x = rank)) + 
  geom_blank() + 
  geom_vline(aes(xintercept=rank), size=0.2) +
  theme_bw()+coord_cartesian(expand=F) +
  facet_grid(name ~ ., switch="y") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
)
  
width <- 0.5
#first bar - narrow - GO group
resultsOfInterestGOGroup$internalRank <- rank(resultsOfInterestGOGroup$rank) + (median(resultsOfInterestGOGroup$rank)-nrow(resultsOfInterestGOGroup)/2) #center internal rank
resultsOfInterestGOGroup$x <- resultsOfInterestGOGroup$internalRank
resultsOfInterestGOGroup$y <- 2 - width/2
resultsOfInterestGOGroup$alpha = 0.5
forPlotGO <- resultsOfInterestGOGroup
resultsOfInterestGOGroup$y <- 2 + width/2
resultsOfInterestGOGroup$alpha = 0.5
forPlotGO <- bind_rows(forPlotGO, resultsOfInterestGOGroup)

#second bar - genome wide
resultsOfInterestGOGroup$x <- resultsOfInterestGOGroup$rank
resultsOfInterestGOGroup$y <- 1 - width/2
resultsOfInterestGOGroup$alpha = 0.5
forPlotGO <- bind_rows(forPlotGO, resultsOfInterestGOGroup)
resultsOfInterestGOGroup$y <- 1 + width/2
resultsOfInterestGOGroup$alpha = 0.15
forPlotGO <- bind_rows(forPlotGO, resultsOfInterestGOGroup)

forPlotGO <- arrange(forPlotGO, y)

#third bar - genome wide
resultsOfInterestCellType$x <- resultsOfInterestCellType$rank
resultsOfInterestCellType$y <- 0 - width/2
resultsOfInterestCellType$alpha = 0.5
forPlotCell <- resultsOfInterestCellType
resultsOfInterestCellType$y <- 0 + width/2
resultsOfInterestCellType$alpha = 0.5
forPlotCell <- bind_rows(forPlotCell, resultsOfInterestCellType)

#fourth bar - internal
resultsOfInterestCellType$internalRank <- rank(resultsOfInterestCellType$rank) + (median(resultsOfInterestCellType$rank)-nrow(resultsOfInterestCellType)/2) #center internal rank
resultsOfInterestCellType$x <- resultsOfInterestCellType$internalRank
resultsOfInterestCellType$y <- -1 - width/2
resultsOfInterestCellType$alpha = 0.5
forPlotCell <- bind_rows(forPlotCell, resultsOfInterestCellType)
resultsOfInterestCellType$y <- -1 + width/2
resultsOfInterestCellType$alpha = 0.15
forPlotCell <- bind_rows(forPlotCell, resultsOfInterestCellType)

forPlotCell <- arrange(forPlotCell, y)

forPlotIntersect <- bind_rows(forPlotGO, forPlotCell) %>% filter(humanGene %in% resultsOfInterestCellTypeAndGO$humanGene) %>% arrange(y)

#compute the 4 AUC's
(GOVrGenome <- signif(auroc_analytic(rank(aw_result_mouseFilter$agingPValuesWithDirection), aw_result_mouseFilter$gene_symbol %in% goTable[[goGroupName]] ),2))
miminalResultsOfInterestGOGroup <- ungroup(resultsOfInterestGOGroup) %>% dplyr::select( humanGene, agingPValuesWithDirection) %>% distinct()
(GOVrCellType <- signif(auroc_analytic(rank(miminalResultsOfInterestGOGroup$agingPValuesWithDirection), miminalResultsOfInterestGOGroup$humanGene %in% resultsOfInterestCellTypeAndGO$humanGene),2))
(wilcoxP <- wilcox.test(miminalResultsOfInterestGOGroup$agingPValuesWithDirection[miminalResultsOfInterestGOGroup$humanGene %in% resultsOfInterestCellTypeAndGO$humanGene],miminalResultsOfInterestGOGroup$agingPValuesWithDirection[!(miminalResultsOfInterestGOGroup$humanGene %in% resultsOfInterestCellTypeAndGO$humanGene)])$p.value)  #, alternative = "less"? force same direction as whole cell type?


(CellTypeGenome <- signif(auroc_analytic(rank(aw_result_mouseFilter$agingPValuesWithDirection), aw_result_mouseFilter$gene_symbol %in% resultsOfInterestCellType$humanGene ),2))
miminalResultsOfInterestCellType <- ungroup(resultsOfInterestCellType) %>% dplyr::select( humanGene, agingPValuesWithDirection) %>% distinct()
(CellTypeVrGO <- signif(auroc_analytic(rank(miminalResultsOfInterestCellType$agingPValuesWithDirection), miminalResultsOfInterestCellType$humanGene %in% resultsOfInterestCellTypeAndGO$humanGene),2))

showAll= F

if (!showAll) {

  forPlotGO <- ungroup(forPlotGO) %>% filter(y < 1.5)
  forPlotGO$alpha <- 0.5
  ungroup(forPlotCell) %>% dplyr::select(y, alpha) %>% distinct()
  
  forPlotIntersect <- filter(forPlotIntersect, y < 1.5)
  ungroup(forPlotIntersect) %>% dplyr::select(y, alpha) %>% distinct()
}
forPlotGO <- dplyr::select(forPlotGO, -log1ExpressionZ, -CellClassID) %>% distinct()

(linePlot <- ggplot(forPlotGO, aes(x = x, y = y, group = humanGene)) + 
  geom_path(aes(alpha=alpha)) + geom_path(data = forPlotCell, aes(alpha=alpha)) +
  geom_path(data= forPlotIntersect, color="blue", aes(alpha=alpha*2)) +
  theme_bw() +
  annotate("text", x=min(resultsOfInterestCellType$internalRank), y=-1, label=paste(goFullName,"within", targetCellType,"enriched genes  \n  AUROC:", CellTypeVrGO, " "), hjust = 1) +
  #annotate("text", x=0, y=1-width/2, label=paste(" ",goFullName, "genes\n   AUROC:", GOVrGenome), hjust = 0,vjust=1.3) +
  #annotate("text", x=0, y=0+width/2, label=paste(" ", targetCellType, "enriched genes\n   AUROC:", CellTypeGenome), hjust = 0,vjust=-0.3, fill = "green") +
  geom_label(aes(x = 0, y = 1-width/2, label = paste(" ",goFullName, "genes\n   AUROC:", GOVrGenome)), hjust = 0,vjust=1.3, fill = "white", label.size=NA, label.padding=unit(0.25, "lines")) +
  geom_label(aes(x = 0, y = 0+width/2, label = paste(" ", targetCellType, "enriched genes\n   AUROC:", CellTypeGenome)), hjust = 0,vjust=-0.3, fill = "white", label.size=NA, label.padding=unit(0.25, "lines")) +
  annotate("text", x=min(resultsOfInterestCellTypeAndGO$rank), y=0.5, label=" Overlaping genes ", hjust = ifelse(min(resultsOfInterestCellTypeAndGO$rank) > 2000,1,0), color="blue") +
  coord_cartesian(expand=F) +
  theme(#axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        #panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        axis.text=element_text(size=12.5),
        axis.title=element_text(size=13,face="bold")) +
  scale_x_continuous(name = paste0("Age-associated gene ranking (",nrow(aw_result_mouseFilter)," genes)"), breaks= c(850, nrow(aw_result_mouseFilter)-1000), labels = c("Up-regulated", "Down-regulated")) +
  geom_path(data= forPlotIntersect, color="blue", aes(alpha=alpha/2)) #add a light background on the blue for text overlaps
)

if(showAll) {
  (linePlot <- linePlot + annotate("text", x=min(resultsOfInterestGOGroup$internalRank), y=2, label=paste(targetCellType,"enriched within", goFullName,"genes, AUROC:", GOVrCellType, " "), hjust = 1))
} 
linePlot 


ggsave(paste0("/Users/lfrench/Desktop/results/CellTypesAging/Figures/projectionPlot.",targetCellType,".and.",goFullName,".showAll.", showAll, ".pdf"), linePlot, width = 9.5, height = 8)

