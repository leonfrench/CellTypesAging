#devtools::install_github("sachsmc/plotROC")
#devtools::install_github("hadley/ggplot2")
#devtools::install_github("baptiste/gridextra") #use the dev version
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

library(scales)
library(plotROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(cowplot)

load("/Users/lfrench/Desktop/results/CellTypesAging/server results/AW.Tasic.emperical.iter.10000.TranscriptomicName.1479930034/aging_genes_aw_fisher_data.RData",v=T)
load("/Users/lfrench/Desktop/results/CellTypesAging/server results/AW.Tasic.emperical.iter.10000.TranscriptomicName.1479930034/rpkmByType2xHuman.RData", v=T)


ranking <- dplyr::select(aw_result_mouseFilter, gene_symbol, agingPValuesWithDirection)
#show only the top results
rpkmByType2xHuman <- rpkmByType2xHuman %>% mutate(CellClassID = if_else(CellClassID=="CelltypeNonSpecific", "Non Specific",CellClassID))

typesToUse <- c("Oligo Opalin","Astro Aqp4", "Oligo 9630013A20Rik", "L5a Batf3","L2/3 Ptgs2","L4 Arf5","Non Specific")
geneToClass <-  ungroup(rpkmByType2xHuman) %>% dplyr::select(gene_symbol=humanGene,CellClassID) %>% distinct() 


geneToClassAUC <- left_join(ranking, geneToClass) %>% spread(key=CellClassID, value=CellClassID) %>% gather(CellClassID, present, -gene_symbol, -agingPValuesWithDirection) %>% filter(CellClassID != "<NA>")
geneToClassAUC <- geneToClassAUC %>% mutate(present = if_else(is.na(present), 0, 1))
geneToClassAUC$dummy <- "True positive fraction"

length(unique(geneToClassAUC$gene_symbol))

geneToClassAUC.Others <- filter(geneToClassAUC, !(CellClassID %in% typesToUse))
geneToClassAUC.Others$CellClassID <- factor(geneToClassAUC.Others$CellClassID, levels = c(unique(geneToClassAUC.Others$CellClassID), "Cell Types with q > 0.05"))

geneToClassAUC <- filter(geneToClassAUC, CellClassID %in% typesToUse)
z <- geneToClassAUC[1,] #add one datapoint to help with ordering of legend
z$CellClassID <- "Cell Types with q > 0.05"
geneToClassAUC <- rbind(z, geneToClassAUC)
geneToClassAUC$CellClassID <- factor(geneToClassAUC$CellClassID, levels = c(typesToUse[1:3], "Cell Types with q > 0.05", typesToUse[4:7]))

(AUCPlot <- (basicplot <- ggplot(geneToClassAUC, aes(d = present, m = agingPValuesWithDirection, color=CellClassID)) +
               ylab("") + 
               geom_roc(n.cuts=0) + 
               style_roc() + coord_cartesian(expand=F)  +
               theme(legend.position = c(1,0), legend.justification = c(1, 0), legend.background= element_rect(fill = "transparent", colour = "transparent"), plot.margin=unit(c(.5,.5,.5,.5),"cm")) + 
               labs(color='Transcriptomic type') +
               facet_grid(dummy ~ ., switch="y") +
               ylab("") +
               theme(strip.background = element_blank(), strip.placement = "inside", strip.text = element_blank()) 
))

for(targetCell in unique(geneToClassAUC.Others$CellClassID)) {
  print(targetCell)
  AUCPlot <- AUCPlot + geom_roc(data=filter(geneToClassAUC.Others,CellClassID==targetCell) %>% mutate(CellClassID="Cell Types with q > 0.05"), n.cuts=0, linealpha=.25,pointalpha=.25)
}


#color generator from John Colby
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors <- gg_color_hue(7)
colors <- c(colors[1:3], "grey", colors[4:7])
(AUCPlot <- AUCPlot + scale_color_manual(values=colors) + geom_roc(n.cuts=0)) #overplot the main curves again


#, strip.text = element_text("True positive fraction")

calc_auc(basicplot)$AUC #just to check

geneToClassAUC <- geneToClassAUC %>% filter(CellClassID != "Cell Types with q > 0.05") #remove the point for forcing legend

geneToClassAUC <- geneToClassAUC %>% group_by(CellClassID) %>% mutate(rank = rank(-1*agingPValuesWithDirection, ties.method = "random"))
#plot as ranks
geneToClassAUC$CellClassID <- gsub(" ", "\n", geneToClassAUC$CellClassID)
geneToClassAUC$CellClassID <- factor(geneToClassAUC$CellClassID, levels = gsub(" ", "\n",typesToUse))

(rasterPlot <- ggplot(geneToClassAUC, aes(x = rank, y = present, color= CellClassID)) + 
  geom_blank() + 
  geom_vline(data = filter(geneToClassAUC, present == 1), aes(xintercept=rank),color="black", size=0.07) + 
  theme_bw()+coord_cartesian(expand=F) +
  ylab("Transcriptomic cell type") + 
  facet_grid(CellClassID ~ ., switch="y") + #, switch = "both"
  theme(strip.background = element_blank(), strip.placement = "inside") + #, strip.text.y = element_text(angle = 180)) +
  theme(axis.title.y = element_blank(),  axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.ticks.x=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(name = paste0("Age-associated gene ranking (",length(unique(geneToClassAUC$gene_symbol))," genes)"), breaks= c(min(geneToClassAUC$rank)+700, max(geneToClassAUC$rank)-700), labels = c("Up-regulated", "Down-regulated")))

#check col counts for lining up diagrams
ncol(ggplotGrob(rasterPlot))
ncol(ggplotGrob(AUCPlot))

(bothPlots <- plot_grid(AUCPlot, rasterPlot, labels = c("A", "B"), nrow = 2, align = "v", rel_heights=c(1,0.8)))
save_plot("/Users/lfrench/Desktop/results/CellTypesAging/Figures/AUCandRaster.plusGrey.pdf", bothPlots, base_width = 8, base_height = 12)
