library(magrittr)
#run plotROC to setup data

geneToClassAUC$CellClassID <- gsub("\n"," ", geneToClassAUC$CellClassID)
geneToClassAUC %<>% filter(CellClassID != "Non Specific")
geneToClassAUC %<>% filter(CellClassID != "L5a Batf3")
geneToClassAUC %<>% filter(CellClassID != "L2/3 Ptgs2")


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

(rasterPlot <- ggplot(geneToClassAUC, aes(x = rank, y = present, color= CellClassID)) + 
   geom_blank() + 
   geom_vline(data = filter(geneToClassAUC, present == 1), aes(xintercept=rank,color= CellClassID),size=0.1) + 
   theme_bw()+coord_cartesian(expand=F) +
   theme(axis.title.y = element_blank(),  axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.ticks.x=element_blank()) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   scale_x_continuous(name = "", breaks= c(min(geneToClassAUC$rank)+700, max(geneToClassAUC$rank)-700), labels = c("Up-regulated", "Down-regulated")) +
  scale_colour_manual(name="", values = c(gg_color_hue(7)[2],gg_color_hue(7)[6],gg_color_hue(7)[3],gg_color_hue(7)[1]),  guide = guide_legend(override.aes = list(size = 5))) +
  theme(legend.position="bottom") 
)



(AUCPlot <- (basicplot <- ggplot(geneToClassAUC, aes(d = present, m = agingPValuesWithDirection, color=CellClassID)) +
               ylab("") + 
               geom_roc(n.cuts=0) + 
               style_roc() + coord_cartesian(expand=F)  +
               theme(legend.position = c(1,0), legend.justification = c(1, 0), legend.background= element_rect(fill = "transparent", colour = "transparent"), plot.margin=unit(c(.5,.5,.5,.5),"cm")) + 
               labs(color='Transcriptomic type') +
               facet_grid(dummy ~ ., switch="y") +
               ylab("") +theme(axis.text.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.ticks.x=element_blank()) +
               scale_colour_manual(name="", values = c(gg_color_hue(7)[2],gg_color_hue(7)[6],gg_color_hue(7)[3],gg_color_hue(7)[1])) +
               guides(color=FALSE) +
               theme(strip.background = element_blank(), strip.placement = "inside", strip.text = element_blank()) 
))

#scale_colour_manual(values = c("purple", "green", "blue", "yellow", "magenta","orange", "cyan", "red", "black"),
#                    guide = guide_legend(override.aes = list(
#                      linetype = c(rep("blank", 7), "solid", "dashed"),
#                      shape = c(rep(16, 7), NA, NA))))
