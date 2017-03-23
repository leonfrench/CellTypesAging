########################
#Venn Diagram
#install.packages('venneuler')
library(venneuler)

#input is a tmod result table and the number of top groups - don't send the whole result
getEulerDiagram <- function(goResult) {
  combined <- NULL
  for(group in goResult$ID) {
    combined <- rbind(combined, data.frame(elements= unlist(geneSetsGO$MODULES2GENES[group]), sets=Term(group)))
  }
  
  v <- venneuler(combined)
  plot(v)
  
  (vForGGplot <- tbl_df(data.frame(diameter=v$diameters, v$centers, color=v$colors, goName=v$labels, stringsAsFactors = F)))
  
  #credited to user5061 and baptiste via stackoverflow.com
  circularise <- function(d, n=360){
    angle <- seq(-pi, pi, length = n)
    make_circle <- function(x,y,r,goName){
      data.frame(x=x+r*cos(angle), y=y+r*sin(angle), goName)
    }
    lmat <- mapply(make_circle, goName = d[,"goName"], 
                   x = d[,"x"], y=d[,"y"], r=d[,"diameter"]/2, SIMPLIFY = FALSE)
    do.call(rbind, lmat)
  }
  
  circles <- circularise(vForGGplot)
  vForGGplot <- as.data.frame(vForGGplot)
  diagram <- ggplot(vForGGplot) + geom_blank(aes(x, y)) + 
    geom_polygon(aes(x,y, group=goName, fill=goName), data=circles, alpha=0.4) +
    geom_label_repel(aes(x, y, fill=goName, label = goName), color="black", segment.alpha=0) +
    coord_fixed() + theme_void() + theme(legend.position="none")
  diagram
}


