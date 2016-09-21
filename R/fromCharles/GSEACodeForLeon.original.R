######## GSEA Code for Leon ####
#> Sys.Date()
#[1] "2016-08-24"

setwd('/Users/matianzhou/Documents/Pitt Spring 2016/Etienne/Singe Cell RNAseq/FromLeon160720')
load('AgingFullResults_UniqueGenes.RData')
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

########## prepare the input data 
unique.aging.genes <- rownames(es.dat.unique.genes)[which(sign(es.dat.unique.genes[,1])*sign(es.dat.unique.genes[,2])==1)]
unique.up.genes <- unique.aging.genes[which(sign(es.dat.unique.genes[unique.aging.genes,1])==1)]
unique.down.genes <- unique.aging.genes[which(sign(es.dat.unique.genes[unique.aging.genes,1])== -1)]
dat <- data.frame(geneNames = c(unique.up.genes,unique.down.genes),
                  direction=c(rep(1,length(unique.up.genes)),rep(-1,length(unique.down.genes))),
                  agingP = c(-log10(aw.out[unique.up.genes,4]),-log10(aw.out[unique.down.genes,4])), 
                  cellMarker=rep("BG",length(unique.up.genes)+length(unique.down.genes)),stringsAsFactors = F)
rownames(dat) <- dat[,1]
dat$agingPValuesWithDirection <- c(dat[dat$direction >0, "agingP"], -dat[dat$direction <0, "agingP"]) 

library(snowfall)
B <- 1000 ## number of permutations
bobs.vector <- rep(NA,length(enriched_genes))
names(bobs.vector) <- names(enriched_genes)
bperm.mat <- matrix(,nrow=B,ncol=length(enriched_genes))
colnames(bperm.mat) <- names(enriched_genes)
p.vector <- rep(NA,length(enriched_genes))
names(p.vector) <- names(enriched_genes)

for (i in 1:length(enriched_genes)) {
  print(paste('Cell line ',i,sep=""))
  start.time <- proc.time()
  ## observed part
  ingenes <- intersect(enriched_genes[[i]],rownames(dat))
  outgenes <- setdiff(rownames(dat),ingenes)
  path = dat[ingenes,"agingPValuesWithDirection"]
  nonpath = dat[outgenes,"agingPValuesWithDirection"]
  g1 <- length(path)
  g2 <- length(nonpath)
  total=data.frame(c(path,nonpath))
  total[,2]=c(rep(1,g1),rep(0,g2 ))
  total[,3]=rank(-total[,1],ties.method = "first")
  inpath=total[1:g1,]
  outpath=total[(g1+1):(g1+g2),]
  J=as.numeric(rep(0,g1+g2))
  
  for (j in 1:(length(J))) {
    
    d=inpath[inpath[,3]<=j,]
    b1=sum(abs(d[,1]))/sum(abs(inpath[,1]))
    t2=outpath[outpath[,3]<=j,]
    b2=nrow(t2)/g2
    J[j]=abs(b1-b2)
  }
  
  bobs.vector[i] = bobsw= max(J)  # observed weighted KS statistics
  
  ## permuted part 
  
  chunks = 10
  each <- B/chunks
  
  sfInit(parallel=T,cpus=chunks,type="SOCK")
  perm_functions <- function(xxx) {
    bdistw=as.numeric(rep(0,each))  
    for (b in 1:each)   {
      index <- sample(1:nrow(total),size=nrow(total))
      total_perm=total[index,]
      total_perm[,2]=c(rep(1,g1),rep(0,g2 ))
      total_perm[,3]=rank(-total_perm[,1],ties.method = "first")
      inpath=total_perm[1:g1,]
      outpath=total_perm[(g1+1):(g1+g2),]
      J=as.numeric(rep(0,g1+g2))
      
      for (j in 1:(length(J))) {
        
        d=inpath[inpath[,3]<=j,]
        b1=sum(abs(d[,1]))/sum(abs(inpath[,1]))
        t2=outpath[outpath[,3]<=j,]
        b2=nrow(t2)/g2        
        J[j]=abs(b1-b2)
       }      
      bdistw[b]=max(J)   
    }
    return(bdistw)
  }
  
  sfExport("perm_functions") 
  sfExport("total")
  sfExport("g1")
  sfExport("g2")
  sfExport("each")
  result<-sfLapply(1:chunks, perm_functions) 
  sfStop()
  
  bperm.mat[,i] = bjoin = c(unlist(result))
#p.vector[i]=p_ksw=(length(bjoin[bjoin>=bobsw])+1)/(length(bjoin)+1) # p-value for weighted KS test  
  end.time <- proc.time()
  print(end.time-start.time)
}

bperm <- c(bperm.mat)
p.share <- rep(NA,length(bobs.vector))
for (i in 1:length(bobs.vector)) {
  b.i <- bobs.vector[i]
  p.share[i] <- (length(bperm[bperm>=b.i])+1)/(length(bperm)+1)
}
names(p.share) <- names(p.vector) <- names(enriched_genes)

