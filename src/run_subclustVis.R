#!/usr/bin/env Rscript

library(openxlsx)
library(ggplot2)
library(feather)
library(edgeR)
library(tibble)

getwd()

if(file.exists(paste(Sys.getenv("SRCPATH"), "hugePalette.R"))) {
  source(paste(Sys.getenv("SRCPATH"), "hugePalette.R"))
} else if(file.exists("../src/hugePalette.R")) {
  source("../src/hugePalette.R")
} else {
  print("Cant find hugePalette.R")
  quit("no", status=1)
}


args = commandArgs(trailingOnly=TRUE)

kasp.pc <- read.table(args[1], sep=",", row.names=1, header=T) #load(args[1]) #should have kasp results

coldata <- read.table(args[2], sep="\t", row.names=1, header=T) #load(args[3]) #should have coldata

coldata <- coldata[rownames(kasp.pc),]

countdata <- read_feather(args[3])
countdata <- column_to_rownames(data.frame(countdata, check.names=F))
countdata <- t(countdata)
#countdata <- edgeR::cpm(countdata, log = TRUE, prior.count = 1)

rownames(countdata)[rownames(countdata)=="HLA-DRA"] <- "HLA.DRA"

outdir <- args[4]

#load(args[2])
#if(length(args==3){

#} else if(length(args==4){
#  magic.data <- read.table(args[3], sep=",", row.names=1, header=T) #load(args[2]) #magic data
#  outdir <- args[4] #directory for plots, new coldata
#}

c("KRT5", "KRT10", "HLA.DRA", "PMEL") %in% rownames(countdata)

coldata[rownames(kasp.pc), "cluster"] <- paste0("cl", kasp.pc$clust_ID + 1)

write.table(coldata, file=paste(outdir, "coldata_clust.csv", sep="/"), quote=F,
           col.names=NA, row.names=T, sep=",")
write.xlsx(list(clustCount = table(coldata$cluster), sampClust=table(coldata$cluster, coldata$sample), tissClust=table(coldata$cluster, coldata$tissue)),
           file=paste(outdir, "cell_counts.xlsx", sep="/"),)

markers <- c("KRT5", "KRT10", "HLA.DRA", "PMEL", "FLG", "LOR")

plot.df <- cbind(kasp.pc[,1:2], coldata[rownames(kasp.pc),c("sample", "tissue", "cluster")], t(countdata[markers, rownames(kasp.pc)])) 
kasp.pc[rownames(coldata), "cluster"] <- coldata$cluster

#p <- ggplot(plot.df, aes_string("PC1", "PC2", col="cluster")) + geom_point(alpha=.5) + theme_bw() +
#guides(colour = guide_legend(override.aes = list(alpha = 1, size=4))) # + scale_color_manual(values = cluster.colors)
#ggsave(p, file=paste(outdir, paste0("pca_cluster.png"), sep="/"))



for(x in c("sample", "cluster", "tissue")){

  plot.df[,x] <- as.factor(plot.df[,x])
  if(nlevels(plot.df[,x]) == 1){
      next
  }
  print(paste("factor table: ", x))
  
  print(table(plot.df[,x]))
  plot.colors <- hugePalette[1:nlevels(plot.df[,x])]
  p <- ggplot(plot.df, aes_string("PC1", "PC2", col=x)) + geom_point(alpha=.5) + theme_bw() +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size=4)))
  if(nlevels(plot.df[,x]) <= 11) {
    p <- p + scale_color_manual(values = plot.colors)
  }
  ggsave(p, file=paste(outdir, paste0("pca_",x,".png"), sep="/"))
}

#for(x in markers){
#  p <- ggplot(plot.df, aes_string("PC1", "PC2", col=x)) + geom_point(alpha=.5) + theme_bw() 
#  ggsave(p, file=paste(outdir, paste0("pca_",x,".png"), sep="/"))
#}
