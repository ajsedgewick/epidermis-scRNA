#!/usr/bin/env Rscript

library(openxlsx)
library(ggplot2)
library(feather)
library(edgeR)
library(tibble)

#source("/home/jovyan/work/code/DE_functions.R")
#source("/home/jovyan/work/analysis_pipeline/plot_functions.R")

args = commandArgs(trailingOnly=TRUE)

kasp.pc <- read.table(args[1], sep=",", row.names=1, header=T) #load(args[1]) #should have kasp results

coldata <- read.table(args[2], sep="\t", row.names=1, header=T) #load(args[3]) #should have coldata

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

findClusterNumbers <- function(kasp.pc, exp.data){
    marker.means <- aggregate(t(exp.data[c("KRT5", "KRT14", "KRT10", "KRT1", "HLA.DRA", "PMEL", "FLG"), rownames(kasp.pc)]),
                          by=list(kasp.pc[, "clust_ID"]), FUN = mean)
    print(marker.means)
    mel.inds <- order(marker.means$PMEL, decreasing=T)[1:2] #which.max(marker.means$PMEL)
    non.krt.cl <- c(mel.inds, which.max(marker.means[,"HLA.DRA"]))
    non.krt.cl <- unique(non.krt.cl) # probably bad news if there are dups here

    # when FLG/LOR cluster is separated out, set as cluster 8
    if(dim(marker.means)[1] == 11){
        non.krt.cl <- c(which.max(marker.means[,"FLG"]), non.krt.cl)
    }

    # use highest krt5/14 for clust 1/2 and KRT1/10 for clust 3-7
    krt.cl <- order(marker.means$KRT5 + marker.means$KRT14, decreasing = T)[1:2]
    krt.cl <- unique(c(krt.cl, order(marker.means$KRT10 + marker.means$KRT1, decreasing = F)))
    new.levels <- marker.means$Group.1[c(krt.cl[!(krt.cl %in% non.krt.cl)], non.krt.cl)]
    new.cluster <- factor(kasp.pc$clust_ID, levels=new.levels)
    levels(new.cluster) <- 1:nlevels(new.cluster)
    return(new.cluster)
}

new.clusts <- findClusterNumbers(kasp.pc, countdata)
coldata[rownames(kasp.pc), "cluster"] <- new.clusts

write.table(coldata, file=paste(outdir, "coldata_clust.csv", sep="/"), quote=F,
           col.names=NA, row.names=T, sep=",")
write.xlsx(list(clustCount = table(coldata$cluster), sampClust=table(coldata$cluster, coldata$sample), tissClust=table(coldata$cluster, coldata$tissue)),
           file=paste(outdir, "cell_counts.xlsx", sep="/"),)

markers <- c("KRT5", "KRT10", "HLA.DRA", "PMEL", "FLG", "LOR")
cluster.colors <- c('#008B00FF', '#7CAE00FF', '#FFA500FF', '#FFD700FF', '#00BFC4FF',
                    '#308ED2FF', '#7771B5FF', '#912CEEFF',  '#8B5A2BFF','#FF5347FF', '#CD853FFF')

plot.df <- cbind(kasp.pc[,1:2], coldata[rownames(kasp.pc),c("sample", "tissue", "cluster")], t(countdata[markers, rownames(kasp.pc)])) 
kasp.pc[rownames(coldata), "cluster"] <- coldata$cluster

#p <- ggplot(plot.df, aes_string("PC1", "PC2", col="cluster")) + geom_point(alpha=.5) + theme_bw() +
#guides(colour = guide_legend(override.aes = list(alpha = 1, size=4))) # + scale_color_manual(values = cluster.colors)
#ggsave(p, file=paste(outdir, paste0("pca_cluster.png"), sep="/"))

for(x in c("sample", "cluster", "tissue")){
  plot.df[,x] <- as.factor(plot.df[,x])
  p <- ggplot(plot.df, aes_string("PC1", "PC2", col=x)) + geom_point(alpha=.5) + theme_bw() +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size=4)))
  if(nlevels(plot.df[,x]) <= 11) {
    p <- p + scale_color_manual(values = cluster.colors)
  }
  ggsave(p, file=paste(outdir, paste0("pca_",x,".png"), sep="/"))
}

for(x in markers){
  p <- ggplot(plot.df, aes_string("PC1", "PC2", col=x)) + geom_point(alpha=.5) + theme_bw() 
  ggsave(p, file=paste(outdir, paste0("pca_",x,".png"), sep="/"))
}
