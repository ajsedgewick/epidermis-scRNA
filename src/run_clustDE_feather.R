#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args) !=3){
    print("Usage: run_clustDE_feather.R countdata.feather coldata.tsv output_directory")
    quit("no", status=1)
}

library(feather)
library(tibble)

source("/single_cell/src/DE_functions.R")
#source("/netapp/home/labjbc/src/DE_functions.R")
#source("/home/jovyan/work/analysis_pipeline/plot_functions.R")

# MAKE SURE TO KNOW IF COUNTS ARE (LOG) CPM
countdata <- read_feather(args[1]) # should have countdata
countdata <- column_to_rownames(data.frame(countdata, check.names=F))

coldata <- read.table(args[2], sep=",", header=T, row.names=1, check.names=F) #should have coldata 

outdir <- args[3]
# Some code to read cluster assignments from a separate file
#
#if(length(args)==3) {
#  outdir <- args[3]
#} else if(length(args)>=4){
#  clusts <- read.table(args[3], sep=",", row.names=1, header=F)
#  if(min(clusts)==0){
#    clusts <- clusts + 1
#  }
#
#  countdata <- countdata[,rownames(clusts)]
#  coldata <- coldata[rownames(clusts),]
#  coldata$cluster <- as.factor(clusts[,1])
#  outdir <- args[4]
#}

normal.cells <- rownames(subset(coldata, tissue == "trunk"))

coldata$cluster <- factor(coldata$cluster)

coldata <- coldata[normal.cells,]
coldata$sample <- droplevels(coldata$sample)
coldata$tissue <- droplevels(coldata$tissue)

countdata <- countdata[,normal.cells]

#DE of each cluster versus the rest
clust.de <- clustDE(coldata, countdata, logcpm=F,
                    outfile=paste(outdir, "clusterDE.xlsx", sep="/"))



#other DE methods from 1st paper
#cl.FC <- make.FC.mat(clust.de)
#cluster.plot.fc(cl.FC, outfile="clusterFCplot.pdf")

#For each cluster perform DE on cells from each cluster vs the rest
#tiss.de <- runTissClustDE(coldata, countdata,
#           outfile=paste(outdir, "clusterTissueDE.xlsx", sep="/"))


#cl.tiss.FC <- make.FC.mat(tiss.de)

#foreskin.FC <- format.FC(cl.tiss.FC, seq(1, ncol(cl.tiss.FC), 3))
#tissue.plot.fc(foreskin.FC,
#               outfile=paste(outdir, "foreskinFCHeatmap.pdf", sep="/"))

#scalp.FC <- format.FC(cl.tiss.FC, seq(2, ncol(cl.tiss.FC), 3))
#tissue.plot.fc(scalp.FC,
#               outfile=paste(outdir, "scalpFCHeatmap.pdf", sep="/"))

#trunk.FC <- format.FC(cl.tiss.FC, seq(3, ncol(cl.tiss.FC), 3))
#tissue.plot.fc(trunk.FC,
#               outfile=paste(outdir, "trunkFCHeatmap.pdf", sep="/"))


