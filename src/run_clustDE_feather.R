#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
usage =   "Usage: run_clustDE_feather.R [--no-pso] countdata.feather coldata.tsv output_directory"
if(length(args) <3){
    print(usage)
    quit("no", status=1)
}

no.pso = F
if(length(args) ==4){
    if(args[1] == "--no-pso"){
        no.pso = T
        args = args[-1]
    } else {
        print(usage)
        quit("no", status=1)
    }
}

library(feather)
library(tibble)


if(file.exists(paste(Sys.getenv("SRCPATH"), "DE_functions.R"))) {
  source(paste(Sys.getenv("SRCPATH"), "DE_functions.R"))
} else if(file.exists("../src/DE_functions.R")) {
  source("../src/DE_functions.R")
} else {
  print("Can't find DE_functions.R")
  quit("no", status=1)
}

#source("/single_cell/src/DE_functions.R")
#source("/netapp/home/labjbc/src/DE_functions.R")
#source("/home/jovyan/work/analysis_pipeline/plot_functions.R")

# MAKE SURE TO KNOW IF COUNTS ARE (LOG) CPM
countdata <- read_feather(args[1]) # should have countdata
countdata <- column_to_rownames(data.frame(countdata, check.names=F))

#filter genes
is_quality <- rowSums(countdata >= 3) >= 100 #3 umi in 100 cells
sum(is_quality)

genes.use <- rownames(countdata)[is_quality]
countdata <-  countdata[genes.use,]


coldata <- read.table(args[2], sep=",", header=T, row.names=1, check.names=F) #should have coldata 

outdir <- args[3]

cells.use <- rownames(coldata)

if(no.pso) {
    cells.use <- rownames(subset(coldata, tissue != "psoriasis"))
}

coldata$cluster <- factor(coldata$cluster)

coldata <- coldata[cells.use,]
coldata$sample <- droplevels(coldata$sample)
coldata$tissue <- droplevels(coldata$tissue)

countdata <- countdata[,cells.use]

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
