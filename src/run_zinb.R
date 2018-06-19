#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
if(length(args) !=3){
    print("Usage: run_zinb.R countdata.feather coldata.tsv output_directory")
    quit("no", status=1)
}

library(Seurat)
library(Matrix)
library(SummarizedExperiment)
library(zinbwave)
library(BiocParallel)
library(edgeR)
library(feather)
library(tibble)


countdata <- read_feather(args[1]) #as.matrix(kcyte@raw.data[,kcyte@cell.names])
countdata <- as.matrix(column_to_rownames(data.frame(countdata, check.names=F)))

coldata <- read.table(args[2], header=T, row.names=1, check.names=F, sep="\t")

outdir <- args[3]

#cpm <- cpm(countdata, log = F, prior.count = 0)
#is_quality <- rowSums(cpm  >= 500 ) >= .001*ncol(cpm) # equivalent of 5 cp10k

is_quality <- rowSums(countdata >= 5) >= 100 #5 umi in 100 cells
sum(is_quality)

genes.use <- rownames(countdata)[is_quality]
countdata <-  countdata[genes.use,]
dim(countdata)

#coldata$lognUMI <- log2(coldata$nUMI + 1)

cur.se <- SummarizedExperiment(countdata,
                               colData=coldata)

cores <- as.numeric(Sys.getenv("NSLOTS"))
			       
zbc <- zinbFit(cur.se, K=20, X="~nUMI + sample + percent.mito", epsilon=1000, BPPARAM=MulticoreParam(cores))
save(zbc, file=paste(outdir, "zinb_fit.RData", sep="/"))

se_norm <- zinbwave(cur.se, fitted_model=zbc, normalizedValues=TRUE,
                    residuals = TRUE, BPPARAM=MulticoreParam(cores))
save(se_norm, file=paste(outdir, "zinb_wavese.RData", sep="/"))
W <- getW(zbc)
colnames(W) <- paste0("W", 1:20)
rownames(W) <- colnames(assays(se_norm)$normalizedValues)
write.table(round(W, 5), file=paste(outdir, "zinb_W.csv", sep="/"),
            sep=",", col.names=NA, row.names=T, quote=F)

write.table(genes.use, file=paste(outdir, "genes_use.txt", sep="/"), col.names=F, row.names=F, quote=F)