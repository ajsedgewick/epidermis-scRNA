#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 2){
    print("Loads scRNA expression data from cellranger mtx directory, QC filter cells, \n does minimal gene filtering and saves countdata and coldata.")
    print("output_directory must contain a tab-delimited file named sample_key.txt with the columns index, library_id and tissue.")
    print("Usage: run_loadQC.R matrix_directory output_directory")
    quit("no", status=1)
}

library(Seurat)
library(feather)
library(tibble)
library(Matrix)

minCells <- 10            # for Seurat object construction
matdir <- args[1] #"/home/jovyan/work/de_aggr/outs/filtered_gene_bc_matrices_mex/GRCh38/"
outdir <- args[2] #"/home/jovyan/work/de_aggr/" OUTDIR MUST HAVE 'sample_key.txt' IN IT

kcyte <- Read10X(data.dir = matdir)
dim(kcyte)
kcyte <- CreateSeuratObject(raw.data = kcyte, min.cells = minCells, project="sc_project", names.delim= "\\-")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = kcyte@data), value = TRUE)
percent.mito <- Matrix::colSums(kcyte@data[mito.genes, ]) / Matrix::colSums(kcyte@data)
kcyte <- AddMetaData(object = kcyte, metadata = percent.mito, col.name = "percent.mito")

bc.sample <- as.numeric(matrix(unlist(strsplit(colnames(kcyte@data), "-")), ncol=2, byrow=T)[,2])
tiss.key <- read.table(paste(outdir, 'sample_key.txt', sep='/'), header=T, sep="\t")
cell.tiss <- tiss.key[bc.sample, 'tissue']
names(cell.tiss) <- kcyte@cell.names

cell.samp <- tiss.key[bc.sample, 'library_id']
names(cell.samp) <- kcyte@cell.names

kcyte <- AddMetaData(kcyte, cell.samp, 'sample')
kcyte <- AddMetaData(kcyte, cell.tiss, 'tissue')

ncell <- c(500, 1000, 3000, 5000, 10000)
mrate <- c(.44, .8, 2.3, 3.9, 7.6)
mrate.fit <- lm(mrate~ncell)
samp.ncell <- data.frame(table(kcyte@meta.data$sample))
mrate.pred <- predict(mrate.fit, data.frame(ncell=samp.ncell$Freq))
names(mrate.pred) <- samp.ncell$Var1

head(kcyte@meta.data)
table(kcyte@meta.data$sample)
table(kcyte@meta.data$tissue)

nGene.samp.filt <- lapply(levels(kcyte@meta.data$sample), 
                            function(x) {
                                cur.cells <- rownames(subset(kcyte@meta.data, sample==x))
                                cur.ng <- kcyte@meta.data[cur.cells, "nGene"]
                                qu <- quantile(cur.ng, 
                                                 probs=c(.005, 1-.02*mrate.pred[x]))
                                filt <-  (cur.ng > qu[1]) & (cur.ng < qu[2])
                                names(filt) <- cur.cells
                                return(filt)
                            })

lapply(nGene.samp.filt, function(x) sum(x) / length(x))
lapply(nGene.samp.filt, length)
nGene.samp.filt <- unlist(nGene.samp.filt)       
nGene.samp.filt <- nGene.samp.filt[rownames(kcyte@meta.data)]

pct.mito.thresh <- quantile(kcyte@meta.data$percent.mito, probs=c(.95))
pct.mito.filt <- kcyte@meta.data$percent.mito < pct.mito.thresh

pct.mito.thresh
sum(pct.mito.filt & nGene.samp.filt)

table(kcyte@meta.data$sample[pct.mito.filt==0])
table(kcyte@meta.data$sample[nGene.samp.filt==0])

cells.use <- rownames(kcyte@meta.data)[pct.mito.filt & nGene.samp.filt]

countdata <- data.frame(as.matrix(kcyte@raw.data[,cells.use]), check.names=F)
countdata <- rownames_to_column(countdata)

write_feather(countdata, paste(outdir, "counts_filt_cells.feather", sep="/"))

write.table(kcyte@meta.data[cells.use,], file=paste(outdir, "coldata.tsv", sep="/"), col.names=NA, row.names=T, sep="\t", quote=F)

#save(kcyte, file=paste(outdir, "kcyte_raw.RData", sep="/"))
#countdata <- as.matrix(kcyte@raw.data)
#write.table(countdata, sep=",", row.names=T, quote=F, col.names=NA, file=paste(outdir, "kcyte_raw_counts.csv", sep="/"))
#coldata <- kcyte@meta.data
#save(coldata, file="kcyte_coldata.RData")

