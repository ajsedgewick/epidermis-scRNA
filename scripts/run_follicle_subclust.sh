#!/bin/bash

set -euo pipefail

export NSLOTS=80 #number of cores to use this should be set AUTOMATICALLY ON QUEUE SYSTEMS

export DATA="/single_cell/run_dir/data"

export SRCPATH="/single_cell/src"
export PYTHONPATH=$SRCPATH
export PATH="${SRCPATH}:/single_cell/scripts:$PATH"

WORKDIR="/single_cell/run_dir"
cd $WORKDIR

# .feather file: raw countdata stored in the feather file format for rapid loading, rownames are in the first column labeled "rownames"
export FEATHER=`ls $DATA/*.feather`

# .tsv file    : metadata about the cells from the countdata, rownames(coldata) should match colnames(countdata)
export TSV=`ls $DATA/*.tsv`

echo FEATHER=$FEATHER TSV=$TSV
export PREFIX="follicle"

# Specify cells to use exclusively for PCA/subclustering
# IMPORTANT this assumes run_clustVis successfully changed the immune cluster to cluster index 11 and melanocytes to cluster 9 and 10
mkdir -p inputs
awk -F "," 'NR==1 {for(i=1;i<=NF;i++){if ($i ~ /cluster/) {clcol=i} else if ($i ~ /tissue/) {tiscol=i}}}
            NR>1 && $clcol<9 && $tiscol ~ /scalp/ {print $1}' ./cl_immune/coldata_clust.csv > inputs/${PREFIX}_cell_names.txt


fi_expr="magic_counts_t10.feather"
nPC=50            # get this many PC vectors
logT_pseudo=1     ## a value >=0 indicates log2(x  + pseudo) transform
fi_cellNames="./inputs/${PREFIX}_cell_names.txt"   ## indicates use all cells
fi_geneNames="./zinb/genes_use.txt"  ## for PCA represent cells by expression of only these genes
# fi_project=""  ## Exclude cell names in this file from PCA fit. Then Project onto PCs --fi_project "$fi_project" 
fo="./pca_imputed/${PREFIX}_magic_counts_t10.csv"

# --randomized flag may speed up
mkdir -p "./pca_imputed"
run_PCA.py --fi_expr "$fi_expr" --nPC "$nPC" --logT_pseudo "$logT_pseudo" --fi_cellNames "$fi_cellNames" --fi_geneNames "$fi_geneNames" --fo "$fo" 


##(5) perform KASP clustering on the top 20 PCs  ##############################################
nfeat=20  ## use the 1st 20 columns of f_in (this is 1st 20 PCs)
nClust=15
f_in="./pca_imputed/${PREFIX}_magic_counts_t10-PCcomps.csv"
fo="./kasp_imputed/${PREFIX}_nFeat${nfeat}_nClust${nClust}.csv"
fo_runDat="./runData/${PREFIX}_kasp.runData"
predictOnly="" ## inputs/cellNames_exclude.txt
alpha=2         ## the data reduction factor (see KASP reference (../README.txt)). For larger data sets 10 is gives good runtime .
kmeans_frac=0.7  ## Fraction of the data used on which Kmeans round 1 is run (0.5 good for larger data sets). If not provided then Kmeans round 2 is skipped (see ../README.txt)
kmeans_nJobs=1
ka=10            ## MAGIC parameter for adaptive cell-cell similarity
k=20             ### MAGIC parameter for adaptive cell-cell similarity
N_nearest=10     ## Nearest neighbor method used to assign cells in predictOnly File to clusters (after fitting the clustering). This parameter sets number of nearest neighbors to consider for cluster assignment

mkdir -p kasp_imputed
run_specCluster.py --f_in "$f_in" --fo "$fo" --fo_runDat "$fo_runDat" --predictOnly "$predictOnly" --nfeat "$nfeat" --nClust "$nClust" --alpha "$alpha" --kmeans_nJobs "$kmeans_nJobs" --ka "$ka" --k "$k" --N_nearest "$N_nearest" --kmeans_frac "$kmeans_frac"

# run_clustVis.R orders clusters by KRT5 expression, with the highest PMEL at second to last and highest HLA-DRA cluster second to last
# generates PCA plots colored by cluster, sample, tissue and common markers

CLDIR="cl_${PREFIX}"
mkdir -p $CLDIR
run_subclustVis.R kasp_imputed/${PREFIX}_nFeat${nfeat}_nClust${nClust}.csv $TSV magic_counts_t10.feather $CLDIR

#Use coldata with sorted cluster assignments for clusterDE
run_clustDE_feather.R $FEATHER $CLDIR/coldata_clust.csv $CLDIR
