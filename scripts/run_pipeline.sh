#!/bin/bash

set -euo pipefail

export NSLOTS=80 #number of cores to use this should be set AUTOMATICALLY ON QUEUE SYSTEMS

export DATA="/single_cell/data"

export SRCPATH="/single_cell/src"
export PYTHONPATH=$SRCPATH
export PATH="${SRCPATH}:/single_cell/scripts:$PATH"

WORKDIR="/single_cell"
cd $WORKDIR

# .feather file: raw countdata stored in the feather file format for rapid loading, rownames are in the first column labeled "rownames"
export FEATHER=`ls $DATA/*.feather`

# .tsv file    : metadata about the cells from the countdata, rownames(coldata) should match colnames(countdata)
export TSV=`ls $DATA/*.tsv`

echo FEATHER=$FEATHER TSV=$TSV

# zinb: an output directory
if [ ! -e "./zinb/zinb_W.csv" ]; then
    mkdir -p zinb
    run_zinb.R $FEATHER $TSV zinb
else
    echo "zinb output file found, skipping"
fi

# Specify cells (usually disease/treated) to exclude from PCA/clustering and then project onto normal/control
# clusters after the fact
mkdir -p inputs
awk -F "\t" 'NR > 1 && /psoriasis/ {print $1}' $TSV > inputs/cellNames_exclude.txt


fi_expr="$FEATHER" #$1 ## csv file with raw counts (genes are rows, columns are cells )
fi_distData="./zinb/zinb_W.csv"               ## csv file with coordinates used used to calcuated cel-cell distances according to MAGIC algorithm (rows are cells , columns are coordinates (e.g PC coords for rows of zinb-wave matrix)  )
fo_imputed="magic_counts_t{}.feather"        ## name for output files {} is filled with corresponding value of MAGIC t parameter
fo_details="runData/impute.runData"          # information on run including runtimes
minFracExpr=0.001                             ## fraction of cells for which raw UMI > 0 needed to preform imputation
tList=10                                      ## comma separated list of t values for which imputation is performed
k=30
ka=10

## use --geneRows to indicate that fi_expr has genes in rows (otherwise it's expected that genes are in columns)
## use --libNormalize to indicate perform MAGIC library nomalization of RAW counts prior to starting imputation (advide always passing this flag)
mkdir -p runData
run_MAGIC.py --fi_expr "$fi_expr" --fi_distData "$fi_distData" --fo_imputed "$fo_imputed"  --fo_details "$fo_details" --geneRows --minFracExpr $minFracExpr --libNormalize --tList $tList -k "$k" --ka "$ka"

f_i="magic_counts_t10.feather"
standLibSize=10000
fo="$f_i"
fo_origLibSizes="origLibSize_magic_counts_t10.csv"
libNorm.py --f_i  "$f_i" --standLibSize "$standLibSize" --fo "$fo" --fo_origLibSizes "$fo_origLibSizes"

fi_expr="$fo"
nPC=50            # get this many PC vectors
logT_pseudo=1     ## a value >=0 indicates log2(x  + pseudo) transform
fi_cellNames=""   ## indicates use all cells
fi_geneNames="./zinb/genes_use.txt"  ## for PCA represent cells by expression of only these genes
fi_project="./inputs/cellNames_exclude.txt"  ## Exclude cell names in this file from PCA fit. Then Project onto PCs
fo="./pca_imputed/magic_counts_t10.csv"

# --randomized flag may speed up
mkdir -p "./pca_imputed"
run_PCA.py --fi_expr "$fi_expr" --nPC "$nPC" --logT_pseudo "$logT_pseudo" --fi_cellNames "$fi_cellNames" --fi_geneNames "$fi_geneNames" --fi_project "$fi_project" --fo "$fo"


##(5) perform KASP clustering on the top 20 PCs  ##############################################
nfeat=20  ## use the 1st 20 columns of f_in (this is 1st 20 PCs)
nClust=10
f_in="./pca_imputed/magic_counts_t10-PCcomps.csv"
fo="./kasp_imputed/nFeat${nfeat}_nClust${nClust}.csv"
fo_runDat="./runData/kasp.runData"
predictOnly="./inputs/cellNames_exclude.txt"
alpha=10         ## the data reduction factor (see KASP reference (../README.txt)). For larger data sets 10 is gives good runtime .
kmeans_frac=0.5  ## Fraction of the data used on which Kmeans round 1 is run (0.5 good for larger data sets). If not provided then Kmeans round 2 is skipped (see ../README.txt)
kmeans_nJobs=-1
ka=10            ## MAGIC parameter for adaptive cell-cell similarity
k=30             ### MAGIC parameter for adaptive cell-cell similarity
N_nearest=10     ## Nearest neighbor method used to assign cells in predictOnly File to clusters (after fitting the clustering). This parameter sets number of nearest neighbors to consider for cluster assignment

mkdir -p kasp_imputed
run_specCluster.py --f_in "$f_in" --fo "$fo" --fo_runDat "$fo_runDat" --predictOnly "$predictOnly" --nfeat "$nfeat" --nClust "$nClust" --alpha "$alpha" --kmeans_frac "$kmeans_frac" --kmeans_nJobs "$kmeans_nJobs" --ka "$ka" --k "$k" --N_nearest "$N_nearest"

# run_clustVis.R orders clusters by KRT5 expression, with the highest PMEL at second to last and highest HLA-DRA cluster second to last
# generates PCA plots colored by cluster, sample, tissue and common markers

CLDIR="cl"
mkdir -p CLDIR
run_clustVis.R kasp_imputed/nFeat20_nClust10.csv $TSV magic_counts_t10.feather $CLDIR

#Use coldata with sorted cluster assignments for clusterDE
run_clustDE_feather.R --no-pso $FEATHER $CLDIR/coldata_clust.csv $CLDIR
