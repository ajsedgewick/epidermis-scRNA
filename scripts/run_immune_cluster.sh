#!/bin/bash


export SRCPATH="/single_cell/src"
export PYTHONPATH=$SRCPATH
export PATH="${SRCPATH}:/single_cell/scripts:$PATH"

export WORKINGDIR="/single_cell/run_dir"
cd $WORKINGDIR

export DATA="${WORKINGDIR}/data"

# .feather file: raw countdata stored in the feather file format for rapid loading, rownames are in the first column labeled "rownames"
export FEATHER=`ls $DATA/*.feather`

# .tsv file    : metadata about the cells from the countdata, rownames(coldata) should match colnames(countdata)
export TSV=`ls $DATA/*.tsv`


# pull out PCs of immune cluster 10
awk -F "," 'NR==FNR && NR==1 {for(i=1; i <=NF; i++) head[$i] = i }
            NR==FNR && NR>1 && $head["cluster"] == 10 {cells[$1]=$1}
            NR!=FNR && FNR==1 {print $0}
            NR!=FNR && $1 in cells {print $0}' ./cl/coldata_clust.csv ./pca_imputed/magic_counts_t10-PCcomps.csv > ./pca_imputed/immune_PCcomps.csv

# psoriasis cells excluded from clustering
sort <(cut -d "," -f1 ./pca_imputed/immune_PCcomps.csv) <(cat ./inputs/cellNames_exclude.txt) | uniq -d > ./inputs/immune_cells_exclude.txt


##(5) perform KASP clustering on the top 20 PCs  ##############################################
nfeat=20  ## use the 1st 20 columns of f_in (this is 1st 20 PCs)
nClust=2
f_in="./pca_imputed/immune_PCcomps.csv"
fo="./kasp_imputed/nFeat${nfeat}_nClust${nClust}_immune.csv"
fo_runDat="./runData/kasp_immune.runData"
predictOnly="./inputs/immune_cells_exclude.txt"
alpha=1         ## the data reduction factor (see KASP reference (../README.txt)). For larger data sets 10 is gives good runtime .
#kmeans_frac=0.5  ## Fraction of the data used on which Kmeans round 1 is run (0.5 good for larger data sets). If not provided then Kmeans round 2 is skipped (see ../README.txt)
kmeans_nJobs=-1
ka=3            ## MAGIC parameter for adaptive cell-cell similarity
k=10             ### MAGIC parameter for adaptive cell-cell similarity
N_nearest=10     ## Nearest neighbor method used to assign cells in predictOnly File to clusters (after fitting the clustering). This parameter sets number of nearest neighbors to consider for cluster assignment

run_specCluster.py --f_in "$f_in" --fo "$fo" --fo_runDat "$fo_runDat" --predictOnly "$predictOnly" --nfeat "$nfeat" --nClust "$nClust" --alpha "$alpha" --kmeans_nJobs "$kmeans_nJobs" --ka "$ka" --k "$k" --N_nearest "$N_nearest"

# merge cluster results together
awk -F "," 'BEGIN {OFS=","}
            NR==1 {print $0; next}
            NR==FNR {cells[$1]=$1; $NF = $NF + 10; print $0}
            NR > FNR && FNR > 1 && !($1 in cells) {print $0}' ./kasp_imputed/nFeat20_nClust2_immune.csv ./kasp_imputed/nFeat20_nClust10.csv > ./kasp_imputed/nFeat20_nClust10_immune2_merged.csv

mkdir -p cl_immune
run_clustVis.R ./kasp_imputed/nFeat20_nClust10_immune2_merged.csv $TSV magic_counts_t10.feather cl_immune

#Use coldata with sorted cluster assignments for clusterDE
run_clustDE_feather.R $FEATHER ./cl_immune/coldata_clust.csv cl_immune
