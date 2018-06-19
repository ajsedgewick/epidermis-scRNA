#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 21:45:57 2018

@author: afinneg2
"""
from __future__ import print_function
import argparse
import numpy as np
import pandas as pd
import sys
import os
import time
from collections import OrderedDict
from sklearn import  neighbors
## Replace /Users/afinneg2/projects/keratinocyteRegulators/processedData/scripts/share with /path/to/this/directory
sys.path.insert(0 , "/src")  ## put cluster_spectral.py on PYTHONPATH
import cluster_spectral as sClust
from utils import readListFromFile

#### Notes: #########################3
### TO DO: Add notes here with references especially indicating the type of spectral clustering used 

#####

##### Parse Command line Arguments ################################################# 
parser = argparse.ArgumentParser(description='Run kasp spectral clustering with adapitve gaussian similarity')
parser.add_argument("--f_in" , required = True , 
                        help = "csv. Samples (e.g cells) are rows columns are features (e.g PCcoords)")
parser.add_argument("--fo" , required =True )
parser.add_argument("--fo_runDat" , default = "" ,help= "name for output file with metaData for run (e.g. run times)"  )
parser.add_argument("--predictOnly" ,default = "" , help = "File name listing samples of f_in that are excluded from fitting the clustering. \
After clustering samples are assigned to cluster of nearest centroid.")
parser.add_argument("--nfeat" , type = int, help = "use the 1st nfeat columnsnof f_in as features")
parser.add_argument("--nClust" , type = int, help = "number clusters")
parser.add_argument("--alpha" ,default = 10.0 ,  type = float, 
  help = "Reduction factor for samples used in spectral clustering. Nsamples_spectral = int(Nsamples/ alpha)")
parser.add_argument("--kmeans_frac" ,  type= float , 
   help = "Fraction of samples to use for 1st round of kmeans. If not provided all samples used and 1 round kmeans")
parser.add_argument("--kmeans1_nInit", default = 20 , type = int, 
 help = "parameter passed to sklearn.cluster  KMeans object" )
parser.add_argument("--kmeans1_maxIter", default = 300, type = int ,
help = "parameter passed to sklearn.cluster  KMeans object" )
parser.add_argument("--kmeans_nJobs" , default = 1, type = int,
help = "parameter passed to sklearn.cluster  KMeans object")
parser.add_argument("--kmeans2_maxIter" , default = 300 , type = int,
help = "parameter passed to sklearn.cluster  KMeans object")
parser.add_argument("--ka" , default =10 , type = int, help = "ka parameter for MAGIC affinity")
parser.add_argument("--k" , default =30 , type = int, help = "k parameter for MAGIC affinity")
parser.add_argument("--N_nearest", default = 10 , type = int ,
                    help = "number of nearest neighbors to use in knn classification of predictOnly samples")
parser.add_argument("--seed" , default =54321 , type = int, help = "specify random seed for reproducibility")
args = parser.parse_args() 


###### Load Data ############################################################
print("Loading Data")
data = pd.read_csv(args.f_in, sep = "," , index_col = 0 )
if args.nfeat is not None:
    data = data.iloc[: , 0:args.nfeat].copy()
if args.predictOnly:
    predictOnly = readListFromFile(args.predictOnly)
    predictMask = data.index.isin(predictOnly )
    print("Using {:d} of {:d} samples for prediction only".format(np.count_nonzero(predictMask) ,
                                                                  len(predictMask)))
    data_fit =data.loc[~predictMask, :  ].copy()
    data_predict = data.loc[predictMask, :].copy()

print("Setting numpy random seed for thie session to: {:d}".format(args.seed))
np.random.seed(args.seed)
    
### run the kasp clustering ####################
print("Running kasp with adaptive gaussian similarity")
sClust_out= sClust.kasp(X = data_fit.values  , nClust = args.nClust,
                                 alpha = args.alpha,
                                 kmeans_frac= args.kmeans_frac, 
                                 kmeans1_kwargs = {"n_init" : args.kmeans1_nInit , "max_iter" : args.kmeans1_maxIter, 'n_jobs' : args.kmeans_nJobs},  
                                 kmeans2_maxIter = args.kmeans2_maxIter ,  
                                 simFunc =   sClust.compute_magicAffinity,
                                 simFunc_kwargs =  {'k':args.k, 'epsilon':1, 'distance_metric' :'euclidean', 'ka' : args.ka})
labels_X, evals, evecs, labels_Xkm, Xkm, _ , runTimes = sClust_out  ## finalLabels, evals , evecs or L_rw , labels assiged to keamns centers , kmeans centers , _ , runtimes
    
#### Construct output files, classifying predictOnly data if necessary
Xkm_data = pd.DataFrame( labels_Xkm, columns= ["clust_ID" ] ).join(
           [ pd.DataFrame( Xkm,  columns =  data.columns.values ) ,
            pd.DataFrame(  evecs , columns = ["evec{}".format(i) for i in range(0,evecs.shape[1])  ] )
           ])

if args.predictOnly:
    ### Assign class to predictOnly simples using knn classifier
    print("predicting with KNeighborsClassifier")
    ti= time.time()
    knn_c = neighbors.KNeighborsClassifier(n_neighbors= args.N_nearest)
    knn_c.fit(X = data_fit.values,  y = labels_X)
    data_predict_labels = knn_c.predict(data_predict.values)
    tf =  time.time()
    runTimes["knn_classifier (s)"] = tf - ti
    ### build output dataframe
    data.loc[~predictMask , "clust_ID"]  =  labels_X
    data.loc[predictMask , "clust_ID"] =  data_predict_labels
else:
    data["clust_ID"] =  labels_X
    
runDetails = OrderedDict([ (key , value) for key ,value in vars(args ).items() ])
runDetails.update(runTimes )

##### Write to file ##################
fo_data = args.fo
data  = data.astype( dtype = {'clust_ID' : np.int32}  )
data.to_csv(fo_data  )

fo_Xkm_data = os.path.splitext(args.fo )[0] + "-Xkms" + os.path.splitext(args.fo)[1]
Xkm_data.to_csv(fo_Xkm_data,  float_format = "%.6f" , index = False  )

fo_evals = os.path.splitext(args.fo )[0] + "-eigVals" + os.path.splitext(args.fo)[1]
np.savetxt(fo_evals  , evals  )

runDat_str ="\t".join([ key +"="+ str(value) for key, value in  runDetails.items() ])
f = open( args.fo_runDat , 'a')
f.write(runDat_str  + "\n" )
f.close()






