#################################################################################
## Specifications to set up graph analysis
## This should be sourced by scripts that import files or compute metrics
if(!exists("basedir")) {
  stop("You must define the basedir for the project")
} else if (!file.exists(basedir)) {
  stop("basedir: ", basedir, " does not exist.")
}

source("functions/get_subj_info.R")
source("functions/calcGraph_binary.R")
source("functions/import_adj_mats.R")
source("functions/setup_graphs.R")
source("functions/setup_community.R")
source("functions/graph_util_redux.R")
source("functions/run_parse_deltacon.R")
source("functions/wibw_module_degree.R")
source("functions/analyze_nodal_metrics_PCA_func.R")
source("Infomap/run_infomap.R")

nnodes <- 269  #varying this will require you to change the .txt coordinate file you read in for roiMat (should be 4 X nnodes txt file labelled:"x", "y", "z", "roi") and the masterlookup atlas
parcellation <- "power269"
roiFile <- file.path(basedir, "data", "bb264coordinate_appended_shift_nate_culled.txt")
roiMat <- read.table(roiFile, header=FALSE, col.names=c("x", "y", "z", "roi"))
atlasFile <- paste0(parcellation, "_atlas.csv")
atlas <- read.csv(file.path(basedir, "data", atlasFile), header = TRUE, stringsAsFactors = FALSE) #eventually would be nice to have a lockstep name here and re-use parcellationName
use.infomap <- 1

preproc_pipeline <- "aroma" #method for data preprocessing. This corresponds to mni_5mm_aroma data

#conn_method <- "pearson" #Jun2017: aroma preprocessing, pearson correlations
#conn_method <- "cor.shrink" #Jun2017: aroma preprocessing, shrinkage estimator of correlation
#conn_method <- "pcor.shrink_partial" #Jun2017: aroma preprocessing, shrinkage estimator of *partial* correlation
conn_method <- "ridge.net_partial" #Jun2017: aroma preprocessing, shrinkage estimator of *partial* correlation

#conn_method <- "pearson_fisherz" #uses older files (I believe from wavelet)


#conn_method  <- "scotmi"
roi.dist <- 20 #distance in mm to zero connections when setting up graphs (cf. Power 2011)
metricstorun.nodal <- c("eigen.cent","degree", "closeness", "betweenness.node", "page.rank",  "part.coeff", "within.module.deg.zscore", "local.clustering", "gateway.coeff.btw", "gateway.coeff.degree", "between.module.deg.zscore")
metricstorun.global <- c("characteristic_path_length", "clustering_coefficient", "small_worldness", "modularity")
adjmats_base <- file.path(basedir, "adjmats")
#adjmats_base <- file.path("/Users/mnh5174/Data_Analysis/bpd_rest", "adjmats") #should be above, but working from a local copy to avoid messing up existing things

densities_desired <- seq(.01, .2, .01) #used globally for binary graphs

message("Initializing RS graph analysis with the following settings: ")
message("Parcellation: ", parcellation)
message("roiMat: ", roiFile)
message("atlas: ", atlasFile)
message("preprocessing pipeline: ", preproc_pipeline)
message("adj mat calculation conn_method: ", conn_method)
message("adjmats directory: ", adjmats_base)

#estimate ROI-ROI Euclidean distance (use cache if available since the double for loop takes ~4 seconds)
expectFile <- file.path(basedir, "cache", paste0("distmat_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
if (file.exists(expectFile)) {
  load(expectFile)
} else {
  ################################################################################
  ### setup Euclidean distance
  roiDist <- matrix(NA, nrow=nnodes, ncol=nnodes) 
  
  system.time(for (i in 1:nnodes) {
        for (j in 1:nnodes) {
          #populate lower triangle only for speed, then transpose onto upper triangle
          if (i > j) { roiDist[i,j] <- sqrt((roiMat[i,1] - roiMat[j,1])^2 + (roiMat[i,2] - roiMat[j,2])^2 + (roiMat[i,3] - roiMat[j,3])^2) }
        }
      })
  
  #copy lower triangle to upper triangle
  roiDist[upper.tri(roiDist)] <- t(roiDist)[upper.tri(roiDist)]
  diag(roiDist) <- 0
  
  ############################quick QA
  # hist(roiDist)
  # vecDist <- roiDist[lower.tri(roiDist)]
  # sum(vecDist < 20)/length(vecDist)
  
  #roi.dist can be changed at the front end of the script
  rmShort <- roiDist > roi.dist
  #creates binary matrix, in which 0 denotes a short distanced connection that is to be removed. 
  rmShort <- apply(rmShort, c(1,2), as.numeric)
  save(file=expectFile, rmShort, roiDist)
}

## Load necessary libraries
suppressMessages(library(igraph))
suppressMessages(library(RGtk2))
suppressMessages(library(brainGraph))
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(gdata))
suppressMessages(library(beepr))
suppressMessages(library(parallel))
suppressMessages(library(data.table))
suppressMessages(library(psych))
suppressMessages(library(tidyr))
library(reshape2)
suppressMessages(library(broom))
suppressMessages(library(gtools))
suppressMessages(library(Hmisc))
