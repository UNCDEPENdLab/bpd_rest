#################################################################################
## Specifications to set up graph analysis
## This should be sourced by scripts that import files or compute metrics
if(!exists("basedir")) {
  stop("You must define the basedir for the project")
} else if (!file.exists(basedir)) {
  stop("basedir: ", basedir, " does not exist.")
}

######################################################################
###SOURCE CUSTOM FUNCTIONS
######################################################################
source("functions/get_subj_info.R")
source("functions/calcGraph_binary.R")
source("functions/import_adj_mats.R")
source("functions/setup_graphs.R")
source("functions/setup_community.R")
source("functions/graph_util_redux.R")
source("functions/run_parse_deltacon.R")
source("functions/wibw_module_degree.R")
source("functions/analyze_nodal_metrics_PCA_func.R")
source("functions/edge_functions.R")
source("Infomap/run_infomap.R")
source("functions/interrogate_node.R")
source("functions/generate_heatmaps.R")
source("functions/PCA_a_priori.R")
source("functions/PCA_all.R")
source("functions/run_group_comparisons_nodal.R")
source("functions/plot_sig_results_nodal.R")
source("functions/plot_significant_ageeffects.R")
source("functions/plot_significant_groupeffects.R")
#source("functions/plot_significant_ixn.R") ##issue, since this calls on subj_info it must be read in later
source("functions/export_bpd_sig_nodefiles.R")
source("functions/generate_agg.g.R")


######################################################################
###SPECIFY PIPELINE INPUTS
######################################################################
nnodes <- 422  #varying this will require you to change the .txt coordinate file you read in for roiMat (should be 4 X nnodes txt file labelled:"x", "y", "z", "roi") and the masterlookup atlas
parcellation <- "schaefer422"
roiFile <- file.path(basedir,  "data", "schaefer422_roiMat.txt")
roiMat <- read.table(roiFile, header=FALSE, col.names=c("x", "y", "z", "roi"))
atlasFile <- paste0(parcellation, "_atlas.csv")
atlas <- read.csv(file.path(basedir, "data", atlasFile), header = TRUE, stringsAsFactors = FALSE) #eventually would be nice to have a lockstep name here and re-use parcellationName
use.infomap <- 0 #from earlier community detection efforts, has since been decommissioned
use.yeo <-1
PCA.method <- "a_priori" #can be set to "all", "a_priori", or "metrics". This specifies how PCA is performed across densities and metrics
figure.dir <- file.path(basedir, "figures/")
node.file.dir <- file.path(basedir, "Viz_files")

preproc_pipeline <- "aroma" #method for data preprocessing. This corresponds to mni_5mm_aroma data. 
##UPDATE 8/8/17. no 5mm spatial smoothing re: Alakorkko et al 2017


#conn_method <- "pearson" #Jun2017: aroma preprocessing, pearson correlations
#conn_method <- "cor.shrink" #Jun2017: aroma preprocessing, shrinkage estimator of correlation
#conn_method <- "pcor.shrink_partial" #Jun2017: aroma preprocessing, shrinkage estimator of *partial* correlation
#conn_method <- "ridge.net_partial" #Jun2017: aroma preprocessing, shrinkage estimator of *partial* correlation
conn_method <- "dens.clime_partial" #Sept2017 aroma preprocessing, density-based approach for partial correlation estimation (uses Constrained L1-Minimization (CLIME))
#conn_method <- "quic" #coming soon
#conn_method <- "pearson_fisherz" #uses older files (from wavelet 5mm on power 269)
#conn_method  <- "scotmi"  #decommissioned due to not enough time points in our data


roi.dist <- 20 #distance in mm to zero connections when setting up graphs (cf. Power 2011)
metricstorun.nodal <- c("eigen.cent","degree", "closeness", "betweenness.node", "page.rank",  "part.coeff", "within.module.deg.zscore", "local.clustering", "gateway.coeff.btw", "gateway.coeff.degree", "between.module.deg.zscore")
metricstorun.global <- c("characteristic_path_length", "clustering_coefficient", "small_worldness", "modularity")
pcametrics <- c("eigen.cent","degree", "closeness", "betweenness.node", "page.rank",  "part.coeff", "within.module.deg.zscore", "gateway.coeff.btw", "gateway.coeff.degree")

adjmats_base <- file.path(basedir, "adjmats")
#adjmats_base <- file.path("/Users/mnh5174/Data_Analysis/bpd_rest", "adjmats") #should be above, but working from a local copy to avoid messing up existing things

densities_desired <- seq(.01, .25, .01) #used globally for binary graphs

message("Initializing RS graph analysis with the following settings: ")
message("Parcellation: ", parcellation)
message("roiMat: ", roiFile)
message("atlas: ", atlasFile)
message("preprocessing pipeline: ", preproc_pipeline)
message("adj mat calculation conn_method: ", conn_method)
message("adjmats directory: ", adjmats_base)
message("PCA method: ", PCA.method)

######################################################################
#estimate ROI-ROI Euclidean distance (use cache if available since the double for loop takes ~4 seconds)
expectFile <- file.path(basedir, "cache", paste0("distmat_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
if (file.exists(expectFile)) {
  load(expectFile)
} else {
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
  
  ####quick QA
  # hist(roiDist)
  # vecDist <- roiDist[lower.tri(roiDist)]
  # sum(vecDist < 20)/length(vecDist)
  
  #roi.dist can be changed at the front end of the script
  rmShort <- roiDist > roi.dist
  #creates binary matrix, in which 0 denotes a short distanced connection that is to be removed. 
  rmShort <- apply(rmShort, c(1,2), as.numeric)
  save(file=expectFile, rmShort, roiDist)
}

######################################################################
## Load necessary libraries
suppressMessages(library(igraph))
# suppressMessages(library(RGtk2))
# suppressMessages(library(brainGraph))
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(gdata))
suppressMessages(library(beepr))
suppressMessages(library(parallel))
suppressMessages(library(data.table))
suppressMessages(library(psych))
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))
suppressMessages(library(broom))
suppressMessages(library(gtools))
suppressMessages(library(Hmisc))
suppressMessages(library(wle))
