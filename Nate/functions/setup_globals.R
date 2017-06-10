#################################################################################
## Specifications to set up pipeline
## This should be sourced by scripts that import files or compute metrics
nnodes <- 269  #varying this will require you to change the .txt coordinate file you read in for roiMat (should be 4 X nnodes txt file labelled:"x", "y", "z", "roi") and the masterlookup atlas
parcellation <- "power269"
roiMat <- read.table(file.path(basedir, "bb264coordinate_appended_shift_nate_culled.txt"), header=FALSE, col.names=c("x", "y", "z", "roi"))
atlas <- read.csv(paste0(parcellation, "_atlas.csv"), header = TRUE, stringsAsFactors = FALSE) #eventually would be nice to have a lockstep name here and re-use parcellationName
use.infomap <- 1
use.scotmi <- 0
if (use.scotmi == 0) { pipeline <- "pearson_fisherz" } else { pipeline  <- "scotmi" }
roi.dist <- 20 #distance in mm to zero connections when setting up graphs (cf. Power 2011)
metricstorun.nodal <- c("eigen.cent","degree", "closeness", "betweenness.node", "page.rank",  "part.coeff", "within.module.deg.zscore", "local.clustering", "gateway.coeff.btw", "gateway.coeff.degree", "between.module.deg.zscore")
metricstorun.global <- c("characteristic_path_length", "clustering_coefficient", "small_worldness", "modularity")
#adjmats_base <- file.path(basedir, "adjmats")
adjmats_base <- file.path("/Users/michael/Data_Analysis/bpd_rest", "adjmats") #should be above, but working from a local copy to avoid messing up existing things
densities_desired <- seq(.01, .2, .01) #used globally for binary graphs

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

## Load necessary libraries
suppressMessages(library(igraph))
suppressMessages(library(RGtk2))
suppressMessages(library(brainGraph))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(gdata))
suppressMessages(library(beepr))
suppressMessages(library(parallel))
suppressMessages(library(data.table))
suppressMessages(library(psych))
suppressMessages(library(tidyr))