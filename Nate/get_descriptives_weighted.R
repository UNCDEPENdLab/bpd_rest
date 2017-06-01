################Weighted graph descriptive stats
##degree dist across individuals
####read in package dependencies and custom functions
setwd("~/Box Sync/RS_BPD_graph")
basedir <- getwd()

source("Graph_util_redux.R")
library(igraph)
library(RGtk2)
library(brainGraph)
library(ggplot2)
library(dplyr)
library(gdata)
library(beepr)
library(parallel)
library(data.table)
library(psych)

#################################################################################
##Specifications to set up descriptives
nnodes <- 269  #varying this will require you to change the .txt coordinate file you read in for roiMat (should be 4 X nnodes txt file labelled:"x", "y", "z", "roi") and the masterlookup atlas
roiMat <- read.table("~/Box Sync/RS_BPD_graph/bb264coordinate_appended_shift_nate_culled.txt", header=FALSE, col.names=c("x", "y", "z", "roi"))
atlas <- read.csv("/Users/nth7/Box Sync/RS_BPD_graph/power269_masterlookup_shift_nate.csv", header = TRUE)
use.infomap <- 1
use.scotmi <- 0 
if(use.scotmi == 0){pipeline <- "pearson"} else{pipeline  <- "scotmi"}

#################################################################################
################################################################################
#######read in already processed Rdata files for faster run throughs:
#all graphs  including community assignment for infomap communities at 10% density
if(file.exists(paste0(basedir, "/cachedRfiles/allg.weighted.infomap10.", pipeline, ".RData")) == TRUE) {
  allg_infomap10 <- get(load(file = paste0(basedir, "/cachedRfiles/allg.weighted.infomap10.", pipeline, ".RData")))
}
#computed nodal metrics
if(file.exists(paste0(basedir, "/cachedRfiles/allmetrics.weighted.", pipeline, ".RData")) == TRUE) {
  allmetrics  <- get(load(paste0(basedir, "/cachedRfiles/allmetrics.weighted.", pipeline, ".RData")))
} 
#significant nodal comparisons
if(file.exists(paste0(basedir, "/cachedRfiles/node.metrics.weighted.",pipeline,".RData")) == TRUE) {
  node.metrics <- get(load(paste0(basedir, "/cachedRfiles/node.metrics.weighted.",pipeline,".RData")))
}

################################################################################
#################read in demographic info and combine with proper files in correct directory
SPECC_rest <- read.csv("SPECC_info_trimmed.csv", header = TRUE)
SPECC_rest$SPECC_ID <- as.character(SPECC_rest$SPECC_ID)
adjbase <- file.path(paste0("~/Box Sync/RS_BPD_graph/adjmats_269_", pipeline, "/"))
SPECC_rest$file <- NA_character_

for (i in 1:nrow(SPECC_rest)) {
  fname <- NULL
  if(use.scotmi == 0){
    #use regex for pearson fisherz
    if (SPECC_rest[i,"LunaMRI"] == 1) {
      #look for file under Luna ID
      fname <- list.files(path=adjbase, pattern=paste0(SPECC_rest[i, "Luna_ID"], "_20\\d{6}_adjmat",nnodes,"_pearson_fisherz.txt.gz"), ignore.case = TRUE, full.names=TRUE)
    } else {
      #look for file under SPECC ID
      fname <- list.files(path=adjbase, pattern=paste0(SPECC_rest[i, "SPECC_ID"], "_\\d{2}\\w{3}20\\d{2}_adjmat",nnodes,"_pearson_fisherz.txt.gz"), ignore.case = TRUE, full.names=TRUE)
    }} else{
      #use regex for scotmi
      if (SPECC_rest[i,"LunaMRI"] == 1) {
        #look for file under Luna ID
        fname <- list.files(path=adjbase, pattern=paste0(SPECC_rest[i, "Luna_ID"], "_20\\d{6}_S-raw.txt"), ignore.case = TRUE, full.names=TRUE)
      } else {
        #look for file under SPECC ID
        fname <- list.files(path=adjbase, pattern=paste0(SPECC_rest[i, "SPECC_ID"], "_\\d{2}\\w{3}20\\d{2}_S-raw.txt"), ignore.case = TRUE, full.names=TRUE)
      }}
  if (length(fname) != 1) {
    print(SPECC_rest[i,])
    stop("Cannot locate a file (or located more than 1) for the above record")
  } else {
    SPECC_rest[i,"file"] <- fname
  }
}

SPECC_rest <- filter(SPECC_rest, pr_over5mm <= .15)
#################################################################################
##graph degree distribution by subject 
# allg_infomap10 <- apply(allg_infomap10, 1, function(sub){
# 
#   set_graph_attr(sub, "ID", attributes(sub)$names)
# })
# 
# 
# get.graph.attribute(allg_infomap10[[79]])
# 
# allg_infomap10 <- lapply(allg_infomap10, function(g) {
#   browser()
#   set_graph_attr(g, "ID", attributes(g)$names)
# })

for (i in length(allg_infomap10)){
  allg_infomap10[[i]] <- set_graph_attr(allg_infomap10[[i]], "ID", names(allg_infomap10[i]))
}

pdf("SubjDegreeDist.weighted.pdf", width=8.5, height=11)
for (i in length(allg_infomap10)){
  plot_obj <- data.frame(degree(allg_infomap10[[i]]))
  colnames(plot_obj) <- "Degree"
  g <- ggplot(plot_obj, aes(Degree)) + geom_histogram() + ggtitle(paste0("Subject: ", names(allg_infomap10[i]))) + stat_bin(binwidth = 15)
  print(g)
  }
dev.off()



