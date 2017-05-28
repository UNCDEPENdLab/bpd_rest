##########RS_BPD_pipeline for WEIGHTED adjacency matrices
####read in package dependencies and custom functions
setwd("~/Box Sync/RS_BPD_graph")
basedir <- getwd()

source("calcGraph_weighted.R")
source("Graph_util_redux.R")
source("run_parse_deltacon.R")
source("wibw_module_degree.R")
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

#########NOTE: for creating adjacency matrices with an equal number of edges, in this script currently, we do not remove edges with a short euclidian distance.

#################################################################################
##Specifications to set up pipeline
nnodes <- 269  #varying this will require you to change the .txt coordinate file you read in for roiMat (should be 4 X nnodes txt file labelled:"x", "y", "z", "roi") and the masterlookup atlas
roiMat <- read.table("~/Box Sync/RS_BPD_graph/bb264coordinate_appended_shift_nate_culled.txt", header=FALSE, col.names=c("x", "y", "z", "roi"))
atlas <- read.csv("power269_masterlookup_shift_nate.csv", header = TRUE)
use.infomap <- 1
infomap.density <- 10
use.scotmi <- 0
if(use.scotmi == 0){pipeline <- "pearson"} else{pipeline  <- "scotmi"}
roi.dist <- 20
##metrics to run nodal comparisons on
metricstorun.nodal <- c("eigen.cent","degree", "closeness", "betweenness.node", "page.rank",  "part.coeff", "within.module.deg.zscore", "local.clustering", "gateway.coeff.btw", "gateway.coeff.degree", "between.module.deg.zscore")
metricstorun.global <- c("characteristic_path_length", "clustering_coefficient", "small_worldness", "modularity")

################################################################################
#######read in already processed Rdata files for faster run throughs:
#computed global metrics
if(file.exists(paste0(basedir, "/cachedRfiles/allmetrics.global.weighted.", pipeline, ".RData")) == TRUE) {
  allmetrics.global <- get(load(file = paste0(basedir, "/cachedRfiles/allmetrics.global.weighted.", pipeline, ".RData")))
} 
#computed nodal metrics
if(file.exists(paste0(basedir, "/cachedRfiles/allmetrics.weighted.", pipeline, ".RData")) == TRUE) {
  allmetrics  <- get(load(paste0(basedir, "/cachedRfiles/allmetrics.weighted.", pipeline, ".RData")))
} 
#significant global comaprisons

if(file.exists(paste0(basedir, "/cachedRfiles/global.metrics.weighted.",pipeline,".RData")) == TRUE) {
  global.metrics <- get(load(paste0(basedir, "/cachedRfiles/global.metrics.weighted.",pipeline,".RData")))
}
#significant nodal comparisons
if(file.exists(paste0(basedir, "/cachedRfiles/node.metrics.weighted.",pipeline,".RData")) == TRUE) {
  node.metrics <- get(load(paste0(basedir, "/cachedRfiles/node.metrics.weighted.",pipeline,".RData")))
}
#compiled significant nodal comparisons
if(file.exists(paste0(basedir, "/output.files/all.sig.nodal.weighted",pipeline,".rds")) == TRUE) {
  all.sig.nodal <- read.csv(paste0(basedir, "/output.files/all.sig.nodal.weighted",pipeline,".csv"))
}


# #total deltacon stats
# if(file.exists(paste0(basedir, "/output.files/deltacon_total_", pipeline, ".rds")) == TRUE) {
#   deltacon_total <- readRDS(paste0(basedir, "/output.files/deltacon_total_", pipeline, ".rds"))
# }
# #edge attribution deltacon values
# if(file.exists(paste0(basedir, "/output.files/edge_diffs_deltacon_", pipeline, ".rds")) == TRUE) {
#   edge_diffs_deltacon <- readRDS(paste0(basedir, "/output.files/edge_diffs_deltacon_", pipeline,".rds"))
# }
# #nodal attribution deltacon values
# if(file.exists(paste0(basedir, "/output.files/node_stats_deltacon_", pipeline, ".rds")) == TRUE) {
#   node_stats_deltacon <- readRDS(paste0(basedir, "/output.files/node_stats_deltacon_", pipeline, ".rds"))
# }


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

################################################################################
####Framewise displacement
#####filter subjects with over .20 brain volumes displaced .5mm or more

##In progress (make sure ics is mounted): get motion info (notes on how to implement this in RS notes folder in OneNote)
####this should include mean FD and max FD at the very least, standard script removes subjects with proportion of FD >.5mm of 20% or more 
#####currently no safeguard against very large head movements

## MH has now converted all SPECC MR directory names to all lower case to allow for match on case-sensitive filesystem
## and to make the naming consistent
# idfile <- "/gpfs/group/mnh5174/default/SPECC/SPECC_Participant_Info.csv"
# ##idinfo <- gdata::read.xls(idfile)
# idinfo <- read.csv(idfile)
# library(dplyr)
# options(dplyr.width=200)
# idinfo <- idinfo %>% rowwise() %>% mutate(mr_dir=ifelse(LunaMRI==1,
#                                                         paste0("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", Luna_ID, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d")), #convert to Date, then reformat YYYYMMDD
#                                                         paste0("/gpfs/group/mnh5174/default/SPECC/MR_Proc/", tolower(SPECC_ID), "_", tolower(format((as.Date(ScanDate, format="%Y-%m-%d")), "%d%b%Y")))))
# 
# #verify that mr_dir is present as expected
# idinfo$dirfound <- file.exists(idinfo$mr_dir)
# subset(idinfo, dirfound==FALSE)
# table(SPECC_rest[,c(3,5)])


##standard FD script
#SPECC_rest <- filter(SPECC_rest, pr_over5mm <= .15)
SPECC_rest <- filter(SPECC_rest, pr_over5mm <= .2)
table(SPECC_rest[,c(3,5)])

#describe(SPECC_rest[,c(1:6, 8)])

################################################################################
###setup Euclidean distance
nrois <- nrow(roiMat)
roiDist <- matrix(NA, nrow=nrois, ncol=nrois) 
#equiv for loop approach (easier to read)
#system.time(for (i in 1:nrois) {
for (j in 1:nrois) {
  roiDist[i,j] <- sqrt((roiMat[i,1] - roiMat[j,1])^2 + (roiMat[i,2]-roiMat[j,2])^2 + (roiMat[i,3] - roiMat[j,3])^2)
}
#})

############################quick QA
# hist(roiDist)
# vecDist <- roiDist[lower.tri(roiDist)]
# sum(vecDist < 20)/length(vecDist)

#roi.dist can be changed at the front end of the script
rmShort <- roiDist > roi.dist
#creates binary matrix, in which 0 denotes a short distanced connection that is to be removed. 
rmShort <- apply(rmShort, c(1,2), as.numeric)

################################################################################
##preallocate empty array to read adjacency matrices into
allmats <- array(NA, c(length(SPECC_rest[,1]), nnodes, nnodes), dimnames=list(id = SPECC_rest$SPECC_ID, roi1=paste0("V", 1:nnodes), roi2=paste0("V", 1:nnodes)))

#convert to matrix to be read as igraph obj and remove short distanced connections by multiplying by binary rmShort matrix
if(use.scotmi ==0){
  for (f in 1:length(SPECC_rest[,1])) {
    m <- as.matrix(read.table(as.character(SPECC_rest[f,7])))
   # m <- m * rmShort
    allmats[f,,] <- m
  }} else {
    for (f in 1:length(SPECC_rest[,1])) {
      m <- as.matrix(read.csv(as.character(SPECC_rest[f,7])))
      m <- rbind(array(NaN, ncol(m)), m)
      m[upper.tri(m)] <- t(m)[upper.tri(m)]
      #m <- m * rmShort
      allmats[f,,] <- m
    }}

################################################################################
#################create igraph object from adjacency matrix, include graph attr for mean strength of weighted edges and label nodes V1, V2,... 
allg <- apply(allmats, 1, function(sub) {
  g <- graph.adjacency(sub, mode="undirected", weighted=TRUE, diag=FALSE)
  g <- set_graph_attr(g, "mean.stren", mean(lowerTriangle(sub))) 
  V(g)$name <- paste0("V", 1:nrow(sub))
  g
})


save(allg, file = paste0(basedir,"/cachedRfiles/allg.weighted.", pipeline,".RData"))

################################################################################
#####POINT OF DEPARTURE FROM BINARY PIPELINE: AT THIS POINT IN BINARY PIPELINE GRAPHS ARE DENSITY THRESHOLDED
################################################################################
######overall strength differences
#Assign BPD and Age attributes to be used for group differences tests
attr(allg, "bpd") <- SPECC_rest[,"BPD"]
attr(allg, "age") <- SPECC_rest[,"AgeAtScan"]
attr(allg, "strength") <- sapply(allg, function(sub){
  sub$mean.stren
})
names(attributes(allg)$strength) <- NULL

stren.df <- data.frame(group = factor(attr(allg, "bpd"), levels = c(0,1), labels = c("control", "bpd")), age = attr(allg, "age"), strength = attr(allg, "strength"))
group.stren <- t.test(strength~group, stren.df, na.rm = TRUE) 
age.stren <- lm(strength~age*group, stren.df, na.action = "na.exclude")
summary(age.stren) ##non-sig p > .05

##remove negative correlations between nodes, if this is run on pearson, the number of edges remaining will be different across subjs
allg_noneg <- lapply(allg, function(g) {
  delete.edges(g, which(E(g)$weight < 0))
})

################################################################################
######COMMUNITY STRUCTURE
if (use.infomap == 1){
  mean.g.infomap <- readRDS(paste0(basedir, "/cachedRfiles/infomap_communitylist.rds"))
  commun.weight <- mean.g.infomap[[infomap.density]]
} else {
mean.graph <- apply(allmats, c(2,3), mean, na.rm = TRUE)
mean.g <- graph.adjacency(mean.graph, mode = "lower", weighted = TRUE, diag = FALSE)
mean.g <- delete.edges(mean.g, which(E(mean.g)$weight < 0))
commun.weight <- cluster_louvain(mean.g)
}


##########condense bad communities into community 99: "badcomm"
#V(mean.g.community[[10]])$community ##use this if using louvain

df <- as.data.frame(table(commun.weight$membership))
badcomm <- df$Var1[df$Freq < 4]
goodcomm <- as.numeric(df$Var1[df$Freq >= 4])
commun.weight$membership[commun.weight$membership %in% badcomm] <- length(goodcomm) + 1   #where highest community number is extraneous nodes
table(commun.weight$membership)
# plot(mean.g.infomap[[10]], mean.g)

#Assign community from mean graph back into subjects
allg_noneg <- lapply(allg_noneg, function(subj){
  V(subj)$community <- commun.weight$membership
  return(subj)
})

save(allg_noneg, file = paste0(basedir,"/cachedRfiles/allg_noneg.weighted.infomap10.", pipeline,".RData"))

################################################################################
##compute global metrics

if(!exists("allmetrics.global")) {
    allmetrics.global <- lapply(allg_noneg, function(subj) {
      calcGraph_weighted_global(subj)
    })
  
  save(allmetrics.global, file = paste0(basedir, "/cachedRfiles/allmetrics.global.weighted.", pipeline, ".RData"))
} 

################################################################################
##compute nodal metrics
#bombs out at foreach (line 273)

if(!exists("allmetrics")) {
  library(foreach)
  library(doSNOW)
  
  setDefaultClusterOptions(master="localhost") #move away from 10187 to avoid collisions
  clusterobj <- makeSOCKcluster(4)
  registerDoSNOW(clusterobj)
  
  allg_noneg <- lapply(1:length(allg_noneg), function(i) { allg_noneg[[i]]$pos <- i; return(allg_noneg[[i]]) } )
  #allg_noneg <- allg_noneg[-54]
  #allmetrics.nodal.weighted <- lapply(allg_noneg, function(subj) {
  allmetrics <- foreach(subj=iter(allg_noneg), .inorder = TRUE, .packages=c("brainGraph", "igraph")) %do% {
    gout <- tryCatch(calcGraph_weighted_nodal(subj), error=function(e) { print(e); cat("Died on subject: ", subj$pos, "\n"); return(NULL) })
    gout
  }

  save(allmetrics, file = paste0(basedir, "/cachedRfiles/allmetrics.weighted.", pipeline, ".RData"))
 
  } 

################################################################################
#IN PROGRESS
deg_0 <- array(NA, c(length(allmetrics.nodal.weighted), length(densities_desired), nnodes )) #3-D FALSE that's subjects x densities x nodes. 
for(s in 1:length(allmetrics)) {
  for(d in 1:length(densities_desired)){
    deg_0[s,d,] <- allmetrics[[s]][[d]]$degree == 0
  }
}
badnodes <- apply(deg_0, 2, function(mat) { colSums(mat) }) 