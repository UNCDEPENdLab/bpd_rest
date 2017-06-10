##########RS_BPD_pipeline PCA 
####read in package dependencies and custom functions
setwd("~/Box Sync/RS_BPD_graph")
basedir <- getwd()

source("calcGraph_binary.R")
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
library(tidyr)

##WARNING: currently reading FD info requires ics to be mounted on your computer, make sure either ics is mounted or move FD files to another folder

#################################################################################
##Specifications to set up pipeline
nnodes <- 269  #varying this will require you to change the .txt coordinate file you read in for roiMat (should be 4 X nnodes txt file labelled:"x", "y", "z", "roi") and the masterlookup atlas
roiMat <- read.table("~/Box Sync/RS_BPD_graph/bb264coordinate_appended_shift_nate_culled.txt", header=FALSE, col.names=c("x", "y", "z", "roi"))
atlas <- read.csv("power269_masterlookup_shift_nate.csv", header = TRUE)
use.infomap <- 1
use.scotmi <- 0
if(use.scotmi == 0){pipeline <- "pearson"} else{pipeline  <- "scotmi"}
roi.dist <- 20
metricstorun.nodal <- c("eigen.cent","degree", "closeness", "betweenness.node", "page.rank",  "part.coeff", "within.module.deg.zscore", "local.clustering", "gateway.coeff.btw", "gateway.coeff.degree", "between.module.deg.zscore")
metricstorun.global <- c("characteristic_path_length", "clustering_coefficient", "small_worldness", "modularity")
#FD motion specs
ics_dir <- "/mnt/ics"
thresh <- 0.5
spikes.exclude <- 0.2
max.exclude <- 10

################################################################################
#######read in already processed Rdata files for faster run throughs:
#all graphs at 1-20% density thresholds
if(file.exists(paste0(basedir, "/cachedRfiles/allg_density.", pipeline, ".RData")) == TRUE) {
  allg_density <- get(load(file = paste0(basedir, "/cachedRfiles/allg_density.", pipeline, ".RData")))
}
#computed global metrics
if(file.exists(paste0(basedir, "/cachedRfiles/allmetrics.global.binary.", pipeline, ".RData")) == TRUE) {
  allmetrics.global <- get(load(file = paste0(basedir, "/cachedRfiles/allmetrics.global.binary.", pipeline, ".RData")))
} 
#computed nodal metrics
if(file.exists(paste0(basedir, "/cachedRfiles/allmetrics.binary.", pipeline, ".RData")) == TRUE) {
  allmetrics  <- get(load(paste0(basedir, "/cachedRfiles/allmetrics.binary.", pipeline, ".RData")))
} 
#significant global comaprisons

if(file.exists(paste0(basedir, "/cachedRfiles/global.metrics.binary.",pipeline,".RData")) == TRUE) {
  global.metrics.bin <- get(load(paste0(basedir, "/cachedRfiles/global.metrics.binary.",pipeline,".RData")))
}
#significant nodal comparisons
if(file.exists(paste0(basedir, "/cachedRfiles/node.metrics.binary.",pipeline,".RData")) == TRUE) {
  node.metrics <- get(load(paste0(basedir, "/cachedRfiles/node.metrics.binary.",pipeline,".RData")))
}
#compiled significant nodal comparisons
if(file.exists(paste0(basedir, "/output.files/all.sig.nodal.",pipeline,".rds")) == TRUE) {
  all.sig.nodal <- read.csv(paste0(basedir, "/output.files/all.sig.nodal.",pipeline,".csv"))
}
#total deltacon stats
if(file.exists(paste0(basedir, "/output.files/deltacon_total_", pipeline, ".rds")) == TRUE) {
  deltacon_total <- readRDS(paste0(basedir, "/output.files/deltacon_total_", pipeline, ".rds"))
}
#edge attribution deltacon values
if(file.exists(paste0(basedir, "/output.files/edge_diffs_deltacon_", pipeline, ".rds")) == TRUE) {
  edge_diffs_deltacon <- readRDS(paste0(basedir, "/output.files/edge_diffs_deltacon_", pipeline,".rds"))
}
#nodal attribution deltacon values
if(file.exists(paste0(basedir, "/output.files/node_stats_deltacon_", pipeline, ".rds")) == TRUE) {
  node_stats_deltacon <- readRDS(paste0(basedir, "/output.files/node_stats_deltacon_", pipeline, ".rds"))
}

################################################################################
#################read in demographic info and combine with proper files in correct directory
SPECC_rest <- read.csv("SPECC_info_trimmed.csv", header = TRUE)

SPECC_rest$SPECC_ID <- as.character(SPECC_rest$SPECC_ID)
adjbase <- file.path(paste0(basedir,"/adjmats_269_", pipeline, "/"))
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

SPECC_rest <- SPECC_rest[order(SPECC_rest$SPECC_ID),] 

################################################################################
####Framewise displacement scrubbing
#####filter subjects with high prop FD over given threshold and with large maximum movements
FDInfo <- MotionInfoSPECC(ics_dir, thresh, spikes.exclude, max.exclude)
####6/8/17 come back to 083 and 136.
FDInfo$fd.info <- FDInfo$fd.info[-which(FDInfo$fd.info == "083DB"),]
FDInfo$fd.info <- FDInfo$fd.info[-which(FDInfo$fd.info == "142MW"),]
FDInfo <- FDInfo$fd.info 
FDInfo$Subj <- as.character(FDInfo$Subj)

if (all(FDInfo$Subj==SPECC_rest$SPECC_ID)) {cat("FDInfo and SPECC_rest have matching subjects")} else {cat("FDInfo and SPECC_rest do not have matching subjects")}

SPECC_rest$FD.exclude <- as.numeric(FDInfo$maxExclude | FDInfo$thresh.exclude)
#final FD cuts, table, and dexcriptives
SPECC_rest <- SPECC_rest[which(SPECC_rest$FD.exclude == 0),]
table(SPECC_rest[,c(3,5)])
describe(SPECC_rest[,c(1:6, 8)])


if(file.exists(paste0(basedir, "/SPECC_rest.csv")) == FALSE) {
  write.csv(SPECC_rest, file = "SPECC_rest.csv")
}

################################################################################
###setup Euclidean distance
nrois <- nrow(roiMat)
roiDist <- matrix(NA, nrow=nrois, ncol=nrois) 

#system.time(for (i in 1:nrois) {
for (i in 1:nrois) {
  for (j in 1:nrois) {
    roiDist[i,j] <- sqrt((roiMat[i,1] - roiMat[j,1])^2 + (roiMat[i,2]-roiMat[j,2])^2 + (roiMat[i,3] - roiMat[j,3])^2)
  }
}#)

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

#Read in raw adjacency files and convert to matrix to be read as igraph obj and remove short distanced connections by multiplying by binary rmShort matrix
if(use.scotmi ==0){
  for (subj in 1:length(SPECC_rest[,1])) {
    m <- as.matrix(read.table(as.character(SPECC_rest[subj,7])))
    m <- m * rmShort
    allmats[subj,,] <- m
  }} else {
    for (subj in 1:length(SPECC_rest[,1])) {
      m <- as.matrix(read.csv(as.character(SPECC_rest[subj,7])))
      m <- rbind(array(NaN, ncol(m)), m)
      m[upper.tri(m)] <- t(m)[upper.tri(m)]
      m <- m * rmShort
      allmats[subj,,] <- m
    }}

################################################################################
#################create igraph object from adjacency matrix and label nodes V1, V2,...
allg <- apply(allmats, 1, function(sub) {
  g <- graph.adjacency(sub, mode="undirected", weighted=TRUE, diag=FALSE)
  V(g)$name <- paste0("V", 1:nrow(sub))
  g
})

##remove negative correlations between nodes, if this is run on pearson, the number of edges remaining will be different across subjs 
allg_noneg <- lapply(allg, function(g) {
  delete.edges(g, which(E(g)$weight < 0))
})

#check num of edges across subjs
for (subj in 1:length(allg)){
  print(length(E(allg_noneg[[subj]])))
}

str(allg$`001RA`)

################################################################################
#DENSITY THRESHOLDING: binarize and threshold graphs at densities ranging from 1-20%
densities_desired <- seq(.01, .2, .01) 

library(foreach)
library(doSNOW)

if(!exists("allg_density")){
  setDefaultClusterOptions(master="localhost") #move away from 10187 to avoid collisions
  clusterobj <- makeSOCKcluster(4)
  registerDoSNOW(clusterobj)
  
  allg_density <- foreach(g=allg_noneg, .packages=c("igraph", "tidyr")) %dopar% {
    nnodes <- length(V(g))
    maxedges <- (nnodes*(nnodes-1))/2
    
    #much less clunky version than the density thresholding below
    dgraphs <- lapply(densities_desired, function(d) {
      #Obtains desired density given graph diameter
      weights <- sort(E(g)$weight, decreasing=TRUE)
      threshold <- weights[length(V(g))*(length(V(g))-1)/2 * d]
      gthresh <- delete.edges(g, which(E(g)$weight < threshold))
      gthresh <- remove.edge.attribute(gthresh, "weight")
      return(gthresh)
    })
    
    
    names(dgraphs) <- paste0("d", densities_desired)
    return(dgraphs)
  }
  save(allg_density, file = paste0(basedir,"/cachedRfiles/allg_density.", pipeline,".RData"))
}

#each element of allg_density is a list of 20 binary graphs for that subject at 1-20% density

################################################################################
######COMMUNITY STRUCTURE
####rerun outside of allg_density for use in community detection
###density thresholding and community detection (louvain)

#allmats_log <- apply(allmats, c(1,2,3), function(x) { log(x+.05)})
mean.graph <- apply(allmats, c(2,3), mean, na.rm = TRUE)
mean.g <- graph.adjacency(mean.graph, mode = "lower", weighted = TRUE, diag = FALSE)
mean.g <- delete.edges(mean.g, which(E(mean.g)$weight < 0))
mean.g.adj <- as_adjacency_matrix(mean.g, attr = "weight")

# hist(E(mean.g)$weight)
# E(mean.g)
# dev.off()

if (use.infomap == 1){
  mean.g.infomap <- readRDS(paste0(basedir, "/cachedRfiles/infomap_communitylist.rds"))
} else{
  ##output will be 20 binary graphs AVERAGED ACROSS SUBJECTS (at 1-20% density) with community membership
  mean.g.community <- lapply(densities_desired, function(d) {
    weights <- sort(E(mean.g)$weight, decreasing=TRUE)
    threshold <- weights[length(V(mean.g))*(length(V(mean.g))-1)/2 * d]
    gthresh <- delete.edges(mean.g, which(E(mean.g)$weight < threshold))
    gthresh <- remove.edge.attribute(gthresh, "weight")
    louv.graph <- cluster_louvain(gthresh)
    V(gthresh)$community <- louv.graph$membership
    
    return(gthresh)
  })
}



##########condense bad communities into the last community: community-dependent results from the last module should not be interpreted
#V(mean.g.community[[10]])$community ##use this if using louvain

df <- as.data.frame(table(mean.g.infomap[[10]]$membership))
badcomm <- df$Var1[df$Freq < 4]
goodcomm <- as.numeric(df$Var1[df$Freq >= 4])
mean.g.infomap[[10]]$membership[mean.g.infomap[[10]]$membership %in% badcomm] <- length(goodcomm) + 1   #where highest community number is extraneous nodes
mean.g$community <- mean.g.infomap[[10]]$membership
table(mean.g.infomap[[10]]$membership)
# plot(mean.g.infomap[[10]], mean.g)

#Assign community from mean graph at a given density (10) back into subjects
#each element of allg_density is a list of subjects with 20 binary graphs for each subject at 1-20% density, now with communities assigned to them
allg_density <- lapply(allg_density, function(subj) {
  subj_graphs <- lapply(1:length(densities_desired), function(d) {
    V(subj[[d]])$community <- mean.g.infomap[[10]]$membership
    return(subj[[d]])
  })
})


mean.g.info.df <- data.frame(mean.g.infomap[[10]]$membership)
colnames(mean.g.info.df) <- NULL
#if(!file.exists(paste0(basedir,"/cachedRfiles/allg_density_infomap10.RData"))){save(allg_density, file = paste0(basedir,"/cachedRfiles/allg_density_infomap10.RData"))}

#example: what is community assignment for subject [[1]] at density [[20]]..will be the same across densities
#V(allg_density[[1]][[20]])$community

################################################################################
##compute global metrics 
if(!exists("allmetrics.global")) {
  allmetrics.global <- foreach(subj=allg_density, .packages = c("igraph", "brainGraph")) %dopar% {
    #for (subj in allg_density) { #put here for more fine-grained debugging
    #require(igraph)
    dl <- lapply(subj, function(dgraph) {
      calcGraph_binary_global(dgraph)
    })
    names(dl) <- paste0("d", densities_desired)
    dl
  }
  save(allmetrics.global, file = paste0(basedir,"/cachedRfiles/allmetrics.global.binary.", pipeline,".RData"))
} 

################################################################################
###compute nodal metrics 
if(!exists("allmetrics")) {
  registerDoSEQ()
  # registerDoSNOW(clusterobj)
  allmetrics <- foreach(subj=allg_density, .packages = c("igraph", "brainGraph")) %dopar% {
    #for (subj in allg_density) { #put here for more fine-grained debugging
    #require(igraph)
    dl <- lapply(subj, function(dgraph) {
      
      calcGraph_binary_nodal(dgraph)
    })
    names(dl) <- paste0("d", densities_desired)
    dl
  }
  save(allmetrics, file = paste0(basedir,"/cachedRfiles/allmetrics.binary.", pipeline,".RData"))
} 

###################################
####PCA pearson
node.table.pearson <- data.frame()

names(allmetrics) <- SPECC_rest$SPECC_ID #assumes elements of allmetrics in same order as this DF
attr(allmetrics, "bpd") <- SPECC_rest[,"BPD"]
attr(allmetrics, "age") <- SPECC_rest[,"AgeAtScan"]

node.table.pearson <- do.call(rbind, lapply(1:length(allmetrics), function (s) {
  subj <- allmetrics[[s]]
  
  ddf <- lapply(1:length(subj), function(d) {
    density <- subj[[d]]
    
    metrics_to_grab <- names(density)[!names(density) %in% c("origgraph")]
    nvec <- lapply(1:length(metrics_to_grab), function(m) { 
      density[[metrics_to_grab[m]]]
    })
    
    names(nvec) <- metrics_to_grab
    ndf <- as.data.frame(do.call(cbind, nvec))
    ndf$density <- names(subj)[d]
    ndf$node <- atlas$vname
    ndf$anat_label <- atlas$anat_label
    ndf$scotmi <- 0
    return(ndf)
  })
  
  #names(ddf) <- names(subj)
  
  
  ddf <- do.call(rbind, ddf)
  ddf$subj <- names(allmetrics)[s]
  ddf$bpd <- attr(allmetrics, "bpd")[s]
  ddf$age <- attr(allmetrics, "age")[s]
  
  return(ddf)
  
}))
metrics <- names(node.table.pearson[, c(2:6, 8:12)])
######################################################################
#PCA on latent degree/ eigenvector across densities

metrics.raw <- node.table.pearson[,c(metrics, "subj", "node", "density")]
#metrics.raw.noNA <- metrics.raw[,c(1:7)]
metrics.raw$density_num <- as.numeric(sub("d(.*)", "\\1", metrics.raw$density, perl=TRUE))

metrics.raw <- metrics.raw %>% select(subj, node, density, density_num, degree, eigen.cent, betweenness.node, within.module.deg.zscore, between.module.deg.zscore) %>%
  filter(density_num > .04) %>% select(-density_num) %>% gather(key="variable", value="value", degree, eigen.cent, betweenness.node, within.module.deg.zscore, between.module.deg.zscore)


#dfm <- melt(metrics.raw, id.vars = c("subj", "node", "density"))

metrics.raw_pca <- dcast(metrics.raw, subj + node ~ density + variable, value.var = "value")

metrics.pca <- prcomp(select(metrics.raw_pca, -subj, -node), center = TRUE, scale. = TRUE)


pcaout1 <- pca(select(metrics.raw_pca, -subj, -node), nfactors=1)
pcaout2 <- pca(select(metrics.raw_pca, -subj, -node), nfactors=2)
pcaout3 <- pca(select(metrics.raw_pca, -subj, -node), nfactors=3)
pcaout4 <- pca(select(metrics.raw_pca, -subj, -node), nfactors=4)
pcaout5 <- pca(select(metrics.raw_pca, -subj, -node), nfactors=5)
pcaout6 <- pca(select(metrics.raw_pca, -subj, -node), nfactors=6)

str(pcaout2$scores)

print(pcaout4$loadings, cutoff = 0.5)

pcasolution <- data.frame(metrics.pca$x[,1:4])
names(pcasolution) <- c("central", "between.node", "within.mod", "between.mod")

toanalyze <- cbind(select(metrics.raw_pca, subj, node), pcasolution)

merge.bpdage <- SPECC_rest[,c("SPECC_ID", "BPD", "AgeAtScan")]
colnames(merge.bpdage) <- c("subj", "BPD", "Age")

toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "subj")
write.csv(toanalyze, file = paste0(basedir, "/output.files/nodal.pca.toanalyze.",pipeline,".csv"))

##################################################################################################
##############
####run tests on nodal metrics collapsed across densities
abst <- 2.64 ##about p < .005

metrics.toanalyze <- names(select(toanalyze, -subj, -node, -BPD, -Age))
nodes <- unique(select(toanalyze, -central, -between.node, -within.mod, -between.mod, -subj, -BPD, -Age))
nodes <- as.numeric(nodes[,1])
atlas <- read.csv("power269_masterlookup_shift_nate.csv", header = TRUE)
# length(nodes)
# max(nodes)
#OR: nodes <- 1:269
sigres <- 1
results <- list()

for (m in metrics.toanalyze) {
    for (n in nodes) {
      # browser()
      thismetric <- toanalyze[toanalyze$node ==n, c(m, "BPD", "Age")]
      thismetric$BPD <- factor(thismetric$BPD, levels = c(0,1), labels = c("control", "BPD"))
      colnames(thismetric) <- c("metric", "group", "age")
      
      node.test <- tryCatch(t.test(metric~group, thismetric, na.rm = TRUE), error = function(errorname) { print(errorname); return(NULL) })
      if (is.null(node.test)) { message("Error occurred for metric: ", m, " and node: ", n) }

      age.test <- tryCatch(lm(metric~age*group, thismetric, na.action = "na.exclude"), error=function(e) { print(e); return(NULL) })
      if (is.null(age.test)) { pvec <- NULL
      } else { pvec <- broom::tidy(age.test)$p.value[-1] } #p-values of age, bpd, and age x bpd. -1 to drop off intercept
      
      if ((!is.null(node.test) && !is.nan(node.test$statistic) && abs(node.test$statistic) > abst) ||
          (!is.null(age.test) && !all(is.nan(pvec)) && any(pvec < .01))) {
        
        results[[sigres]] <- list(nodename=as.character(atlas$anat_label[atlas$vname == n]), ttest=node.test, 
                                  metric = m, agetest=age.test, nodenum=n)
        sigres <- sigres+1
      }
    }
}

#results

####BPD main effects 
results.ttest <- lapply(results, function(node){
  df <- broom::tidy(node$ttest)
  df$nodename <- node$nodename
  df$nodenum <- node$nodenum
  df$metric <- node$metric
  return(df)
})

results.ttest.df <- do.call(rbind, results.ttest)
a <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "central" & p.value < .01)
b <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "between.node" & p.value < .01)
c <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "within.mod" & p.value < .01)
d <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "between.mod" & p.value < .01)

####Age and agexbpd lm 
results.agetest <- lapply(results, function(node){
  df <- broom::tidy(node$agetest)
  df$nodename <- node$nodename
  df$nodenum <- node$nodenum
  df$metric <- node$metric
  return(df)
})
results.agetest.df <- do.call(rbind, results.agetest)
###age main effects
e <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "central" & term == "age" & p.value < .01)
f <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "between.node" & term == "age" & p.value < .01)
g <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "within.mod" & term == "age" & p.value < .01)
h <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "between.mod" & term == "age" & p.value < .01)

###age x bpd interactions
i <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "central" & term == "age:groupBPD" & p.value < .01)
j <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "between.node" & term == "age:groupBPD" & p.value < .01)
k <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "within.mod" & term == "age:groupBPD" & p.value < .01)
l <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "between.mod" & term == "age:groupBPD" & p.value < .01)

library(gtools)
all.sig.nodal.pca <- bind_rows(a,b,c,d,e,f,g,h,i,j,k,l)
write.csv(all.sig.nodal.pca, file = paste0(basedir, "/output.files/all.sig.nodal.pca.",pipeline,".csv"))

