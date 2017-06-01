################Binary graph descriptive stats
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
library(data.table)

#################################################################################
##Specifications to set up descriptives
nnodes <- 269  #varying this will require you to change the .txt coordinate file you read in for roiMat (should be 4 X nnodes txt file labelled:"x", "y", "z", "roi") and the masterlookup atlas
roiMat <- read.table("~/Box Sync/RS_BPD_graph/bb264coordinate_appended_shift_nate_culled.txt", header=FALSE, col.names=c("x", "y", "z", "roi"))
use.infomap <- 1
use.scotmi <- 0 
if(use.scotmi == 0){pipeline <- "pearson"} else{pipeline  <- "scotmi"}
atlas <- read.csv("power269_masterlookup_shift_nate.csv", header = TRUE)
roi.dist <- 20
mean.g.infomap <- readRDS(paste0(basedir, "/cachedRfiles/infomap_communitylist.rds"))

#################################################################################
################################################################################
#######read in already processed Rdata files for faster run throughs:
#all graphs at 1-20% density thresholds including community assignment at infomap communities at 10% density
if(file.exists(paste0(basedir, "/cachedRfiles/allg_density_infomap10.", pipeline, ".RData")) == TRUE) {
  allg_density_infomap10 <- get(load(file = paste0(basedir, "/cachedRfiles/allg_density_infomap10.", pipeline, ".RData")))
}
#computed nodal metrics
if(file.exists(paste0(basedir, "/cachedRfiles/allmetrics.binary.", pipeline, ".RData")) == TRUE) {
  allmetrics  <- get(load(paste0(basedir, "/cachedRfiles/allmetrics.binary.", pipeline, ".RData")))
} 
#significant nodal comparisons
if(file.exists(paste0(basedir, "/cachedRfiles/node.metrics.binary.",pipeline,".RData")) == TRUE) {
  node.metrics <- get(load(paste0(basedir, "/cachedRfiles/node.metrics.binary.",pipeline,".RData")))
}
#compiled significant nodal comparisons
if(file.exists(paste0(basedir, "/output.files/all.sig.nodal.",pipeline,".rds")) == TRUE) {
  all.sig.nodal <- read.csv(paste0(basedir, "/output.files/all.sig.nodal.",pipeline,".csv"))
}

densities_desired <- seq(.01, .2, .01) 

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


################################################################################
#####filter subjects with over .15 brain volumes displaces .5mm or more
table(SPECC_rest[,c(3,5)])
SPECC_rest <- filter(SPECC_rest, pr_over5mm <= .15)
table(SPECC_rest[,c(3,5)])

describe(SPECC_rest[,c(1:6, 8)])
################################################################################
###setup Euclidean distance
nrois <- nrow(roiMat)
roiDist <- matrix(NA, nrow=nrois, ncol=nrois) 
#equiv for loop approach (easier to read)
system.time(for (i in 1:nrois) {
  for (j in 1:nrois) {
    roiDist[i,j] <- sqrt((roiMat[i,1] - roiMat[j,1])^2 + (roiMat[i,2]-roiMat[j,2])^2 + (roiMat[i,3] - roiMat[j,3])^2)
  }
})

############################quick QA
# hist(roiDist)
# vecDist <- roiDist[lower.tri(roiDist)]
# sum(vecDist < 20)/length(vecDist)

rmShort <- roiDist > roi.dist
rmShort <- apply(rmShort, c(1,2), as.numeric)

################################################################################
##preallocate empty array to read adjacency matrices into
allmats <- array(NA, c(length(SPECC_rest[,1]), nnodes, nnodes), dimnames=list(id = SPECC_rest$SPECC_ID, roi1=paste0("V", 1:nnodes), roi2=paste0("V", 1:nnodes)))

#convert to matrix to be read as igraph obj and remove short distanced connections by multiplying by binary rmShort matrix
if(use.scotmi ==0){
  for (f in 1:length(SPECC_rest[,1])) {
    m <- as.matrix(read.table(as.character(SPECC_rest[f,7])))
    m <- m * rmShort
    allmats[f,,] <- m
  }} else {
    for (f in 1:length(SPECC_rest[,1])) {
      m <- as.matrix(read.csv(as.character(SPECC_rest[f,7])))
      m <- rbind(array(NaN, ncol(m)), m)
      m[upper.tri(m)] <- t(m)[upper.tri(m)]
      m <- m * rmShort
      allmats[f,,] <- m
    }}

#################################################################################
#trim communities to density of 10
df <- as.data.frame(table(mean.g.infomap[[10]]$membership))
badcomm <- df$Var1[df$Freq < 4]
goodcomm <- as.numeric(df$Var1[df$Freq >= 4])
mean.g.infomap[[10]]$membership[mean.g.infomap[[10]]$membership %in% badcomm] <- length(goodcomm) + 1   #where highest community number is extraneous nodes
table(mean.g.infomap[[10]]$membership)

module.membership <- mean.g.infomap[[10]]$membership




#################################################################################
######IN PROGRESS: high degree nodes in density of 10
den10deg <- sapply(allmetrics, function(sub){
  return(sub[[10]]$degree)
})
  
###calculate degree based on mean graph  
mean.graph <- apply(allmats, c(2,3), mean, na.rm = TRUE)
mean.g <- graph.adjacency(mean.graph, mode = "lower", weighted = TRUE, diag = FALSE)
mean.g <- delete.edges(mean.g, which(E(mean.g)$weight < 0))

mean.g.deg10 <- list()
mean.g.deg10[["degree"]] <- apply(den10deg, 1, mean)

topdeg10 <- plotMetricQuant(mean.g.deg10, .90, "degree", atlas = atlas)

deg10_node <- NodeFile(mean.g, atlas = atlas, nodestp = topdeg10, name = "practice1", nodestat = mean.g.deg10[["degree"]], community = module.membership)

names(topdeg10)

#################################################################################
######plot nodes with high mean metrics
density <- .10
metric <- "Betweenness Centrality"  #still have to change sub$metric in metric.node call

metric.node <- sapply(allmetrics_den10, function(sub){

  return(sub$betweenness.node)
})


metric.mean <- apply(metric.node, 1, mean)
# metric.sd <- apply(metric.node, 1, sd)
# mean.sd.metric <- cbind(metric.mean, metric.sd)

atlas <- read.csv("/Users/nth7/Box Sync/RS_BPD_graph/power269_masterlookup_shift_nate.csv", header = TRUE)
atlas$label <- paste0("V",atlas$vname, "_", atlas$anat_label)

metric.sort <- factor(atlas$label, levels = atlas$label[order(metric.mean)])
df.plot <- data.frame(metric = metric.mean, roi = atlas$label, roisort = metric.sort)

library(ggplot2)
pdf(paste0(metric, "_Total.pdf"), width=8.5, height=35)
ggplot(df.plot, aes(x = roisort, y = metric)) +geom_point() + coord_flip() + labs(title = paste0("Overall ", metric, " by node"), y = metric, x = "Node Number and Anatomical Label")
dev.off()


node.file.betw <- NodeFile( allmetrics_den10[["origgraph"]])
#################################################################################
####IN PROGRESS split by group
metric <- "Degree"  #still have to change sub$metric in metric.node call

metric.node <- sapply(allmetrics_den10, function(sub){
  return(sub$degree)
})

ID <- SPECC_rest$SPECC_ID
BPD <- SPECC_rest$BPD

metric.node <- data.frame(cbind(BPD, t(metric.node)))

group.agg <- aggregate( .~BPD, metric.node, mean)
rownames(group.agg) <- c("Control", "BPD")
group.agg <- group.agg[,-1]
atlas <- read.csv("/Users/natehall/Box Sync/RS_BPD_graph/power269_masterlookup_shift_nate.csv", header = TRUE)
atlas$label <- paste0("V",atlas$vname, "_", atlas$anat_label)
group.agg

control.metric <- group.agg["Control",]
row.names(control.metric) <- NULL
control.sort <- factor(atlas$label, levels = atlas$label[order(control.metric)])
# control.plot <- cbind(control.metric, atlas$label, control.sort)
control.plot <- data.frame(metric = t(control.metric), roi = atlas$label, roisort = control.sort)

gg.control <- ggplot(df.plot, aes(x = roisort, y = metric)) +geom_point() + coord_flip() + labs( y = metric, x = "Node Number and Anatomical Label")
gg.control
# metric.node.long <- reshape(data = metric.node, direction = "long", varying = seq(V1, V269, 1)) 

bpd.metric <- group.agg["BPD",]
row.names(bpd.metric) <- NULL
bpd.sort <- factor(atlas$label, levels = atlas$label[order(bpd.metric)])
# control.plot <- cbind(control.metric, atlas$label, control.sort)
bpd.plot <- data.frame(metric = t(bpd.metric), roi = atlas$label, roisort = bpd.sort)

gg.bpd <- ggplot(df.plot, aes(x = roisort, y = metric)) +geom_point() + coord_flip() + labs( y = metric, x = "Node Number and Anatomical Label")
gg.bpd

library(cowplot)
# pdf(paste0(metric, "_Total_groupSplit.pdf"), width=10, height=60)
gg.split <- plot_grid(gg.control, gg.bpd), labels = c("Control", "BPD"), nrow = 2, align = "v")
# dev.off()

save_plot(paste0(metric, "_Total_groupSplit.pdf"), gg.split,
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
           base_height = 30, base_width = 8.5, limitsize = FALSE
)

metric.mean <- apply(metric.node, 2, mean)
# metric.sd <- apply(metric.node, 1, sd)
# mean.sd.metric <- cbind(metric.mean, metric.sd)

atlas <- read.csv("/Users/nth7/Box Sync/RS_BPD_graph/power269_masterlookup_shift_nate.csv", header = TRUE)
atlas$label <- paste0("V",atlas$vname, "_", atlas$anat_label)

metric.sort <- factor(atlas$label, levels = atlas$label[order(metric.mean)])
df.plot <- data.frame(metric = metric.mean, roi = atlas$label, roisort = metric.sort)

library(ggplot2)
pdf(paste0(metric, "_Total.pdf"), width=8.5, height=35)
ggplot(df.plot, aes(x = roisort, y = metric)) +geom_point() + coord_flip() + labs(title = paste0("Overall ", metric, " by node"), y = metric, x = "Node Number and Anatomical Label")
dev.off()


#################################################################################
######IN PROGRESS: degree dists in density of 10
#all graphs  including community assignment for infomap communities at 10% density
if(file.exists(paste0(basedir, "/cachedRfiles/allg.weighted.infomap10.", pipeline, ".RData")) == TRUE) {
  allg_infomap10 <- get(load(file = paste0(basedir, "/cachedRfiles/allg.weighted.infomap10.", pipeline, ".RData")))
}

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


#################################################################################
##IN PROGRESS: graph degree distribution by subject and by density
allg_density <- lapply(allg_density_infomap10, function(subj) {
  subj_graphs <- lapply(1:length(densities_desired), function(d) {
    V(subj[[d]])$community <- mean.g.infomap[[10]]$membership
    return(subj[[d]])
  })
})

length(allmetrics)
pdf("SubjDegreeDist.pdf", width=8.5, height=11)
for (subj in length(allmetrics)){
 for(den in length(subj)){
   browser()
   plot <- den$degree
   hist(plot)
   
 }
}
#################################################################################

for (subj in divideBySubject) {
  browser()
  subj <- subset(subj, blockcode=="animalblock")
  #cat("Subject is: ", subj$subject[1], "RT range is: ", range(subj$latency), "\n")
  #cat("Subject is: ", subj$subject[1], "RT range (correct) is: ", range(subj$latency[which(subj$correct==1)]), "\n")
  cat("Subject is: ", subj$subject[1], "Accuracy %: ", round(100*sum(subj$correct)/length(subj$correct), 2), "\n")
  g <- ggplot(subset(subj, correct==1), aes(x=latency)) + geom_histogram(binwidth = 250) + ggtitle(paste0("Subject: ", subj$subject[1]))
  #xlim(c(0,3700))
  print(g)
}
dev.off()
