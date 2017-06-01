##########RS_BPD_pipeline
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

################################################################################
####Framewise displacement
#####filter subjects with over .20 brain volumes displaced .5mm or more

##In progress (make sure ics is mounted): 
## MH has now converted all SPECC MR directory names to all lower case to allow for match on case-sensitive filesystem
## and to make the naming consistent
idfile <- "/gpfs/group/mnh5174/default/SPECC/SPECC_Participant_Info.csv"
##idinfo <- gdata::read.xls(idfile)
idinfo <- read.csv(idfile)
library(dplyr)
options(dplyr.width=200)
idinfo <- idinfo %>% rowwise() %>% mutate(mr_dir=ifelse(LunaMRI==1,
                                                        paste0("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", Luna_ID, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d")), #convert to Date, then reformat YYYYMMDD
                                                        paste0("/gpfs/group/mnh5174/default/SPECC/MR_Proc/", tolower(SPECC_ID), "_", tolower(format((as.Date(ScanDate, format="%Y-%m-%d")), "%d%b%Y")))))

#verify that mr_dir is present as expected
idinfo$dirfound <- file.exists(idinfo$mr_dir)
subset(idinfo, dirfound==FALSE)
table(SPECC_rest[,c(3,5)])


##standard FD script
#SPECC_rest <- filter(SPECC_rest, pr_over5mm <= .15)
SPECC_rest <- filter(SPECC_rest, pr_over5mm <= .2)
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

################################################################################
#################create igraph object from adjacency matrix and label nodes V1, V2,...
allg <- apply(allmats, 1, function(sub) {
  g <- graph.adjacency(sub, mode="undirected", weighted=TRUE, diag=FALSE)
  V(g)$name <- paste0("V", 1:nrow(sub))
  g
})

for (subj in 1:length(allg)){
print(length(E(allg_noneg[[subj]])))
}

##remove negative correlations between nodes
allg_noneg <- lapply(allg, function(g) {
  delete.edges(g, which(E(g)$weight < 0))
})

################################################################################
#binarize and threshold graphs at densities ranging from 1-20%
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
##output will be 20 binary graphs AVERAGED ACROSS SUBJECTS (at 1-20% density) with community membership
#allmats_log <- apply(allmats, c(1,2,3), function(x) { log(x+.05)})
mean.graph <- apply(allmats, c(2,3), mean, na.rm = TRUE)
mean.g <- graph.adjacency(mean.graph, mode = "lower", weighted = TRUE, diag = FALSE)
mean.g <- delete.edges(mean.g, which(E(mean.g)$weight < 0))
mean.g.adj <- as_adjacency_matrix(mean.g, attr = "weight")

hist(E(mean.g)$weight)

E(mean.g)
dev.off()
if (use.infomap == 1){
  mean.g.infomap <- readRDS(paste0(basedir, "/cachedRfiles/infomap_communitylist.rds"))
} else{
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


##########condense bad communities into community 99: "badcomm"
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

data.frame(mean.g.infomap[[10]]$membership)
if(!file.exists(paste0(basedir,"/cachedRfiles/allg_density_infomap10.RData"))){save(allg_density, file = paste0(basedir,"/cachedRfiles/allg_density_infomap10.RData"))}
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

################################################################################
###IN PROGRESS: get rid of nodes with degree of 0
deg_0 <- array(NA, c(length(allmetrics), length(densities_desired), nnodes )) #3-D FALSE that's subjects x densities x nodes. 
for(s in 1:length(allmetrics)) {
  for(d in 1:length(densities_desired)){
    deg_0[s,d,] <- allmetrics[[s]][[d]]$degree == 0
  }
}
badnodes <- apply(deg_0, 2, function(mat) { colSums(mat) }) 

################################################################################
####Test group differences in network coherence and integration
#assign density of interest to dtarget and create cvec, which is numeric vector of community assignment for an arbitrary subject(1) at the desired density
#IN PROGRESS: expand across densities
dtarget <- 10
cvec <- allmetrics[[1]][[dtarget]]$community.membership

#WITHIN MODULE: first call loops over networks and subjects and outputs the mean within module degree z score for each subject and network at the predefined density
meanWithinCon <- sapply(sort(unique(cvec)), function(network) {
  sapply(allmetrics, function(subject) {
    mean(subject[[dtarget]]$within.module.deg.zscore[subject[[dtarget]]$community.membership == network])
  })
})

meanWithinCon <- cbind(meanWithinCon, SPECC_rest[, "BPD"], SPECC_rest[,"AgeAtScan"])
df.names <- c(sort(unique(cvec)),"BPD", "Age")
meanWithinCon <- data.frame(meanWithinCon)
colnames(meanWithinCon) <- df.names

within_mod_test <- list() 
for(n in 1:length(unique(cvec))){
  group.test <- t.test(meanWithinCon[,n] ~ meanWithinCon$BPD)
  age.test <- lm(meanWithinCon[,n] ~ meanWithinCon$BPD*meanWithinCon$Age)
  within_mod_test[[n]] <- list(group.test = group.test, age.test = summary(age.test))
}

#BETWEEN MODULE: first call loops over networks and subjects and outputs the mean within module degree z score for each subject and network at the predefined density
meanBetweenCon <- sapply(sort(unique(cvec)), function(network) {
  sapply(allmetrics, function(subject) {
    mean(subject[[dtarget]]$between.module.deg.zscore[subject[[dtarget]]$community.membership == network])
  })
})

meanBetweenCon <- cbind(meanBetweenCon, SPECC_rest[, "BPD"], SPECC_rest[,"AgeAtScan"])
df.names <- c(sort(unique(cvec)),"BPD", "Age")
meanBetweenCon <- data.frame(meanBetweenCon)
colnames(meanBetweenCon) <- df.names

between_mod_test <- list() 
for(n in 1:length(unique(cvec))){
  group.test <- t.test(meanBetweenCon[,n] ~ meanBetweenCon$BPD)
  age.test <- lm(meanBetweenCon[,n] ~ meanBetweenCon$BPD*meanBetweenCon$Age)
  between_mod_test[[n]] <- list(group.test = group.test, age.test = summary(age.test))
}

##plot of the age x BPD ixn in FPN between network connectivity
btw.mod.vals.fpn <- array(NA, length(meanBetweenCon[,1]), dimnames = list("BPD" = SPECC_rest[,"BPD"]))
for (i in 1:length(meanBetweenCon[,1])){
  btw.mod.vals.fpn[i] <- meanBetweenCon[i,4]
}
bw.plotobj <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], between_module_connectivity=btw.mod.vals.fpn)
bw.plotobj$BPD <- factor(bw.plotobj$BPD, levels=c(0,1), labels=c("Control", "BPD"))
# ss.metric.1$roiname <- "L Rolandic Operculum" 
ggplot(bw.plotobj, aes(x = AgeAtScan, color=BPD, y = between_module_connectivity)) + geom_point() + stat_smooth(method="lm", se=FALSE)+ 
  labs(x = "Age", y = "Average Between Network Connectivity", title = "Lateral Fronto-parietal Network") +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))


################################################################################
####Test group differences in global and nodal statistics
#Assign BPD and Age attributes to be used for group differences tests
attr(allmetrics, "bpd") <- SPECC_rest[,"BPD"]
attr(allmetrics, "age") <- SPECC_rest[,"AgeAtScan"]

attr(allmetrics.global, "bpd") <- SPECC_rest[,"BPD"]
attr(allmetrics.global, "age") <- SPECC_rest[,"AgeAtScan"]

################################################################################
###Testing global differences in shortest path length, and clustering coefficient by patient group, age, and group x age interactions
abst <- 2.38 #about p = .01
if(!exists("global.metrics")){ 
  global.metrics <- foreach (d = 1:length(densities_desired)) %do% {
	#for(d in 1:length(densities_desired)) {
    sigres <- 1
    results <- list()
    for (m in metricstorun.global) {
        this.metric.glob <- do.call(rbind, sapply(allmetrics.global, function(subj){subj[[d]][m]}))
        thisdf <- data.frame(group=factor(attr(allmetrics.global, "bpd"), levels=c(0,1), labels=c("control", "bpd")), 
                             age=attributes(allmetrics.global)$age, metric=this.metric.glob)
        colnames(thisdf) <- c("group", "age", "metric")

      group.test <- t.test(metric~group, thisdf, na.action = na.exclude)
      age.test <- lm(metric~group*age, thisdf, na.action = na.exclude)
      #remove snapshot of environment from lm object: https://blogs.oracle.com/R/entry/is_the_size_of_your
      #rm(list=ls(envir = attr(age.test$terms, ".Environment")), envir = attr(age.test$terms, ".Environment")) #THIS WILL REMOVE EVERYTHING IN GLOBAL ENVIRONMENT
			attr(age.test$terms, ".Environment") <- new.env() #add dummy environment
      pvec <- broom::tidy(age.test)$p.value[-1]
     
      if ((!is.null(group.test) && !is.nan(group.test$statistic) && abs(group.test$statistic) > abst) ||
          (!is.null(age.test) && !all(is.nan(pvec)) && any(pvec < .01))) {
        
        results[[sigres]] <- list(ttest=group.test, metric = m, density=densities_desired[d], agetest=age.test)
        sigres <- sigres+1
        } 
      }
    return(results)
  }
  try(stopCluster(clusterobj))
  
  
  save(global.metrics, file = paste0(basedir, "/cachedRfiles/global.metrics.binary.",pipeline,".RData"))
  
}

################################################################################
###Testing nodal differences by patient group, age, and group x age interactions
abst <- 2.64 ##about p < .005
if(!exists("node.metrics")){
  
  ##define empty vector of length 20 (# of densities) and an atlas based on the 269 dataset
  node.metrics <- vector(mode="list", length=length(densities_desired))
  names(node.metrics) <- densities_desired
  nodes <- 1:nnodes #for now
  metrics <- names(allmetrics[[1]][[1]]) #for now: metrics of first subject and first density
    
    node.metrics <- foreach(d=1:length(densities_desired), .packages="broom") %do% {
      sigres <- 1
      results <- list()
      
      #conduct two-sample t-tests for each node and each centrality measure and each density 
      #build a matrix of subjects and nodes
      
      for (m in metrics) {
        if(m == "eigen.cent" | m == "degree" | m == "closeness" | m == "betweenness.node" | m == "page.rank" | m == "part.coeff" | m == "within.module.deg.zscore" | m == "local.clustering" | 
           m == "between.module.deg.zscore" | m == "gateway.coeff.btw" | m == "gateway.coeff.deg" | m =="between.module.deg.zscore")
          
          for (n in nodes) {
            thismetric <- do.call(c, lapply(allmetrics, function(subj) {
              subj[[d]][[m]][n] #NB!!! If you inadvertently pass in 0.1, 0.2 etc., R rounds up to 1. Therefore, use d (index), not densities_desired[d].
            }))
            
            #combine these values with identifying variables -- age and BPD
            thisdf <- data.frame(group=factor(attr(allmetrics, "bpd"), levels=c(0,1), labels=c("control", "bpd")), 
                                 age=attributes(allmetrics)$age, metric=thismetric)
            
            
            #skip test if there are more than a few NaNs in metric (not defined)
            # if (sum(is.na(thisdf$metric)) > 5) {
            #   #should issue a message() here
            #   next
            # }
            
            node.test <- tryCatch(t.test(metric~group, thisdf, na.rm = TRUE), error = function(errorname) { print(errorname); return(NULL) })
            if (is.null(node.test)) { message("Error occurred for density: ", densities_desired[d], ", metric: ", m, " and node: ", n) }
            
            age.test <- tryCatch(lm(metric~age*group, thisdf, na.action = "na.exclude"), error=function(e) { print(e); return(NULL) })
            if (is.null(age.test)) { pvec <- NULL
            } else { pvec <- broom::tidy(age.test)$p.value[-1] } #p-values of age, bpd, and age x bpd. -1 to drop off intercept
            
            if ((!is.null(node.test) && !is.nan(node.test$statistic) && abs(node.test$statistic) > abst) ||
                (!is.null(age.test) && !all(is.nan(pvec)) && any(pvec < .01))) {
              
              results[[sigres]] <- list(nodename=as.character(atlas[n,5]), ttest=node.test, 
                                        metric = m, density=densities_desired[d],
                                        agetest=age.test, nodenum=n)
              sigres <- sigres+1
              
              
            }
          }
      }
      
      return(results)
    }
    
    try(stopCluster(clusterobj))
    
    save(node.metrics, file = paste0(basedir, "/cachedRfiles/node.metrics.binary.",pipeline,".RData"))
  }

##output yields a list with 20 elements (1 for each density) with significant results for a number of 
#centrality measures as well as participation coefficient and within-module degree z-score (within network connectivity)

##how many significant comparisons are there across densities
for(i in 1:length(node.metrics)){
  print(length(node.metrics[[i]]))
}
node.metrics[[10]]


#######compile significant findings
if(!exists("all.sig.nodal")){
    library(broom)
    
    allnodestats <- lapply(node.metrics, function(node.metrics_d) {
      x <- lapply(node.metrics_d, function(node) {
        df <- broom::tidy(node$ttest)
        df$nodename <- node$nodename
        df$nodenum <- node$nodenum
        df$density <- node$density
        df$metric <- node$metric
        return(df)
      })
      
      return(do.call(rbind, x))
    })
    
    alldf <- do.call(rbind, allnodestats)
    
    #alldf %>% group_by(nodename, metric) %>% print()
    
    #setwd("~/Box Sync/RS_BPD_graph/output.files")
    
    ##Main effect of BPD
    a <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="betweenness.node" & p.value < .005) 
    b <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="degree" & p.value < .005) #will have ns ps now that we blend age and group in output
    c <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="eigen.cent" & p.value < .005)
    d <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="local.clustering" & p.value < .005)
    e <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="closeness" & p.value < .005)
    f <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="page.rank" & p.value < .005)
    g <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="within.module.deg.zscore" & p.value < .005)
    h <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="part.coeff" & p.value < .005)
    y <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="gateway.coeff.btw" & p.value < .005)
    z <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="gateway.coeff.degree" & p.value < .005)
    aa <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="between.module.deg.zscore" & p.value < .005)
    
    
    allagestats <- lapply(node.metrics, function(node.metrics_d) {
      x <- lapply(node.metrics_d, function(node) {
        df <- broom::tidy(node$agetest)
        df$nodename <- node$nodename
        df$nodenum <- node$nodenum
        df$density <- node$density
        df$metric <- node$metric
        return(df)
      })
      
      return(do.call(rbind, x))
    })
    
    allagedf <- do.call(rbind, allagestats)
    
    #age main effects
    i <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="betweenness.node" & term=="age" & p.value < .01)
    j <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="degree" & term=="age" & p.value < .01) #will have ns ps now that we blend age and group in output
    k <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="eigen.cent" & term=="age" & p.value < .01)
    l <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="local.clustering" & term=="age" & p.value < .01)
    m <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="closeness" & term=="age" & p.value < .01) #will have ns ps now that we blend age and group in output
    n <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="page.rank" & term=="age" & p.value < .01)
    o <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="part.coeff" & term=="age" & p.value < .01)
    p <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="within.module.deg.zscore" & term=="age" & p.value < .01)
    bb <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="between.module.deg.zscore" & term=="age" & p.value < .01)
    cc <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="gateway.coeff.btw" & term=="age" & p.value < .01)
    dd <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="gateway.coeff.degree" & term=="age" & p.value < .01)
    
    
    
    #age x bpd interactions
    q <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="betweenness.node" & term=="age:groupbpd" & p.value < .01)
    r <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="eigen.cent" & term=="age:groupbpd" & p.value < .01)
    s <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="degree" & term=="age:groupbpd" & p.value < .01)
    t <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="local.clustering" & term=="age:groupbpd" & p.value < .01)
    u <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="closeness" & term=="age:groupbpd" & p.value < .01)
    v <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="page.rank" & term=="age:groupbpd" & p.value < .01)
    w <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="part.coeff" & term=="age:groupbpd" & p.value < .01)
    x <-allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="within.module.deg.zscore" & term=="age:groupbpd" & p.value < .01)
    ee <-allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="between.module.deg.zscore" & term=="age:groupbpd" & p.value < .01)
    ff <-allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="gateway.coeff.btw" & term=="age:groupbpd" & p.value < .01)
    gg <-allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="gateway.coeff.degree" & term=="age:groupbpd" & p.value < .01)
    
    library(gtools)
    all.sig.nodal <- bind_rows(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,bb,cc,dd,ee,ff,gg)
    write.csv(all.sig.nodal, file = paste0(basedir, "/output.files/all.sig.nodal.binary",pipeline,".csv"))
  }


#quick brute force edge analysis on fisher-transformed correlations
#allmats <- atanh(allmats) #fisher z transform
controls <- which(SPECC_rest$BPD == 0)
bpd <- which(SPECC_rest$BPD == 1)

edgepvals <- matrix(NA, nrow=nnodes, ncol=nnodes)
abst <- 3.2 #about p = .0001 (two-tailed)
edge.diffs <- data.frame()
for (i in 1:nnodes) {
  for (j in 1:nnodes) {
    if (i > j) {
      #lower triangle
      controlcorrs <- allmats[controls,i,j]
      bpdcorrs <- allmats[bpd, i, j]
      res <- t.test(atanh(controlcorrs), atanh(bpdcorrs))
      edgepvals[i,j] <- res$p.value
      edgepvals[j,i] <- res$statistic #t stat on upper triangle
      if (!is.nan(res$statistic) && abs(res$statistic) > abst) {
        cat("Significant edge difference: ", i, j, ", t=", round(res$statistic, 3), "p =", round(res$p.value, 3), "means are:", round(res$estimate, 5), "\n")
        sig.edge <- data.frame(edge1 = i, edge2 = j, tval = round(res$statistic, 3), p = round(res$p.value, 3), mean.control = round(res$estimate[1], 5), mean.bpd = round(res$estimate[2], 5))
        edge.diffs <- rbind(edge.diffs, sig.edge)
      }
    }
  }
}
row.names(edge.diffs) <- NULL
edge.diffs

####IN PROGRESSL: FDR correction for multiple comparisons
pvec <- edgepvals[lower.tri(edgepvals)]
#tvec <- edgepvals[upper.tri(edgepvals)]
lookup <- which(lower.tri(edgepvals) == TRUE, arr.ind=TRUE)
padj <- p.adjust(pvec, "fdr")
lookup[which(pvec < .0001),]
padj[which(padj < .05)]
lookup[which(padj < .05),]



###DELTACON analyses including overall network similarity and node and edge attributions 
if(!exists("deltacon_total") & !exists("node_stats_deltacon") & !exists("edge_diffs_deltacon")){
den10 <- lapply(allg_density, function(x) x[[10]])
attr(den10, "bpd") <- SPECC_rest[,"BPD"]
node_stats_deltacon <- c()
edge_diffs_deltacon <- c()
deltacon_total <- c()

for (i in 1:length(den10)) {
  for (j in 1:length(den10)) {
    if (i < j) {
      #only iterate over lower triangle of matrix
      #run delta con on i, j here and put into a N x N list: all_list[i,j] <- delta_con
      
      g1 <- den10[[i]]
      g2 <- den10[[j]]
      dc <- run_parse_deltacon(g1,g2)
      
      thisdc <- dc$deltacon
      thisdc <- data.frame(thisdc)
      thisdc$subj1 <- i
      thisdc$subj2 <- j
      thisdc$subj1.bpd <- attributes(den10)$bpd[i]
      thisdc$subj2.bpd <- attributes(den10)$bpd[j]
      deltacon_total <- rbind(deltacon_total, thisdc)

      ##nodal rankings
      thisnodecomp <- dc$nodestats
      thisnodecomp$subj1 <- i
      thisnodecomp$subj2 <- j
      thisnodecomp$subj1.bpd <- attributes(den10)$bpd[i]
      thisnodecomp$subj2.bpd <- attributes(den10)$bpd[j]
      node_stats_deltacon <- rbind(node_stats_deltacon, thisnodecomp)
      
      ##edge rankings
      thisedgecomp <- dc$edgestats
      thisedgecomp$subj1 <- i
      thisedgecomp$subj2 <- j
      thisedgecomp$subj1.bpd <- attributes(den10)$bpd[i]
      thisedgecomp$subj2.bpd <- attributes(den10)$bpd[j]
      edge_diffs_deltacon <- rbind(edge_diffs_deltacon, thisedgecomp)
      
    }
  }
}
saveRDS(edge_diffs_deltacon, paste0("edge_diffs_deltacon_", pipeline,".rds"))
saveRDS(deltacon_total, paste0("deltacon_total_", pipeline, ".rds"))
saveRDS(node_stats_deltacon, paste0("node_stats_deltacon_", pipeline, ".rds"))
}

####################################################
#####MH workbench deltacon 

# str(node_stats_deltacon)
# 
# nodesplit <- split(node_stats_deltacon, node_stats_deltacon$node)
# lapply(nodesplit, function(node) {
#   wibpd <- dplyr::filter(node, subj1.bpd==1 & subj2.bpd==1)
#   wicontrol <- dplyr::filter(node, subj1.bpd==0 & subj2.bpd==0)
#   
#   #compare whether this node is different in its within-group contribution
#   df <- data.frame(attr_stat=c(wibpd$attr_stat, wicontrol$attr_stat), group=c(rep("bpd", nrow(wibpd)), rep("control", nrow(wicontrol))))
#   t.test(attr_stat ~ group, df)
#   
#   bwgroups <- dplyr::filter(node, subj1.bpd != subj2.bpd)
#   browser()
# })
# library(ggplot2)
# ggplot(node_stats_deltacon, aes(x=factor(node), y=attr_stat)) + geom_boxplot()
# ggplot(node_stats_deltacon, aes(x=factor(node), y=attr_stat)) + stat_summary()
# 
# library(dplyr)
# node_agg <- node_stats_deltacon %>% group_by(node) %>% summarize(attr_m=mean(attr_stat), attr_se=sd(attr_stat)/sqrt(nrow(node_stats_deltacon)))
# 
# node_stats_deltacon$node_sort <- factor(node_stats_deltacon$node, levels=(node_agg$node[order(node_agg$attr_m, decreasing = TRUE)]))
# ggplot(node_stats_deltacon, aes(x=node_sort, y=attr_stat)) + stat_summary() + coord_flip()
# hist(node_agg$attr_m)
