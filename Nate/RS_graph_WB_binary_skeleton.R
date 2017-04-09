##########RS_BPD_pipeline
####read in package dependencies and custom functions
setwd("~/Box Sync/RS_BPD_graph")

source("calcGraph_binary.R")
source("Graph_util_redux.R")
source("run_parse_deltacon.R")
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
nnodes <- 269  #varying this will require you to change the .txt coordinate file you read in, should be 4 X nnodes txt file labelled:"x", "y", "z", "roi"
roiMat <- read.table("~/Box Sync/RS_BPD_graph/bb264coordinate_appended_shift_nate_culled.txt", header=FALSE, col.names=c("x", "y", "z", "roi"))
use.infomap <- 1
use.scotmi <- 0 #currently, SCOTMI still using 265 X 265
if(use.scotmi == 0){pipeline <- "pearson"} else{pipeline  <- "scotmi"}
roi.dist <- 20

################################################################################
#######read in already processed Rdata files for faster run throughs:
#computed global metrics
if(file.exists(paste0("allmetrics.global.", pipeline, ".FDremove.rds")) == TRUE) {
  allmetrics.global <- readRDS(paste0("allmetrics.global.", pipeline, ".FDremove.rds"))
} 
#computed nodal metrics
if(file.exists(paste0("allmetrics.", pipeline, ".FDremove.rds")) == TRUE) {
  allmetrics  <- readRDS(paste0("allmetrics.", pipeline, ".FDremove.rds"))
} 
#significant nodal comparisons
if(file.exists(paste0("node.metrics.",pipeline,".fd.rds")) == TRUE) {
  node.metrics <- readRDS(paste0("node.metrics.",pipeline,".fd.rds"))
}
#compiled significant nodal comparisons
if(file.exists(paste0("all.sig.nodal.",pipeline,".rds")) == TRUE) {
  all.sig.nodal <- read.csv(paste0("all.sig.nodal.",pipeline,".csv"))
}
#total deltacon stats
if(file.exists(paste0("deltacon_total_", pipeline, ".rds")) == TRUE) {
  deltacon_total <- readRDS(paste0("deltacon_total_", pipeline, ".rds"))
}
#edge attribution deltacon values
if(file.exists(paste0("edge_diffs_deltacon_", pipeline, ".rds")) == TRUE) {
  edge_diffs_deltacon <- readRDS(paste0("edge_diffs_deltacon_", pipeline,".rds"))
}
#nodal attribution deltacon values
if(file.exists(paste0("node_stats_deltacon_", pipeline, ".rds")) == TRUE) {
  node_stats_deltacon <- readRDS(paste0("node_stats_deltacon_", pipeline, ".rds"))
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
#####filter subjects with over .15 brain volumes displaces .5mm or more
table(SPECC_rest[,c(3,5)])
SPECC_rest <- filter(SPECC_rest, pr_over5mm <= .15)
table(SPECC_rest[,c(3,5)])

describe(SPECC_rest[,c(1:6, 8)])

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

##remove negative correlations between nodes
allg_noneg <- lapply(allg, function(g) {
  delete.edges(g, which(E(g)$weight < 0))
})

################################################################################
#binarize and threshold graphs at densities ranging from 1-20%
densities_desired <- seq(.01, .2, .01) 

library(foreach)
library(doSNOW)

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


if (use.infomap == 1){
  mean.g.infomap <- readRDS("infomap_communitylist.rds")
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
#example: what is community assignment for subject [[1]] at density [[20]]..will be the same across densities
#V(allg_density[[1]][[20]])$community

################################################################################
##compute global metrics 
if(!exists("allmetrics.global")) {
    allmetrics.global <- foreach(subj=allg_density, .packages = c("igraph", "brainGraph")) %do% {
      #for (subj in allg_density) { #put here for more fine-grained debugging
      #require(igraph)
      dl <- lapply(subj, function(dgraph) {
        calcGraph_binary_global(dgraph)
      })
      names(dl) <- paste0("d", densities_desired)
      dl
    }
    saveRDS(allmetrics.global, paste0("allmetrics.global.", pipeline,".FDremove.rds"))
  } 

################################################################################
###compute nodal metrics 
  if(!exists("allmetrics")) {
    allmetrics <- foreach(subj=allg_density, .packages = c("igraph", "brainGraph")) %dopar% {
      #for (subj in allg_density) { #put here for more fine-grained debugging
      #require(igraph)
      dl <- lapply(subj, function(dgraph) {
        
        calcGraph_binary_nodal(dgraph)
      })
      names(dl) <- paste0("d", densities_desired)
      dl
    }
    
    saveRDS(allmetrics, paste0("allmetrics.", pipeline, ".FDremove.rds"))
  } 

################################################################################
###get rid of nodes with degree of 0
deg_0 <- array(NA, c(length(allmetrics), length(densities_desired), nnodes )) #3-D FALSE that's subjects x densities x nodes. 
for(s in 1:length(allmetrics)) {
  for(d in 1:length(densities_desired)){
    deg_0[s,d,] <- allmetrics[[s]][[d]]$degree == 0
  }
}
badnodes <- apply(deg_0, 2, function(mat) { colSums(mat) }) 

################################################################################
####Test group differences
#Assign BPD and Age attributes to be used for group differences tests
attr(allmetrics, "bpd") <- SPECC_rest[,"BPD"]
attr(allmetrics, "age") <- SPECC_rest[,"AgeAtScan"]

##define empty vector of length 20 (# of densities) and an atlas based on the 269 dataset
node.metrics <- vector(mode="list", length=length(densities_desired))
names(node.metrics) <- densities_desired

atlas <- read.csv("power269_masterlookup.csv", header = TRUE)

nodes <- 1:nnodes #for now
metrics <- names(allmetrics[[1]][[1]]) #for now: metrics of first subject and first density
abst <- 3.0 #about p = .005

if(!exists(node.metrics)){
    metricstorun <- c("eigen.cent","degree", "closeness", "betweenness.node", "page.rank",  "part.coeff", "within.module.deg.zscore", "local.clustering")
    
    registerDoSNOW(clusterobj)
    node.metrics <- foreach(d=1:length(densities_desired), .packages="broom") %dopar% {
      sigres <- 1
      results <- list()
      
      #conduct two-sample t-tests for each node and each centrality measure and each density 
      #build a matrix of subjects and nodes
      
      for (m in metrics) {
        # browser()
        if(m == "eigen.cent" | m == "degree" | m == "closeness" | m == "betweenness.node" | m == "page.rank" | m == "part.coeff" | m == "within.module.deg.zscore" | m == "local.clustering")
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
    saveRDS(node.metrics, paste0("node.metrics.",pipeline,".fd.rds"))
  }
##output yields a list with 20 elements (1 for each density) with significant results for a number of 
#centrality measures as well as participation coefficient and within-module degree z-score (within network connectivity)

##how many significant comparisons are there across densities
for(i in 1:length(node.metrics)){
  print(length(node.metrics[[i]]))
}
node.metrics[[10]][[46]]


#######compile significant findings
if(!exists(all.sig.nodal)){
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
    setwd("~/Box Sync/RS_BPD_graph/output.node")
    
    
    ##Main effect of BPD
    a <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="betweenness.node" & p.value < .005) 
    b <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="degree" & p.value < .005) #will have ns ps now that we blend age and group in output
    c <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="eigen.cent" & p.value < .005)
    d <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="local.clustering" & p.value < .005)
    e <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="closeness" & p.value < .005)
    f <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="page.rank" & p.value < .005)
    g <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="within.module.deg.zscore" & p.value < .005)
    h <- alldf %>% arrange(nodename, metric, density) %>% filter(metric=="part.coeff" & p.value < .005)
    
    
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
    
    #age x bpd interactions
    q <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="betweenness.node" & term=="age:groupbpd" & p.value < .01)
    r <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="eigen.cent" & term=="age:groupbpd" & p.value < .01)
    s <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="degree" & term=="age:groupbpd" & p.value < .01)
    t <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="local.clustering" & term=="age:groupbpd" & p.value < .01)
    u <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="closeness" & term=="age:groupbpd" & p.value < .01)
    v <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="page.rank" & term=="age:groupbpd" & p.value < .01)
    w <- allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="part.coeff" & term=="age:groupbpd" & p.value < .01)
    x <-allagedf %>% arrange(nodename, metric, density) %>% filter(metric=="within.module.deg.zscore" & term=="age:groupbpd" & p.value < .01)
    
    library(gtools)
    all.sig.nodal <- bind_rows(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x)
    write.csv(all.sig.nodal, file = past0("all.sig.nodal.",pipeline,".csv"))
  }


#quick brute force edge analysis on fisher-transformed correlations
#allmats <- atanh(allmats) #fisher z transform
controls <- which(SPECC_rest$BPD == 0)
bpd <- which(SPECC_rest$BPD == 1)

edgepvals <- matrix(NA, nrow=nnodes, ncol=nnodes)
abst <- 3.2 #about p = .0001 (two-tailed)
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
      
      }
    }
  }
}

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



str(node_stats_deltacon)

nodesplit <- split(node_stats_deltacon, node_stats_deltacon$node)
lapply(nodesplit, function(node) {
  wibpd <- dplyr::filter(node, subj1.bpd==1 & subj2.bpd==1)
  wicontrol <- dplyr::filter(node, subj1.bpd==0 & subj2.bpd==0)
  
  #compare whether this node is different in its within-group contribution
  df <- data.frame(attr_stat=c(wibpd$attr_stat, wicontrol$attr_stat), group=c(rep("bpd", nrow(wibpd)), rep("control", nrow(wicontrol))))
  t.test(attr_stat ~ group, df)
  
  bwgroups <- dplyr::filter(node, subj1.bpd != subj2.bpd)
  browser()
})
library(ggplot2)
ggplot(node_stats_deltacon, aes(x=factor(node), y=attr_stat)) + geom_boxplot()
ggplot(node_stats_deltacon, aes(x=factor(node), y=attr_stat)) + stat_summary()

library(dplyr)
node_agg <- node_stats_deltacon %>% group_by(node) %>% summarize(attr_m=mean(attr_stat), attr_se=sd(attr_stat)/sqrt(nrow(node_stats_deltacon)))

node_stats_deltacon$node_sort <- factor(node_stats_deltacon$node, levels=(node_agg$node[order(node_agg$attr_m, decreasing = TRUE)]))
ggplot(node_stats_deltacon, aes(x=node_sort, y=attr_stat)) + stat_summary() + coord_flip()
hist(node_agg$attr_m)
