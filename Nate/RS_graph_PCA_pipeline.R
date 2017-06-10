##########RS_BPD_pipeline
####read in package dependencies and custom functions
#setwd("~/Box Sync/RS_BPD_graph")
setwd("/Users/michael/Data_Analysis/bpd_rest/Nate")
basedir <- getwd()

source("functions/setup_globals.R") #this will setup details of the parcellation, pipeline, and connection distance
source("functions/calcGraph_binary.R")
source("functions/import_adj_mats.R")
source("functions/get_subj_info.R")
source("functions/graph_util_redux.R")
source("functions/run_parse_deltacon.R")
source("functions/wibw_module_degree.R")

#get_subj info here
subj_info <- get_subj_info(adjmats_base, parcellation, pipeline, file_extension=".txt.gz")

#import raw adjacency matrices here (subj_info already contains the identified raw files)
allmats <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

#obtain weighted, non-negative weighted, and density-thresholded binary graphs
gobjs <- setup_graphs(allmats, allowCache=TRUE)

#gobjs contains a list of weighted, non-negative weighted, and binary matrices
#pull these out into single variables for simplicity
allg <- gobjs$allg
allg_noneg <- gobjs$allg_noneg
allg_density <- gobjs$allg_density

rm(gobjs) #remove from environment to save memory

#compute global metrics on density-thresholded graphs
globalmetrics_dthresh <- compute_global_metrics(allg_density, allowCache=TRUE)

###STOPPED HERE: IN PROGRESS




################################################################################
#######read in already processed Rdata files for faster run throughs:

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




#check num of edges across subjs
# for (subj in 1:length(allg)){
#   print(length(E(allg_noneg[[subj]])))
# }



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

mean.g.infomap

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

##################################################################################################
##############
####run tests on nodal metrics collapsed across densities
abst <- 2.64 ##about p < .005

metrics.toanalyze <- names(select(toanalyze, -subj, -node, -BPD, -Age))
sigres <- 1
results <- list()

for (m in metrics.toanalyze) {
    for (n in nodes) {
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
        
        results[[sigres]] <- list(nodename=as.character(atlas[n,5]), ttest=node.test, 
                                  metric = m, agetest=age.test, nodenum=n)
        sigres <- sigres+1
      }
    }
}



results.df <- do.call(rbind, results)

####BPD main effects 
results.ttest <- lapply(results, function(node){
  df <- broom::tidy(node$ttest)
  df$nodename <- node$nodename
  df$nodenum <- node$nodenum
  df$metric <- node$metric
  return(df)
})

results.ttest.df <- do.call(rbind, results.ttest)
a <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "central" & p.value < .005)
b <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "between.node" & p.value < .005)
c <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "within.mod" & p.value < .005)
d <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "between.mod" & p.value < .005)

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

