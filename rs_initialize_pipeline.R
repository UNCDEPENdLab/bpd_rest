########## RS_BPD_pipeline
####read in package dependencies and custom functions
#setwd("~/Box Sync/RS_BPD_graph")
setwd("/Users/michael/Data_Analysis/bpd_rest")
basedir <- getwd()

#this will setup details of the parcellation, conn_method, preproc_pipeline, and connection distance
#it also sources all helper functions for additional analysis
source("functions/setup_globals.R") 

#get_subj info here
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub=TRUE, allowCache=TRUE)

#import raw adjacency matrices here (subj_info already contains the identified raw files)
allmats <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

#obtain weighted, non-negative weighted, and density-thresholded binary graphs
gobjs <- setup_graphs(allmats, allowCache=TRUE)

#gobjs contains a list of weighted, non-negative weighted, and binary matrices
#pull these out into single variables for simplicity
allg <- gobjs$allg; allg_noneg <- gobjs$allg_noneg; allg_density <- gobjs$allg_density

rm(gobjs) #remove from environment to save memory

#estimate and setup community structure (most work on exploring this has moved to determine_communities.R)
# comm_weighted_louvain <- run_community_detection_on_agg(allmats, "louvain")
# comm_weighted_greedy <- run_community_detection_on_agg(allmats, "fast_greedy")
# comm_d15 <- run_community_detection_on_agg(allmats, "louvain", density=0.15)
# comm_d20 <- run_community_detection_on_agg(allmats, "louvain", density=0.20)
# compare(comm_weighted_louvain, comm_d15, method="nmi")
# compare(comm_weighted_louvain, comm_d20, method="nmi")
# compare(comm_weighted_louvain, comm_weighted_greedy, method="nmi")
# comm_infomap_d10 <- run_community_detection_on_agg(allmats, "infomap", density=0.10)
# comm_infomap_weighted <- run_community_detection_on_agg(allmats, "infomap", hierarchical=FALSE, verbose=FALSE)
community <- readRDS(paste0(getwd(), "/cache/d12_louv_n83.rds"))

allg_noneg <- assign_communities(allg_noneg, community, "community")
allg_density <- assign_communities(allg_density, community, "community")

#assign weighted louvain into weighted and density-thresholded structures in attribute "wcomm_louvain"
#allg_noneg <- assign_communities(allg_noneg, comm_weighted_louvain, "wcomm_louvain")
#allg_density <- assign_communities(allg_density, comm_weighted_louvain, "wcomm_louvain")

#compute global metrics on density-thresholded graphs
globalmetrics_dthresh <- compute_global_metrics(allg_density, allowCache=TRUE, community_attr="community") #community_attr determines how global/nodal statistics that include community are computed

#compute nodal metrics on density-thresholded graphs
nodalmetrics_dthresh <- compute_nodal_metrics(allg_density, allowCache=TRUE, community_attr="community") #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame

##calculate group comparisons (PCA pipeline)
sig_PCA_nodes <- analyze_nodal_metrics_PCA(nodalmetrics_dthresh$allmetrics.nodal.df, allowCache = TRUE)


###STOPPED HERE: IN PROGRESS




################################################################################
#######read in already processed Rdata files for faster run throughs:



#significant nodal comparisons
#if(file.exists(paste0(basedir, "/cachedRfiles/node.metrics.binary.",pipeline,".RData")) == TRUE) {
#  node.metrics <- get(load(paste0(basedir, "/cachedRfiles/node.metrics.binary.",pipeline,".RData")))
#}
##compiled significant nodal comparisons
#if(file.exists(paste0(basedir, "/output.files/all.sig.nodal.",pipeline,".rds")) == TRUE) {
#  all.sig.nodal <- read.csv(paste0(basedir, "/output.files/all.sig.nodal.",pipeline,".csv"))
#}
##total deltacon stats
#if(file.exists(paste0(basedir, "/output.files/deltacon_total_", pipeline, ".rds")) == TRUE) {
#  deltacon_total <- readRDS(paste0(basedir, "/output.files/deltacon_total_", pipeline, ".rds"))
#}
##edge attribution deltacon values
#if(file.exists(paste0(basedir, "/output.files/edge_diffs_deltacon_", pipeline, ".rds")) == TRUE) {
#  edge_diffs_deltacon <- readRDS(paste0(basedir, "/output.files/edge_diffs_deltacon_", pipeline,".rds"))
#}
##nodal attribution deltacon values
#if(file.exists(paste0(basedir, "/output.files/node_stats_deltacon_", pipeline, ".rds")) == TRUE) {
#  node_stats_deltacon <- readRDS(paste0(basedir, "/output.files/node_stats_deltacon_", pipeline, ".rds"))
#}

#check num of edges across subjs
# for (subj in 1:length(allg)){
#   print(length(E(allg_noneg[[subj]])))
# }



#mean.g.infomap

##########condense bad communities into the last community: community-dependent results from the last module should not be interpreted
#V(mean.g.community[[10]])$community ##use this if using louvain

#df <- as.data.frame(table(mean.g.infomap[[10]]$membership))
#badcomm <- df$Var1[df$Freq < 4]
#goodcomm <- as.numeric(df$Var1[df$Freq >= 4])
#mean.g.infomap[[10]]$membership[mean.g.infomap[[10]]$membership %in% badcomm] <- length(goodcomm) + 1   #where highest community number is extraneous nodes
#mean.g$community <- mean.g.infomap[[10]]$membership
#table(mean.g.infomap[[10]]$membership)
## plot(mean.g.infomap[[10]], mean.g)
#
##Assign community from mean graph at a given density (10) back into subjects
##each element of allg_density is a list of subjects with 20 binary graphs for each subject at 1-20% density, now with communities assigned to them
#allg_density <- lapply(allg_density, function(subj) {
#  subj_graphs <- lapply(1:length(densities_desired), function(d) {
#    V(subj[[d]])$community <- mean.g.infomap[[10]]$membership
#    return(subj[[d]])
#  })
#})
#
#
#mean.g.info.df <- data.frame(mean.g.infomap[[10]]$membership)
#colnames(mean.g.info.df) <- NULL

#if(!file.exists(paste0(basedir,"/cachedRfiles/allg_density_infomap10.RData"))){save(allg_density, file = paste0(basedir,"/cachedRfiles/allg_density_infomap10.RData"))}

#example: what is community assignment for subject [[1]] at density [[20]]..will be the same across densities
#V(allg_density[[1]][[20]])$community

################################################################################

################################################################################




