########## RS_BPD_pipeline
####read in package dependencies and custom functions
setwd("~/Box Sync/DEPENd/Projects/bpd_rest/")
#setwd("/Users/mnh5174/Data_Analysis/bpd_rest")
basedir <- getwd()

source("functions/setup_globals.R") #this will setup details of the parcellation, conn_method, preproc_pipeline, and connection distance
source("functions/get_subj_info.R")
source("functions/calcGraph_binary.R")
source("functions/import_adj_mats.R")
source("functions/setup_graphs.R")
source("functions/setup_community.R")
source("functions/graph_util_redux.R")
source("functions/run_parse_deltacon.R")
source("functions/wibw_module_degree.R")

#get_subj info here, includes motion exclusion procedure
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz")

#import raw adjacency matrices here (subj_info already contains the identified raw files)
allmats <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

#obtain weighted, non-negative weighted, and density-thresholded binary graphs
gobjs <- setup_graphs(allmats, allowCache=TRUE)

#gobjs contains a list of weighted, non-negative weighted, and binary matrices
#pull these out into single variables for simplicity
allg <- gobjs$allg; allg_noneg <- gobjs$allg_noneg; allg_density <- gobjs$allg_density

rm(gobjs) #remove from environment to save memory

#estimate and setup community structure
comm_weighted_louvain <- run_community_detection_on_agg(allmats, "louvain")
comm_weighted_greedy <- run_community_detection_on_agg(allmats, "fast_greedy")
# comm_weighted_infomap <- run_community_detection_on_agg(allmats, "infomap")

comm_d05 <- run_community_detection_on_agg(allmats, "louvain", density=0.05)
comm_d10 <- run_community_detection_on_agg(allmats, "louvain", density=0.10)
comm_d15 <- run_community_detection_on_agg(allmats, "louvain", density = 0.15)
comm_d20 <- run_community_detection_on_agg(allmats, "louvain", density=0.20)
compare(comm_weighted_louvain, comm_d05, method="nmi")
compare(comm_weighted_louvain, comm_d10, method="nmi")
compare(comm_weighted_louvain, comm_d15, method = "nmi")
compare(comm_weighted_louvain, comm_d20, method = "nmi")
compare(comm_weighted_louvain, comm_weighted_greedy)
compare(comm_d20, comm_d10, method = "nmi")
compare(comm_d05, comm_weighted_greedy, method="nmi")

#assign weighted louvain into weighted and density-thresholded structures in attribute "wcomm_louvain"
allg_noneg <- assign_communities(allg_noneg, comm_weighted_louvain, "comm_weighted_louvain")
#allg_density <- assign_communities(allg_density, comm_weighted_louvain, "wcomm_louvain")

allg_noneg_comm1 <- get.vertex.attribute(allg_noneg$`001RA`)$name[which(get.vertex.attribute(allg_noneg$`001RA`)$comm_weighted_louvain == 1)]
allg_noneg_comm2 <- as.numeric(gsub("V", "", get.vertex.attribute(allg_noneg$`001RA`)$name))[which(get.vertex.attribute(allg_noneg$`001RA`)$comm_weighted_louvain == 2)]
allg_noneg_comm3 <- as.numeric(gsub("V", "", get.vertex.attribute(allg_noneg$`001RA`)$name))[which(get.vertex.attribute(allg_noneg$`001RA`)$comm_weighted_louvain == 3)]
allg_noneg_comm4 <- as.numeric(gsub("V", "", get.vertex.attribute(allg_noneg$`001RA`)$name))[which(get.vertex.attribute(allg_noneg$`001RA`)$comm_weighted_louvain == 4)]
allg_noneg_comm5 <- as.numeric(gsub("V", "", get.vertex.attribute(allg_noneg$`001RA`)$name))[which(get.vertex.attribute(allg_noneg$`001RA`)$comm_weighted_louvain == 5)]


node.file <- NodeFile(atlas = atlas, 
                      # community = get.vertex.attribute(allg_noneg$`001RA`)$comm_weighted_louvain == 1,
                      nnodes = nnodes,
                      nodestp = allg_noneg_comm1,
                      labels = 0,
                      filename = paste0(parcellation, "_", preproc_pipeline, "_", conn_method, "_comm_weighted_louvain_1"),
                      outputdir = paste0(getwd(), "/BNV_nodefiles/weighted_ridge_comms")
)

#compute global metrics on density-thresholded graphs
globalmetrics_dthresh <- compute_global_metrics(allg_density, allowCache=TRUE)

#compute nodal metrics on density-thresholded graphs
nodalmetrics_dthresh <- compute_nodal_metrics(allg_density, allowCache=TRUE) #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame





###STOPPED HERE: IN PROGRESS

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




