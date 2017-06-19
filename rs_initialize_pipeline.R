########## RS_BPD_pipeline
####read in package dependencies and custom functions

setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/")
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
source("Infomap/infomap_communities.R")
source("functions/analyze_nodal_metrics_PCA_func.R")

#get_subj info here, includes motion exclusion procedure
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz")
table(subj_info[,c(4,8)])

#import raw adjacency matrices here (subj_info already contains the identified raw files)
allmats <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

#obtain weighted, non-negative weighted, and density-thresholded binary graphs
gobjs <- setup_graphs(allmats, allowCache=TRUE)

#gobjs contains a list of weighted, non-negative weighted, and binary matrices
#pull these out into single variables for simplicity
allg <- gobjs$allg; allg_noneg <- gobjs$allg_noneg; allg_density <- gobjs$allg_density

rm(gobjs) #remove from environment to save memory

#estimate and setup community structure. note: currently pulling from MH's community structure, unclear how this was generated
# comm_d10_l <- run_community_detection_on_agg(allmats, "louvain", density=0.1)
community <- readRDS(paste0(getwd(), "/cache/d10_louv_MH.rds"))

allg_noneg <- assign_communities(allg_noneg, community, "community")
allg_density <- assign_communities(allg_density, community, "community")

#compute global metrics on density-thresholded graphs
globalmetrics_dthresh <- compute_global_metrics(allg_density, allowCache=TRUE, community_attr="community") #community_attr determines how global/nodal statistics that include community are computed

#compute nodal metrics on density-thresholded graphs
nodalmetrics_dthresh <- compute_nodal_metrics(allg_density, allowCache=TRUE, community_attr="community") #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame

##calculate group comparisons (PCA pipeline)
sig_PCA_nodes <- analyze_nodal_metrics_PCA(nodalmetrics_dthresh$allmetrics.nodal.df, allowCache = FALSE)

##insert weighted and binary comparisons here



###STOPPED HERE: IN PROGRESS

