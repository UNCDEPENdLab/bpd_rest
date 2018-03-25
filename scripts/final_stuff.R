
# Final Pipeline ----------------------------------------------------------

setwd("~/Box Sync/bpd_rest"); basedir <- getwd()
#initialize the graph analysis pipeline (includes sourcing pipeline inputs, creating graphs, thresholding, binarization, community assignment, and calculation and reduction of global and nodal graph metrics)

source("scripts/setup_globals.R")

#specify pipeline inputs.leave blank for defaults, consult with function for deviations from default#RIDGE
inputs <- specify_inputs(thresh_weighted = "binary", 
                         conn_method = "cor.shrink",
                         fc_out_rm = FALSE, 
                         preproc_pipeline = "nosmooth_aroma_bp_nonaggr",
                         reducemetrics =  c("degree", "betweenness.node", "part.coeff", "eigen.cent", "gateway.coeff.btw", "gateway.coeff.degree", "within.module.deg", "leverage.cent", "page.rank", "local.clustering", "closeness"),
                         rs_desired_log = logspace(log10(.5), log10(.65), 15)) 

for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

##if you want to add on an additional tag, make sure to check that this tag is relfected in dir structure
add_tag <- "update2018"
file_tag <- paste0(file_tag, "_", add_tag)
file_tag_nothresh <- paste0(file_tag_nothresh, "_", add_tag)

source("scripts/estimate_euclidean_distance.R") ##creates rmShort which will delete edges close in euclidean distance


# Subject Info and import adjmats -----------------------------------------

#get_subj info, includes motion scrubbing procedure. 003BU and 008JH have had their data truncated to 300 volumes
subj_info <- get_subj_info(adjmats_base, 
                           parcellation, 
                           conn_method, 
                           preproc_pipeline, 
                           add_tag = add_tag,
                           file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)

##import raw adjacency matrices here (subj_info already contains the identified raw files)
##for ridge remove the short euclidean distance removal
allmats <- import_adj_mats(subj_info, rmShort = NULL, allowCache=TRUE)

# hist(allmats_)
# Setup Graphs and Assign Community Structure -----------------------------

#obtain weighted, non-negative weighted, proportional and fc density-thresholded binary, and mean aggregate graphs
gobjs <- setup_graphs(allmats, file_tag = file_tag, file_tag_nothresh = file_tag_nothresh, fc_out_rm = fc_out_rm, allowCache=TRUE)

if(!conn_method == "dens.clime_partial"){
  #gobjs contains a list of weighted, non-negative weighted, and binary matrices
  #pull these out into single variables for simplicity
  allg <- gobjs$allg; allg_noneg <- gobjs$allg_noneg; allg_density <- gobjs$allg_density; agg.g <- gobjs$agg.g; allg_density_fc <- gobjs$allg_density_fc
  agg.g.controls <- gobjs$agg.g.controls; agg.g.bpd <- gobjs$agg.g.bpd
} else {
  #In the case of dens.clime
  allg_density <- gobjs[[1]] #list of weighted graphs. Based on proportional thresholding
  allg_noneg <- gobjs[[2]] # list of weighted non-negative graphs. subjs X densities
  agg.g <- gobjs[[3]] # single aggregate mean graph from ridge 
}

rm(gobjs) #remove from environment to save memory


#community assignment
if (use.yeo == 1) { 
  yeo7 <- yeo7_community(agg.g)
  
  allg_noneg <- assign_communities(allg_noneg, yeo7, "community")
  allg_density <- assign_communities(allg_density, yeo7, "community")
  allg_density_fc <- assign_communities(allg_density_fc, yeo7, "community")
  
  membership_df <- data.frame(node = rownames(do.call(cbind,yeo7)), membership = as.numeric(do.call(cbind,yeo7)[,1]))
  
  community.names <- data.frame(node = rownames(do.call(cbind,yeo7)), membership = as.numeric(do.call(cbind,yeo7)[,1]),
                                community.name = mapvalues(membership_df$membership, from = c("1","2","3","4","5","6","7"), to = c("VIS", "SOMMOT", "DORSATTN", "SALVENTATTN", "LIMBIC", "FPN", "DMN")))
} else {
  ####INSERT COMMUNITY DETECTION PROCEDURE HERE IF DESIRED
  community <- readRDS(paste0(getwd(), "/cache/power269/d12_louv_n83.rds"))
  allg_noneg <- assign_communities(allg_noneg, community, "community")
  allg_density <- assign_communities(allg_density, community, "community")
}

# Compute Thresholded Graph Metrics (nodal)  --------------------------------------------

if (!conn_method == "dens.clime_partial") {
  if(thresh == "fc"){
    #compute global metrics on BINARY density-thresholded graphs
    globalmetrics_dthresh <- compute_global_metrics(allg_density_fc, allowCache=TRUE, community_attr="community") #community_attr determines how global/nodal statistics that include community are computed
    #compute nodal metrics on BINARY density-thresholded graphs
    nodalmetrics_dthresh <- compute_nodal_metrics(allg_density_fc, allowCache=FALSE, community_attr="community") #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame
    
    #nodalmetrics_dthresh <- get(load("/Users/nth7/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/cache/threshnodalmetrics_transformed_blahschaefer422_nosmooth_aroma_bp_nonaggr_ridge.net_partial_fc_binary_all.RData"))
    
    #nodalmetrics_dthresh$allmetrics.nodal.df$density <- rep(rep(rs_desired_log, each = 422),length(allg)) #should be superseded by adding this attribute to each graph in compute_nodal_metrics
    allmetrics.nodal.df <- nodalmetrics_dthresh$allmetrics.nodal.df
    
    globalmetrics_dthresh.df <- globalmetrics_dthresh$allmetrics.global.df
    # globalmetrics_dthresh.df$density <- rep(rs_desired_log, length(allg))

    #what mean densities did we hit with FC thresholding?
    #globalmetrics_dthresh.df %>% group_by(wthresh) %>% dplyr::summarize(mean(edge_density), median(edge_density))
  } else {
    globalmetrics_dthresh <- compute_global_metrics(allg_density, allowCache=TRUE, community_attr="community") #community_attr determines how global/nodal statistics that include community are computed
    #compute nodal metrics on BINARY density-thresholded graphs
    nodalmetrics_dthresh <- compute_nodal_metrics(allg_density, allowCache=TRUE, community_attr="community") #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame
    
    nodalmetrics_dthresh$allmetrics.nodal.df$density <- rep(rep(densities_desired, each = 422),length(allg)) 
    allmetrics.nodal.df <- nodalmetrics_dthresh$allmetrics.nodal.df
  }
} else {
  globalmetrics_dthresh <- compute_global_metrics(allg_density, allowCache = TRUE, community_attr = "community")
  nodalmetrics_dthresh <- compute_nodal_metrics(allg_density, allowCache=TRUE, community_attr="community", weighted = FALSE) 
}

#couple of sanity checks
# qplot(allmetrics.nodal.df$degree)
# ggplot(allmetrics.nodal.df, aes(x=degree)) + geom_histogram() + facet_wrap(~wthresh)

# Transform the nodal graph metrics ---------------------------------------
source(file.path(basedir, "scripts", "transform_graph_metrics.R"))
# source(paste0(basedir, "/scripts/transform_graph_metrics_10rs.R"))

# Data reduce Binary(FA/PCA) ----------------------------------------------------

if (data.reduce == "fa"){
  if(thresh == "fc"){
    ##nodal
    toanalyze_thresh <- reduce_centrality_fa(allmetrics.nodal.df, 
                                             reducemetrics = reducemetrics, 
                                             den = 0, 
                                             allowCache =FALSE, 
                                             browse = TRUE)
    
    #toanalyze_thresh <- reduce_centrality_fa(allmetrics.nodal.df, reducemetrics = reducemetrics, den = 0, allowCache = TRUE, browse = FALSE)
    for(i in 1:length(toanalyze_thresh)) assign(names(toanalyze_thresh)[i], toanalyze_thresh[[i]])
    
    toanalyze <- left_join(toanalyze, membership_df, by = "node")
    toanalyze$id <- as.character(toanalyze$id)
    
    fa.metrics <- colnames(dplyr::select(toanalyze, -id, -node, -BPD, -Age, -membership))  
    
    heatmap_FA(faout)
    
   
  } else{
    toanalyze <- reduce_centrality_fa(allmetrics.nodal.df, reducemetrics = reducemetrics, allowCache = TRUE, weighted = FALSE, browse = FALSE, fc = FALSE)
    
    toanalyze <- left_join(toanalyze, membership_df, by = "node")
    toanalyze$id <- as.character(toanalyze$id)
    
    fa.metrics <- colnames(dplyr::select(toanalyze, -id, -node, -BPD, -Age, -membership))  
  }
} else {
  #run PCA across metrics and densities and pull scores into toanalyze
  #NOTE PCA a_priori and all have only been verified to work with pearson and ridge
  if(PCA.method == "a_priori"){
    pca.out <- PCA_a_priori(nodalmetrics_dthresh$allmetrics.nodal.df, pcametrics, allowCache = TRUE, den = .05)
    ##Note: this will depend on which conn_method you select above
    toanalyze <- pca.out[[1]]
    PClabs <- pca.out[[2]]
    metrics.pca <- pca.out[[3]]
  }
  if(PCA.method == "all"){
    toanalyze <- PCA_all(nodalmetrics_dthresh$allmetrics.nodal.df, pcametrics, allowCache = TRUE, den = .05)
  }
}



# Run binary network-level MLM on reduced data and plot ----------------------------
mlm_res <- run_group_mlm(toanalyze, allowCache = TRUE, browse = FALSE)#, thresh = thresh)
# mlm_plot <- plot_group_mlm(mlm_res, browse = TRUE)


# Run binary nodal comparisons on reduced data and plot -------------------
sig_nodes <- run_group_comparisons_nodal(dplyr::select(toanalyze, -membership), abst = 2.64, allowCache = FALSE, browse = TRUE)

bpd.main.all <- plot_significant_groupeffects(sig_nodes$all.sigttest.nodal, toanalyze)

source("functions/plot_significant_ixn.R")

signod.lm <- sig_nodes$all.siglm.nodal
plot_significant_ixn(signod.lm, toanalyze = toanalyze)

###If you want to plot individual nodes by network to look at which nodes may be contributing to age x bpd ixns:
# plot_nodes_by_network(toanalyze)
#beepr::beep("mario")

###Interrogate nodes of interest
# outputdir.interr <- paste0(basedir, "/BNV_nodefiles/interrogate_node_edgefiles/")
# a <- interrogate_node(allmats, subj_info, 422, t.stat = 2, outputdir = outputdir.interr)


trans.file
