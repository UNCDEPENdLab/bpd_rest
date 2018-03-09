
# Final Pipeline NEGATIVE EDGES----------------------------------------------------------

setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
#initialize the graph analysis pipeline (includes sourcing pipeline inputs, creating graphs, thresholding, binarization, community assignment, and calculation and reduction of global and nodal graph metrics)

source("scripts/setup_globals.R")
#RIDGE
inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, preproc_pipeline = "nosmooth_aroma_bp_nonaggr", neg_edges = TRUE) #leave blank for defaults, consult with function for deviations from default
#PEARSON
#inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, conn_method = "pearson", rs_desired_log = logspace(log10(.2), log10(.4), 20)) #leave blank for defaults, consult with function for deviations from default
for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

source("scripts/estimate_euclidean_distance.R") ##creates rmShort which will delete edges close in euclidean distance


# Subject Info and import adjmats -----------------------------------------

#get_subj info, includes motion scrubbing procedure. 003BU and 008JH have had their data truncated to 300 volumes
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)
#table(subj_info$BPD, subj_info$Female)


##import raw adjacency matrices here (subj_info already contains the identified raw files)
allmats <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

# for (i in 1:83) {
#   thissubj <- allmats[i,,]
#   hist(thissubj)
# }

# Setup Graphs and Assign Community Structure -----------------------------

#obtain weighted, non-negative weighted, proportional and fc density-thresholded binary, and mean aggregate graphs
gobjs <- setup_graphs(allmats, file_tag = file_tag, file_tag_nothresh = file_tag_nothresh, fc_out_rm = fc_out_rm, allowCache=TRUE)

if(!conn_method == "dens.clime_partial"){
  #gobjs contains a list of weighted, non-negative weighted, and binary matrices
  #pull these out into single variables for simplicity
  allg <- gobjs$allg; allg_noneg <- gobjs$allg_noneg; allg_density <- gobjs$allg_density
  agg.g <- gobjs$agg.g; allg_density_fc <- gobjs$allg_density_fc; agg.g.controls <- gobjs$agg.g.controls
  agg.g.bpd <- gobjs$agg.g.bpd; allg_nopos <- gobjs$allg_nopos; allg_density_fc_neg <- gobjs$allg_density_fc_neg
} else {
  #In the case of dens.clime
  allg_density <- gobjs[[1]] #list of weighted graphs. Based on proportional thresholding
  allg_noneg <- gobjs[[2]] # list of weighted non-negative graphs. subjs X densities
  agg.g <- gobjs[[3]] # single aggregate mean graph from ridge 
}
rm(gobjs) #remove from environment to save memory


#community assignment
if(use.yeo == 1){ 
  yeo7 <- yeo7_community(agg.g)
  
  allg_noneg <- assign_communities(allg_noneg, yeo7, "community")
  allg_density <- assign_communities(allg_density, yeo7, "community")
  allg_density_fc <- assign_communities(allg_density_fc, yeo7, "community")
  allg_density_fc_neg <- assign_communities(allg_density_fc_neg, yeo7, "community")
  allg_nopos <- assign_communities(allg_nopos, yeo7, "community")
  
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

if(!conn_method == "dens.clime_partial"){
  if(thresh == "fc"){
    #compute global metrics on BINARY density-thresholded graphs
    globalmetrics_dthresh <- compute_global_metrics(allg_density_fc_neg, allowCache=FALSE, community_attr="community") #community_attr determines how global/nodal statistics that include community are computed
    #compute nodal metrics on BINARY density-thresholded graphs
    nodalmetrics_dthresh <- compute_nodal_metrics(allg_density_fc_neg, allowCache=FALSE, community_attr="community") #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame
    
    nodalmetrics_dthresh$allmetrics.nodal.df$density <- rep(rep(rs_desired_log, each = 422),length(allg)) 
    allmetrics.nodal.df <- nodalmetrics_dthresh$allmetrics.nodal.df
    
    globalmetrics_dthresh.df <- globalmetrics_dthresh$allmetrics.global.df
    globalmetrics_dthresh.df$density <- rep(rs_desired_log, length(allg))
  } else{
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


# Data reduce Binary(FA/PCA) ----------------------------------------------------

if (data.reduce == "fa"){
  if(thresh == "fc"){
    ##nodal
    toanalyze_thresh <- reduce_centrality_fa(allmetrics.nodal.df, reducemetrics = reducemetrics, den = 0, allowCache = FALSE, browse = FALSE)
    for(i in 1:length(toanalyze_thresh)) assign(names(toanalyze_thresh)[i], toanalyze_thresh[[i]])
    
    toanalyze <- left_join(toanalyze, membership_df, by = "node")
    toanalyze$id <- as.character(toanalyze$id)
    
    fa.metrics <- colnames(dplyr::select(toanalyze, -id, -node, -BPD, -Age, -membership))  
    
    qplot(toanalyze$integration)
    
    ###the extracted factor scores are pretty skewed.. try and make them as normal as possible
    toanalyze$integration <- log(2- Winsorize(toanalyze$integration, minval = quantile(toanalyze$integration, .01), maxval = quantile(toanalyze$integration, .99)))
    toanalyze$within.mod <- Winsorize(toanalyze$within.mod, minval = quantile(toanalyze$within.mod, .01), maxval = quantile(toanalyze$within.mod, .99))
    toanalyze$central <- Winsorize(toanalyze$central, minval = quantile(toanalyze$central, .01), maxval = quantile(toanalyze$central, .99))
    
    heatmap_FA(faout)
    
    ##global 
    ##IN PROGRESS
    # toanalyze_global <- reduce_centrality_global(globalmetrics_dthresh.df)
    
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
mlm_res <- run_group_mlm(toanalyze, allowCache = FALSE, browse = FALSE)#, thresh = thresh)
mlm_plot <- plot_group_mlm(mlm_res, browse = FALSE)



# Run binary nodal comparisons on reduced data and plot -------------------
sig_nodes <- run_group_comparisons_nodal(dplyr::select(toanalyze, -membership), abst = 2.64, allowCache = TRUE, browse = FALSE)

bpd.main.all <- plot_significant_groupeffects(sig_nodes$all.sigttest.nodal, toanalyze)

source("functions/plot_significant_ixn.R")

signod.lm <- sig_nodes$all.siglm.nodal
plot_significant_ixn(signod.lm, toanalyze = toanalyze)

###If you want to plot individual nodes by network to look at which nodes may be contributing to age x bpd ixns:
plot_nodes_by_network(toanalyze)
beepr::beep("mario")

###Interrogate nodes of interest
outputdir.interr <- paste0(basedir, "/BNV_nodefiles/interrogate_node_edgefiles/")
a <- interrogate_node(allmats, subj_info, 422, t.stat = 2, outputdir = outputdir.interr)


# self-report data --------------------------------------------------------
selfreports <- read.csv("data/allselfreports.csv")
colnames(selfreports)[1] <- "id"
write_csv(selfreports, "data/allselfreports.csv")

dmn_toanalyze <- toanalyze %>% dplyr::filter(membership == 7) %>% left_join(selfreports, by = "id")

# lm(dmn_toanalyze$)
