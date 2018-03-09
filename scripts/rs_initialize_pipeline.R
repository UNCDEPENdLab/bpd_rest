##################### initialize_pipeline #####################
##note: this script initializes graph analysis up through calculation of global and nodal metrics. 
#Before running this script, the user is required to load libraries and source custom functions (setup_globals) and specify the pipeline inputs

if (!exists("basedir")){setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()}

# Euclidean Distance estimation -------------------------------------------

#estimate ROI-ROI Euclidean distance (use cache if available since the double for loop takes ~4 seconds)
expectFile <- file.path(basedir, "cache", paste0("distmat_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
if (file.exists(expectFile)) {
  load(expectFile)
} else {
  roiDist <- matrix(NA, nrow=nnodes, ncol=nnodes) 
  
  system.time(for (i in 1:nnodes) {
    for (j in 1:nnodes) {
      #populate lower triangle only for speed, then transpose onto upper triangle
      if (i > j) { roiDist[i,j] <- sqrt((roiMat[i,1] - roiMat[j,1])^2 + (roiMat[i,2] - roiMat[j,2])^2 + (roiMat[i,3] - roiMat[j,3])^2) }
    }
  })
  
  #copy lower triangle to upper triangle
  roiDist[upper.tri(roiDist)] <- t(roiDist)[upper.tri(roiDist)]
  diag(roiDist) <- 0
  
  ####quick QA
  # hist(roiDist)
  # vecDist <- roiDist[lower.tri(roiDist)]
  # sum(vecDist < 20)/length(vecDist)
  
  #roi.dist can be changed at the front end of the script
  rmShort <- roiDist > roi.dist
  #creates binary matrix, in which 0 denotes a short distanced connection that is to be removed. 
  rmShort <- apply(rmShort, c(1,2), as.numeric)
  save(file=expectFile, rmShort, roiDist)
}


# Subject Info and import adjmats -----------------------------------------

#get_subj info, includes motion scrubbing procedure
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)
#table(subj_info$BPD, subj_info$Female)


##import raw adjacency matrices here (subj_info already contains the identified raw files)
allmats <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

# Setup Graphs and Assign Community Structure -----------------------------


#obtain weighted, non-negative weighted, proportional and fc density-thresholded binary, and mean aggregate graphs
gobjs <- setup_graphs(allmats, file_tag = file_tag, file_tag_nothresh = file_tag_nothresh, fc_out_rm = fc_out_rm, allowCache=TRUE)


if(!conn_method == "dens.clime_partial"){
  #gobjs contains a list of weighted, non-negative weighted, and binary matrices
  #pull these out into single variables for simplicity
  allg <- gobjs$allg; allg_noneg <- gobjs$allg_noneg; allg_density <- gobjs$allg_density; agg.g <- gobjs$agg.g; allg_density_fc <- gobjs$allg_density_fc
} else {
  #In the case of dens.clime
  allg_density <- gobjs[[1]] #list of weighted graphs. Based on proportional thresholding
  allg_noneg <- gobjs[[2]] # list of weighted non-negative graphs. subjs X densities
  agg.g <- gobjs[[3]] # single aggregate mean graph from ridge 
}
rm(gobjs) #remove from environment to save memory

#create heatmaps on aggregated data and random subjects (typically, this can stay commented out)
#generate_heatmaps(agg.g, graph_list = allg_noneg)

#community assignment
if(use.yeo == 1){ 
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

if(!conn_method == "dens.clime_partial"){
  if(thresh == "fc"){
    #compute global metrics on BINARY density-thresholded graphs
    globalmetrics_dthresh <- compute_global_metrics(allg_density_fc, allowCache=FALSE, community_attr="community") #community_attr determines how global/nodal statistics that include community are computed
    #compute nodal metrics on BINARY density-thresholded graphs
    nodalmetrics_dthresh <- compute_nodal_metrics(allg_density_fc, allowCache=FALSE, community_attr="community") #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame
    
    nodalmetrics_dthresh$allmetrics.nodal.df$density <- rep(rep(rs_desired_log, each = 422),length(allg)) 
    allmetrics.nodal.df <- nodalmetrics_dthresh$allmetrics.nodal.df
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
    toanalyze_thresh <- reduce_centrality_fa(allmetrics.nodal.df, reducemetrics = reducemetrics, den = 0, allowCache = FALSE, browse = TRUE)
    
    toanalyze_thresh <- left_join(toanalyze_thresh, membership_df, by = "node")
    toanalyze_thresh$id <- as.character(toanalyze_thresh$id)
    
    fa.metrics <- colnames(dplyr::select(toanalyze_thresh, -id, -node, -BPD, -Age, -membership))  
  } else{
    toanalyze_thresh <- reduce_centrality_fa(allmetrics.nodal.df, reducemetrics = reducemetrics, allowCache = TRUE, weighted = FALSE, browse = FALSE, fc = FALSE)
    
    toanalyze_thresh <- left_join(toanalyze_thresh, membership_df, by = "node")
    toanalyze_thresh$id <- as.character(toanalyze_thresh$id)
    
    fa.metrics <- colnames(dplyr::select(toanalyze_thresh, -id, -node, -BPD, -Age, -membership))  
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



# Compute Weighted Graph Metrics (nodal) ------------------------------------------

if(include.weighted == 1) {
  if(!conn_method == "dens.clime_partial") {
    #compute global metrics on weighted graphs
    globalmetrics_weighted <- compute_global_metrics(allg_noneg, allowCache=TRUE, community_attr="community", weighted = 1) #community_attr determines how global/nodal statistics that include community are computed
    #compute nodal metrics on weighted graphs
    nodalmetrics_weighted <- compute_nodal_metrics(allg_noneg, allowCache=TRUE, community_attr="community", weighted = 1) #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame
    
    nodalmetrics_weighted_df <- nodalmetrics_weighted$allmetrics.nodal.df #pull data frame for analysis
    nodalmetrics_weighted_df <- subj_info %>% dplyr::select(BPD, Age = AgeAtScan, id = SPECC_ID) %>% inner_join(nodalmetrics_weighted_df, by = "id") #add BPD and age columns
    
    ##data reduction. in progress
    #toanalyze_weighted <- reduce_centrality_fa(nodalmetrics_weighted_df, reducemetrics = reducemetrics, allowCache = TRUE, weighted = TRUE, browse = TRUE)
  }
}

