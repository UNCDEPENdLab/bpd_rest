
# Final Pipeline comparison----------------------------------------------------------

setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
#initialize the graph analysis pipeline (includes sourcing pipeline inputs, creating graphs, thresholding, binarization, community assignment, and calculation and reduction of global and nodal graph metrics)

source("scripts/setup_globals.R")
#RIDGE
inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, preproc_pipeline = "nosmooth_aroma_bp") #leave blank for defaults, consult with function for deviations from default
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
  allg <- gobjs$allg; allg_noneg <- gobjs$allg_noneg; allg_density <- gobjs$allg_density; agg.g <- gobjs$agg.g; allg_density_fc <- gobjs$allg_density_fc; agg.g.controls <- gobjs$agg.g.controls; agg.g.bpd <- gobjs$agg.g.bpd
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
    globalmetrics_dthresh <- compute_global_metrics(allg_density_fc, allowCache=TRUE, community_attr="community") #community_attr determines how global/nodal statistics that include community are computed
    #compute nodal metrics on BINARY density-thresholded graphs
    nodalmetrics_dthresh <- compute_nodal_metrics(allg_density_fc, allowCache=TRUE, community_attr="community") #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame
    
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


##### separate for comparison
allmetrics_1 <- allmetrics.nodal.df


# SECOND PIPELINE ---------------------------------------------------------


inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, preproc_pipeline = "nosmooth_aroma_bp_nonaggr") #leave blank for defaults, consult with function for deviations from default

for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

source("scripts/estimate_euclidean_distance.R") ##creates rmShort which will delete edges close in euclidean distance


# Subject Info and import adjmats -----------------------------------------

#get_subj info, includes motion scrubbing procedure. 003BU and 008JH have had their data truncated to 300 volumes
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)
#table(subj_info$BPD, subj_info$Female)


##import raw adjacency matrices here (subj_info already contains the identified raw files)
allmats_2 <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

# for (i in 1:83) {
#   thissubj <- allmats[i,,]
#   hist(thissubj)
# }

# Setup Graphs and Assign Community Structure -----------------------------

#obtain weighted, non-negative weighted, proportional and fc density-thresholded binary, and mean aggregate graphs
gobjs_2 <- setup_graphs(allmats_2, file_tag = file_tag, file_tag_nothresh = file_tag_nothresh, fc_out_rm = fc_out_rm, allowCache=FALSE)

if(!conn_method == "dens.clime_partial"){
  #gobjs contains a list of weighted, non-negative weighted, and binary matrices
  #pull these out into single variables for simplicity
  allg_2 <- gobjs_2$allg; allg_noneg_2 <- gobjs_2$allg_noneg; allg_density_2 <- gobjs_2$allg_density; agg.g_2 <- gobjs_2$agg.g; allg_density_fc_2 <- gobjs_2$allg_density_fc; agg.g.controls_2 <- gobjs_2$agg.g.controls; agg.g.bpd_2 <- gobjs_2$agg.g.bpd
} else {
  #In the case of dens.clime
  allg_density <- gobjs_2[[1]] #list of weighted graphs. Based on proportional thresholding
  allg_noneg <- gobjs_2[[2]] # list of weighted non-negative graphs. subjs X densities
  agg.g <- gobjs_2[[3]] # single aggregate mean graph from ridge 
}
rm(gobjs_2) #remove from environment to save memory


#community assignment
if(use.yeo == 1){ 
  yeo7 <- yeo7_community(agg.g)
  
  allg_noneg_2 <- assign_communities(allg_noneg_2, yeo7, "community")
  allg_density_2 <- assign_communities(allg_density_2, yeo7, "community")
  allg_density_fc_2 <- assign_communities(allg_density_fc_2, yeo7, "community")
  
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
    globalmetrics_dthresh_2 <- compute_global_metrics(allg_density_fc_2, allowCache=TRUE, community_attr="community") #community_attr determines how global/nodal statistics that include community are computed
    #compute nodal metrics on BINARY density-thresholded graphs
    nodalmetrics_dthresh_2 <- compute_nodal_metrics(allg_density_fc_2, allowCache=TRUE, community_attr="community") #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame
    
    nodalmetrics_dthresh_2$allmetrics.nodal.df$density <- rep(rep(rs_desired_log, each = 422),length(allg)) 
    allmetrics.nodal.df_2 <- nodalmetrics_dthresh_2$allmetrics.nodal.df
    
    globalmetrics_dthresh.df_2 <- globalmetrics_dthresh_2$allmetrics.global.df
    globalmetrics_dthresh.df_2$density <- rep(rs_desired_log, length(allg))
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


# Now, compare ------------------------------------------------------------

##compare degree within a given FC thresh value
allmetrics_2 <- allmetrics.nodal.df_2

#degree_1 <- allmetrics_1 %>% dplyr::select(id, node, density, degree)

#to avoid filtering weirdness
allmetrics_1$density <- as.character(allmetrics_1$density)
allmetrics_2$density <- as.character(allmetrics_2$density)

results <- c()
for(i in subj_info$SPECC_ID){
  this.subj_1 <- allmetrics_1 %>% dplyr::filter(density == "0.0144024653753876") %>% 
    dplyr::filter(id == i)
  
  this.subj_2 <- allmetrics_2 %>% dplyr::filter(density == "0.0144024653753876") %>% 
    dplyr::filter(id == i)
  
  this.compare <- data.frame(s1 = this.subj_1$degree, s2 = this.subj_2$degree)
  
  res <- cor.test( ~ s1 + s2, 
            data=this.compare,
            method = "spearman",
            continuity = FALSE,
            conf.level = 0.95)
  
  results <- c(results, unname(res$estimate))
}

qplot(results, binwidth = .03)
str(allmetrics_1)
unique(allmetrics_1$density)

##compare mean graphs across pipe 1 and 2
mean_1 <- apply(allmats, c(2,3), mean)
mean_2 <- apply(allmats_2, c(2,3), mean)

cormat <- cor(mean_1, mean_2)  

meltcor <- melt(cormat)
ggplot(data = meltcor, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() 

hist(cor(mean_1, mean_2))

##FC comparisons across subjects
x <- array(NA, 83)
for (i in 1:83){
  subj_1 <- allmats[i,,]
  subj_2 <- allmats_2[i,,]
  
  x[i] <- cor(subj_1[lower.tri(subj_1)], subj_2[lower.tri(subj_2)])
}
x
mean(x)
