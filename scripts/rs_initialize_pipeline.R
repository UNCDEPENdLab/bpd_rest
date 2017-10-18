########## RS_BPD_pipeline
####read in package dependencies and custom functions
setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
#setwd("/Users/mnh5174/Data_Analysis/bpd_rest")

source("functions/setup_globals.R") #this will setup details of the parcellation, conn_method, preproc_pipeline, and connection distance, loads helper functions, and required packages

#get_subj info, includes motion scrubbing procedure
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".RData", fd.scrub = TRUE, allowCache = TRUE)
# subj_info
# table(subj_info[,c(4,8)])


##import raw adjacency matrices here (subj_info already contains the identified raw files)
allmats <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

###test: check density dists across densities for 3 random subjects
# pdf("dens.clime_hists_roix.pdf", width =12, height = 8)
# for (sub in sample(1:length(allmats),3)){
#   sub.mats <- allmats[[sub]]
#   for (den in seq(1, dim(sub.mats)[1], 1)){
#     plot.mat <- sub.mats[den,,]
#     plot.df <- data.frame(plot.mat)
#     plot.tidy <- gather(plot.df,  roi, value)
#     gg.obj <- ggplot(plot.tidy, aes(x = roi, y = value)) + geom_histogram(stat = "identity", position = "dodge") + labs(title = paste0("subject:", sub, ", density: ", den))
#     plot(gg.obj)
#   }
# }
# dev.off()

#obtain weighted, non-negative weighted, density-thresholded binary, and mean aggregate graphs
gobjs <- setup_graphs(allmats, allowCache=FALSE)

if(!conn_method == "dens.clime_partial"){
#gobjs contains a list of weighted, non-negative weighted, and binary matrices
#pull these out into single variables for simplicity
allg <- gobjs$allg; allg_noneg <- gobjs$allg_noneg; allg_density <- gobjs$allg_density; agg.g <- gobjs$agg.g
} else {
  #In the case of dens.clime
  allg_density <- gobjs[[1]] #list of weighted graphs. subjs X densities
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
} else {
  ####INSERT COMMUNITY DETECTION PROCEDURE HERE IF DESIRED
  community <- readRDS(paste0(getwd(), "/cache/power269/d12_louv_n83.rds"))
  
  allg_noneg <- assign_communities(allg_noneg, community, "community")
  allg_density <- assign_communities(allg_density, community, "community")
}

if(!conn_method == "dens.clime_partial"){
#compute global metrics on BINARY density-thresholded graphs
globalmetrics_dthresh <- compute_global_metrics(allg_density, allowCache=FALSE, community_attr="community") #community_attr determines how global/nodal statistics that include community are computed
#compute nodal metrics on BINARY density-thresholded graphs
nodalmetrics_dthresh <- compute_nodal_metrics(allg_density, allowCache=FALSE, community_attr="community") #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame
} else {
  globalmetrics_dthresh <- compute_global_metrics(allg_noneg, allowCache = FALSE, community_attr = "community")
  nodalmetrics_dthresh <- compute_nodal_metrics(allg_noneg, allowCache=FALSE, community_attr="community", weighted = TRUE) 
}
#########################PCA Analysis Pipeline. 
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

sig_PCA_nodes <- run_group_comparisons_nodal(toanalyze, abst = 2.64, allowCache = TRUE)

#again, pull sig_PCA_nodes apart for simplicity and remove sig_PCA data
signod.wle <- sig_PCA_nodes$all.sigwlelm.nodal.pca
signod.lm <- sig_PCA_nodes$all.siglm.nodal.pca
signod.bpd <- sig_PCA_nodes$all.sigttest.nodal.pca
rm(sig_PCA_nodes)

##########INSERT WEIGHTED AND BINARY COMPARISONS here if desired, these can be adapted easily from the PCA pipeline

######################################################################
###PLOT SIGNIFICANT FINDINGS
######################################################################
bpd.main.all <- plot_significant_groupeffects(signod.bpd) ##output from this function required for writing .node files for brainnet viewer

plot_significant_ageeffects(signod.lm)
plot_significant_ixn(signod.lm)

export_bpd_sig_nodefiles(bpd.main.all, node.file.dir)

##Post-hoc interrogation of nodes of interest: will export top percentage of significant 
#nodal difference for nodes of interest and exports to an outputdir for plotting in BNV
outputdir.interr <- paste0(basedir, "/BNV_nodefiles/schaefer422_ridge/interrogate_node_edgefiles/")
for (i in 1:nnodes){
  a <- interrogate_node(allmats, subj_info, i, t.stat = 2, outputdir = outputdir.interr)
}
##############################END NODAL FINDINGS PIPELINE##############################



##############################EDGE COMPARISONS##############################

##Conduct edge comparisons and export edge files for vizualization
edge.comp <- run_edge_comparisons(subj_info)

edge.outputdir <-  paste0(basedir, "/BNV_nodefiles/schaefer422_ridge/")
edge.bnv.output <- edge_bnv_files(edge.comp, edge.outputdir)
##

##IN PROGRESS: IMPLEMENTATION OF NBS AND DELTACON-ATTR

##this attempts to use brainGraph code, instead, let's rework for input into MATLAB
# allmats.nbs <- aperm(allmats, c(2,3,1))
# covars.NBS <- data.table(Group = factor(subj_info$BPD))
# con.vec <- c(1,-1)
# nbs.out <- NBS(allmats.df, covars.NBS, con.vec, p.init = .001, N = 1000, symmetric = TRUE)


#reexport files for NBS
subj_info <- arrange(subj_info, BPD, SPECC_ID)
subj_info$NBS_ID <- 1:nrow(subj_info)
subj_info <- as.data.frame(subj_info) #get rid of tibble weirdness

for (i in 1:nrow(subj_info)) {
  x <- read.table(subj_info[i,"file"])
  write.table(file=paste0("/Users/nth7/Documents/MATLAB/NBS1.2/schaefer422_aroma_ridge.net_partial/matrices_numeric/subj_", sprintf("%03d", i), ".txt"), x=x, row.names=FALSE, col.names=FALSE)
}

nbs.mat.ixn <- model.matrix(lm(NBS_ID ~ BPD*AgeAtScan, subj_info))
subj_info$BPD
write.table(file="/Users/nth7/Documents/MATLAB/NBS1.2/schaefer422_aroma_ridge.net_partial/design_matrix_bpdxage.txt", x=nbs.mat.ixn, row.names=FALSE, col.names=FALSE)
model.matrix(lm(NBS_ID~BPD, subj_info))

design_matrix_nbs <- subj_info$NBS_ID
design_matrix_nbs <- cbind(design_matrix_nbs, ifelse(design_matrix_nbs == 0, 1, 0))
write.table(design_matrix_nbs, file ="/Users/nth7/Documents/MATLAB/NBS1.2/schaefer422_aroma_ridge.net_partial/design_matrix_group.txt", row.names = FALSE, col.names = FALSE)

