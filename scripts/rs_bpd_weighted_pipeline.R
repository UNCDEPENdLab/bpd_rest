########## RS_BPD_weighted_pipeline
####read in package dependencies and custom functions
setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
source("functions/setup_globals.R") #this will setup details of the parcellation, conn_method, preproc_pipeline, and connection distance, loads helper functions, and required packages

#get_subj info, includes motion scrubbing procedure
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".RData", fd.scrub = TRUE, allowCache = TRUE)

##import raw adjacency matrices here (subj_info already contains the identified raw files)
allmats <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

#obtain weighted, non-negative weighted, density-thresholded binary, and mean aggregate graphs
gobjs <- setup_graphs(allmats, allowCache=TRUE)
allg <- gobjs$allg
allg_noneg <- gobjs$allg_noneg
agg.g <- gobjs$agg.g

rm(gobjs) #remove from environment to save memory

#community assignment
if(use.yeo == 1){ 
  yeo7 <- yeo7_community(agg.g)
  allg_noneg <- assign_communities(allg_noneg, yeo7, "community")
  allg <- assign_communities(allg, yeo7, "community")
} else {
  warning("COMMUNITY DETECTION NOT CURRENTLY SUPPORTED")
  ####INSERT COMMUNITY DETECTION PROCEDURE HERE IF DESIRED
  community_louv <- readRDS(paste0(getwd(), "/cache/power269/d12_louv_n83.rds"))
  allg_noneg <- assign_communities(allg_noneg, community, "community")
  allg <- assign_communities(allg, community, "community")
}
#compute global metrics on weighted graphs
globalmetrics_weighted <- compute_global_metrics(allg_noneg, allowCache=TRUE, community_attr="community", weighted = TRUE) #community_attr determines how global/nodal statistics that include community are computed
#compute nodal metrics on weighted graphs
nodalmetrics_weighted <- compute_nodal_metrics(allg_noneg, allowCache=TRUE, community_attr="community", weighted = TRUE) #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame

#pull data frame for analysis
nodalmetrics_weighted_df <- nodalmetrics_weighted$allmetrics.nodal.df

#add BPD and age columns
nodalmetrics_weighted_df <- subj_info %>% dplyr::select(BPD, Age = AgeAtScan, id = SPECC_ID) %>% inner_join(nodalmetrics_weighted_df, by = "id")

#conduct group statistical tests
sig_nodes_weighted <- run_group_comparisons_nodal(nodalmetrics_weighted_df, abst = 2.64, allowCache = TRUE, weighted = TRUE)


#again, pull significant results apart and remove the object
signod.wle <- sig_nodes_weighted$all.sigwlelm.nodal
signod.lm <- sig_nodes_weighted$all.siglm.nodal.pca
signod.bpd <- sig_nodes_weighted$all.sigttest.nodal
rm(sig_nodes_weighted)


# Data reduce (FA/PCA) ----------------------------------------------------



######################################################################
###PLOT SIGNIFICANT FINDINGS
######################################################################
bpd.main.all <- plot_significant_groupeffects(signod.bpd) ##output from this function required for writing .node files for brainnet viewer

plot_significant_ageeffects(signod.lm)
plot_significant_ixn(signod.lm)

export_bpd_sig_nodefiles(bpd.main.all, node.file.dir)


##########INSERT WEIGHTED AND BINARY COMPARISONS here if desired, these can be adapted easily from the PCA pipeline
sig_weighted_nodes <- run_group_comparisons_nodal()


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

