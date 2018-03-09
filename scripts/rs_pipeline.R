##################### RS_BPD_pipeline #####################
setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
#initialize the graph analysis pipeline (includes sourcing pipeline inputs, creating graphs, thresholding, binarization, community assignment, and calculation and reduction of global and nodal graph metrics)

source("scripts/setup_globals.R")
inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE) #leave blank for defaults, consult with function for deviations from default
for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

source("scripts/rs_initialize_pipeline.R")

# Run binary network-level MLM on reduced data and plot ----------------------------
mlm_res <- run_group_mlm(toanalyze_thresh, allowCache = TRUE, browse = FALSE, thresh = thresh)
mlm_plot <- plot_group_mlm(mlm_res, browse = TRUE)


# Run binary nodal comparisons on reduced data and plot  --------------------------------------------------
metrics.toanalyze <- names(dplyr::select(toanalyze_binary, -id, -node, -BPD, -Age, -membership))
ttest_bpd <- list()

for (metric in metrics.toanalyze) {
  for (node in levels(toanalyze_binary$node)) {
    thismetric <- toanalyze_binary[toanalyze_binary$node == node, c(metric, "BPD", "Age")]
    thismetric$BPD <- factor(thismetric$BPD, levels = c(0,1), labels = c("control", "BPD"))
    colnames(thismetric) <- c("metric", "group", "age")
    
    node.test <- tryCatch(t.test(metric~group, thismetric, na.rm = TRUE), error = function(errorname) { print(errorname); return(NULL) })
    if (is.null(node.test)) { message("Error occurred for metric: ", metric, " and node: ", node) }
    ttest_bpd[[metric]][[length(ttest_bpd)+1]] <- data.frame(broom::tidy(node.test), nodename=as.character(atlas$anat_label[atlas$name == node]), metric=metric, nodenum=node)
    
  }
}

sig_nodes_reduced <- run_group_comparisons_nodal(dplyr::select(toanalyze_binary, -membership), abst = 2.64, allowCache = TRUE, browse = FALSE)
bpd.main.all <- plot_significant_groupeffects(sig_nodes_reduced$all.sigttest.nodal, toanalyze = toanalyze_binary, fa.metrics)
plot_significant_ixn(signod.wle)


# Run weighted group MLM --------------------------------------------------
nodalmetrics_weighted_df_reduced <- nodalmetrics_weighted_df %>% dplyr::select(BPD, Age, id, node, eigen.cent, strength, betweenness.node, within.module.deg, part.coeff,community.membership) %>% dplyr::rename(membership = community.membership)

mlm_res <- run_group_mlm(nodalmetrics_weighted_df_reduced, allowCache = FALSE, browse = FALSE, weighted = TRUE)
plot_group_mlm(mlm_res, browse = FALSE); beep(sound = 8)


# Run weighted nodal comparisons and plot ---------------------------------

#conduct group statistical tests
sig_nodes_weighted <- run_group_comparisons_nodal(nodalmetrics_weighted_df, abst = 2.64, allowCache = TRUE, weighted = TRUE)

#again, pull significant results apart and remove the object
signod.wle_weighted <- sig_nodes_weighted$all.sigwlelm.nodal
signod.lm_weighted <- sig_nodes_weighted$all.siglm.nodal
signod.bpd_weighted <- sig_nodes_weighted$all.sigttest.nodal

weighted.metrics <- unique(signod.bpd_weighted$metric)

bpd.main.all.wieghted <- plot_significant_groupeffects(signod.bpd_weighted, nodalmetrics_weighted_df, weighted.metrics, weighted = TRUE)
plot_significant_ixn(signod.wle_weighted, nodalmetrics_weighted_df, reducemetrics = weighted.metrics, weighted = TRUE)

# export_bpd_sig_nodefiles(bpd.main.all, node.file.dir)
# Interrogate nodes of interest -------------------------------------------
##Post-hoc interrogation of nodes of interest: will export top percentage of significant 
#nodal difference for nodes of interest and exports to an outputdir for plotting in BNV
outputdir.interr <- paste0(basedir, "/BNV_nodefiles/schaefer422_ridge/interrogate_node_edgefiles/")
for (i in 1:nnodes){
  a <- interrogate_node(allmats, subj_info, i, t.stat = 2, outputdir = outputdir.interr)
}
# Run Edge Comparisons ----------------------------------------------------
##Conduct edge comparisons and export edge files for vizualization
edge.comp <- run_edge_comparisons(subj_info)

edge.outputdir <-  paste0(basedir, "/BNV_nodefiles/schaefer422_ridge/")
edge.bnv.output <- edge_bnv_files(edge.comp, edge.outputdir)
##


