## Load libraries and custom functions

# Stop if basedir doesn’t exist -------------------------------------------


if(!exists("basedir")) {
  stop("You must define the basedir for the project")
} else if (!file.exists(basedir)) {
  stop("basedir: ", basedir, " does not exist.")
} 
# Load necessary libraries ------------------------------------------------


suppressMessages(library(igraph))
suppressMessages(library(brainGraph)) ##current version has been adapted to not allow GUI compatibility
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(gdata))
suppressMessages(library(beepr))
suppressMessages(library(parallel))
suppressMessages(library(data.table))
suppressMessages(library(psych))
suppressMessages(library(tidyr))
suppressMessages(library(reshape2))
suppressMessages(library(broom))
suppressMessages(library(gtools))
suppressMessages(library(Hmisc))
suppressMessages(library(wle))
suppressMessages(library(lsmeans))
suppressMessages(library(pracma))
suppressMessages(library(sjPlot))
if(!file.exists("/gpfs/group/mnh5174/default")){suppressMessages(library(erer))}

suppressMessages(library(DescTools))
suppressMessages(library(taRifx))
suppressMessages(library(scatr))
suppressMessages(library(doMC))
suppressMessages(library(plyr))
suppressMessages(library(gridExtra))
suppressMessages(library(permute))



# from BG manual, not sure if this will come up later ---------------------


OS <- .Platform$OS.type

num.cores <- detectCores()
registerDoMC(num.cores)




# Source Custom Functions -------------------------------------------------

source("functions/get_subj_info.R")
source("functions/calcGraph_nodal.R")
source("functions/calcGraph_global.R")
source("functions/compute_global_metrics.R")
source("functions/compute_nodal_metrics.R")
source("functions/import_adj_mats.R")
source("functions/setup_graphs.R")
source("functions/setup_community.R")
source("functions/graph_util_redux.R")
source("functions/run_parse_deltacon.R")
source("functions/wibw_module_degree.R")
source("functions/analyze_nodal_metrics_PCA_func.R")
source("functions/edge_functions.R")
source("Infomap/run_infomap.R")
source("functions/interrogate_node.R")
source("functions/generate_heatmaps.R")
source("functions/PCA_a_priori.R")
source("functions/PCA_all.R")
source("functions/run_group_comparisons_nodal.R")
source("functions/plot_sig_results_nodal.R")
source("functions/plot_significant_ageeffects.R")
source("functions/plot_significant_groupeffects.R")
#source("functions/plot_significant_ixn.R") ##issue, since this calls on subj_info it must be read in later
source("functions/export_bpd_sig_nodefiles.R")
source("functions/generate_agg.g.R")
source("functions/gateway_coeff_NH.R")
source("functions/run_group_mlm.R")
source("functions/threshold_graph.R")
source("functions/reduce_centrality_fa.R")
source("functions/plot_group_mlm.R")
source("functions/lmerCellMeans.R")
source("functions/table_glht.R")
source("functions/reduce_networks_for_sem.R")
source("functions/get_meanfc_outliers.R")
source("functions/specify_inputs.R")
source("functions/heatmap_fa.R")
source("functions/reduce_centrality_global.R")
source("functions/plot_nodes_by_network.R")
source("functions/BoxCox_extract.R")
source("functions/vertex_roles.R")
source("functions/centr_lev.R")
source("functions/fa.CFI.R")
source("functions/format_to_bg.R")
source("functions/mtpc_NH.R") #copied from github repo for debugging
source("functions/brainGraph_GLM_NH.R") #copied from github repo for debugging