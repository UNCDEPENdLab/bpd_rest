# Extract ridge edge vals for NBS connected subnet ------------------------

# nbs_mnh <- load()
# 
# nbs.struct <- readMat("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/BNV_nodefiles/NBS_controls_greater_t_3.mat")
# nbs <- nbs.struct[["NBS"]]
# 
# nbsnetwork <- as.matrix(nbs[3,,]$con.mat[[1]][[1]])
basedir <- "~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"
setwd(basedir)
source("scripts/setup_globals.R")
inputs <- specify_inputs(thresh_weighted = "binary", 
                         fc_out_rm = FALSE, 
                         preproc_pipeline = "nosmooth_aroma_bp_nonaggr", 
                         reducemetrics = c("degree", "eigen.cent", "betweenness.node", "part.coeff", "gateway.coeff.btw",  "within.module.deg"),
                         rs_desired_log = logspace(log10(.01), log10(.02), 20)) 
#PEARSON
#inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, conn_method = "pearson", rs_desired_log = logspace(log10(.2), log10(.4), 20)) #leave blank for defaults, consult with function for deviations from default
for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

setwd(paste0(basedir, "/NBS_FEB_2018"))


nbsnetwork <- read.table("BPD_greater_t2.85_211nodes_300perms.txt")
edge.mat <- which(nbsnetwork!=0,arr.ind = T) 
#length(unique(edge.mat[,1], edge.mat[,2])) #re: nnodes = nedges.. this seems to not be the case. nedges = 65 nnodes = 41

colnames(edge.mat) <- c("node1","node2")

subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".RData", fd.scrub = TRUE, allowCache = TRUE)
allmats <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

edges_nbs_bpd <- data.frame(BPD = subj_info$BPD)
edge.vals <- as.data.frame(array(NA, c(length(edges_nbs_bpd[,1]), length(edge.mat[,1]))))
edges_nbs_bpd <- cbind(edges_nbs_bpd, edge.vals)
                        
for(s in 1:dim(allmats)[1]){
  new.row <- rep(NA, length(edge.mat[,1]))
  for(e in 1:length(edge.mat[,1])){
     new.row[e] <- allmats[s, edge.mat[e,1],edge.mat[e,2]] 
  }
  edges_nbs_bpd[s,c(1:length(edge.mat[,1]) +1)] <- new.row
}

colnames(edges_nbs_bpd) <- c("BPD", paste0("E_", seq(1,length(edge.mat[,1]),1)))

library(reshape2)
edges_melt <- reshape2::melt(edges_nbs_bpd, id = "BPD")
colnames(edges_melt) <- c("BPD", "Edge_number", "Ridge_r")

range(edges_melt$Ridge_r)
setwd(basedir)
pdf("figures/NBS_edge_boxplots_bpd_greater.pdf", width=8.5, height=11)
for(e in 1:length(edge.mat[,1])){
e_name <- paste0("E_", e)
edges_melt_filtered <- dplyr::filter(edges_melt, Edge_number == e_name)
boxplots_nbsEdges <- ggplot(edges_melt_filtered, aes(x = factor(BPD), y = Ridge_r, fill = factor(BPD))) + 
                              geom_boxplot() + 
                              coord_cartesian(ylim = range(edges_melt_filtered$Ridge_r)) +
                              labs(x = "Group Membership (BPD = 1)", 
                                   y = "Ridge Regression Value",
                                   title = "Distribution of edge weights based off of ridge regression values between nodes:",
                                   subtitle = paste0(edge.mat[e,1], " and ", edge.mat[e,2]))
plot(boxplots_nbsEdges)
}
dev.off()


