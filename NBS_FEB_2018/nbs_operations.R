
# NBS Operations ----------------------------------------------------------


basedir <- "~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"
setwd(basedir)
source("scripts/setup_globals.R")
inputs <- specify_inputs(thresh_weighted = "binary", 
                         # conn_method = "cor.shrink",
                         fc_out_rm = FALSE, 
                         preproc_pipeline = "nosmooth_aroma_bp_nonaggr", 
                         reducemetrics = c("degree", "eigen.cent", "betweenness.node", "part.coeff", "gateway.coeff.btw",  "within.module.deg"),
                         rs_desired_log = logspace(log10(.01), log10(.02), 20)) 
#PEARSON
#inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, conn_method = "pearson", rs_desired_log = logspace(log10(.2), log10(.4), 20)) #leave blank for defaults, consult with function for deviations from default
for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

allmats <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

setwd(paste0(basedir, "/NBS_FEB_2018"))



nbsnetwork <- read.table("BPD_greater_t3_46nodes_5000perms.txt")

nbs_ts <- nbsnetwork[which(nbsnetwork !=0)]

###outputs binary NBS
outputdir <-  paste0("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/NBS_FEB_2018/plotting_files/")
write.table(nbsnetwork, file = file.path(outputdir, "nbs_t3_46nodes.edge"), row.names = FALSE, col.names = FALSE)


####combine with t stats for weighted vizualization
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".RData", fd.scrub = TRUE, allowCache = TRUE)
edge.comps <- run_edge_comparisons(subj_info)
t.stats <- as.matrix(edge.comps[[2]]) #t-values on upper and p vals on lower
t.stats[lower.tri(t.stats)] <- 0


nbs.t.weighted <- data.matrix(t.stats*nbsnetwork)
#0 out any remaining NAs..
for (i in 1:length(nbs.t.weighted[,1])){
  for(j in 1:length(nbs.t.weighted[,1])){
    if(is.na(nbs.t.weighted[i,j])) {nbs.t.weighted[i,j] <- 0}  
  }
}



write.table(nbs.t.weighted, file = file.path(outputdir, "nbs_edge_lowconn_t3_5000perms_weighted.edge"), append = TRUE,row.names = FALSE, col.names = FALSE)

nbs.t.weighted[lower.tri(nbs.t.weighted)] <- t(nbs.t.weighted)[lower.tri(nbs.t.weighted)]
#isSymmetric(nbs.t.weighted) 
nbs.t.weighted <- abs(nbs.t.weighted)

E(nbs.graph)$weight

nbs.graph <- graph.adjacency(nbs.t.weighted, mode = "undirected", weighted = TRUE, diag = FALSE)
V(nbs.graph)$name <- atlas$name
nbs.graph <- tagGraph(nbs.graph, atlas)
nbs.graph.solo <- delete_vertices(nbs.graph, degree(nbs.graph) == 0)


nodal_stats_nbs <- list()
nodal_stats_nbs[["degree"]] <- sort(degree(nbs.graph.solo, v = V(nbs.graph.solo)))
#nodal_stats_nbs[["ev.cent"]] <- sort(round(eigen_centrality(nbs.graph.solo)$vector, 2))
#nodal_stats_nbs[["betweenness"]] <- betweenness(nbs.graph.solo, v = V(nbs.graph.solo))
nodal_stats_nbs[["strength"]] <- sort(strength(nbs.graph.solo, v = V(nbs.graph.solo)))

nodal_stats_nbs_low_df <- do.call(cbind, nodal_stats_nbs)

nodal_stats_nbs_low_df <- tibble::rownames_to_column(as.data.frame(nodal_stats_nbs_low_df), var = "name")

anat_merge <- atlas[,c("name","anat_label")]

nodal_stats_nbs_low_df <- dplyr::inner_join(nodal_stats_nbs_low_df, anat_merge, by = "name")
nodal_stats_nbs_low_df[order(-nodal_stats_nbs_low_df[,"strength"]),]

membership_df <- get(load("/Users/nth7/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/cache_protected/membership_df.RData"))
membership_df <- membership_df %>%
  mutate(membership=factor(membership, levels=1:7, labels=c("VIS", "SOMMOT", "DORSATTN", "SALVENTATTN", "LIMBIC", "FPN", "DMN")))
colnames(nodal_stats_nbs_low_df) <- c("node", "degree", "strength", "region_name")
nodal_stats_nbs_low_df <- nodal_stats_nbs_low_df %>% left_join(membership_df, by = "node")

nodal_stats_nbs_low_df <- nodal_stats_nbs_low_df[order(-nodal_stats_nbs_low_df[,"strength"]),]
resultsdir <- paste0(basedir, "results")
write.csv(nodal_stats_nbs_low_df, file = file.path(resultsdir,  "nodal_stats_nbs_low_df.csv"))
# Hyperconnectivity results -----------------------------------------------
#read nbs HYPERCONNECTIVITY results, subnetwork in which BPD participants are more connected

nbs.struct.hyper <- readMat("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/BNV_nodefiles/NBS_bpd_greater_t_2.95.mat")
nbs.hyp <- nbs.struct.hyper[["NBS"]]

nbsnetwork.hyp <- as.matrix(nbs.hyp[3,,]$con.mat[[1]][[1]])

###outputs binary NBS
outputdir <-  paste0("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/BNV_nodefiles/schaefer422_ridge/")
write.table(nbsnetwork.hyp, file = file.path(outputdir, "hyper.nbs.edge.txt"), row.names = FALSE, col.names = FALSE)

##for weighted analysis
####combine with t stats for weighted vizualization
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".RData", fd.scrub = TRUE, allowCache = TRUE)
edge.comps <- run_edge_comparisons(subj_info)
t.stats <- as.matrix(edge.comps[[2]]) #t-values on upper and p vals on lower
t.stats[lower.tri(t.stats)] <- 0



nbs.hyp.t.weighted <- data.matrix(t.stats*nbsnetwork.hyp)
#0 out any remaining NAs..
for (i in 1:length(nbs.hyp.t.weighted[,1])){
  for(j in 1:length(nbs.hyp.t.weighted[,1])){
    if(is.na(nbs.hyp.t.weighted[i,j])) {nbs.hyp.t.weighted[i,j] <- 0}  
    # message("setting NA value to zero at:", i, " ",j)
  }
}



write.table(nbs.hyp.t.weighted, file = file.path(outputdir, "nbs.edge.high.conn.t.weighted.txt"), row.names = FALSE, col.names = FALSE)


nbs.hyp.t.weighted[lower.tri(nbs.hyp.t.weighted)] <- t(nbs.hyp.t.weighted)[lower.tri(nbs.hyp.t.weighted)]
#isSymmetric(nbs.t.weighted) 

#nbs.t.weighted <- abs(nbs.t.weighted) #this step should be unneccesary given that t values should be positive in this case



nbs.hyp.graph <- graph.adjacency(nbs.hyp.t.weighted, mode = "undirected", weighted = TRUE, diag = FALSE)
V(nbs.hyp.graph)$name <- atlas$name
nbs.hyp.graph <- tagGraph(nbs.hyp.graph, atlas)
nbs.hyp.graph.solo <- delete_vertices(nbs.hyp.graph, degree(nbs.hyp.graph) == 0)


nodal_stats_nbs.hyp <- list()
nodal_stats_nbs.hyp[["degree"]] <- sort(degree(nbs.hyp.graph.solo, v = V(nbs.hyp.graph.solo)))
nodal_stats_nbs.hyp[["ev.cent"]] <- sort(round(eigen_centrality(nbs.hyp.graph.solo)$vector, 2))
nodal_stats_nbs.hyp[["betweenness"]] <- betweenness(nbs.hyp.graph.solo)


nodal_stats_nbs_high_df <- do.call(cbind, nodal_stats_nbs.hyp)

nodal_stats_nbs_high_df <- tibble::rownames_to_column(as.data.frame(nodal_stats_nbs_high_df), var = "name")

anat_merge <- atlas[,c("name","anat_label")]

nodal_stats_nbs_high_df <- dplyr::inner_join(nodal_stats_nbs_high_df, anat_merge, by = "name")
nodal_stats_nbs_high_df[order(-nodal_stats_nbs_high_df[,"ev.cent"]),]

resultsdir <- paste0(basedir, "/results")
write.csv(nodal_stats_nbs_high_df, file = file.path(resultsdir,  "/schaefer422_ridge/nodal_stats_nbs_high_df.csv"))



# Work bench --------------------------------------------------------------

# 
# class(nbs.graph.solo)
# a <- cluster_leading_eigen(nbs.graph.solo)
# plot(a, layout = coords)
# plot(nbs.graph.solo)
# 
# maximal.cliques(nbs.graph.solo)
# 
# ggplot(ggnetwork(nbs.graph.solo)) + 
#   
# ##mess with ggnetwork example
# data(email, package = 'geomnet')
# em.cet <- as.character(
#   email$nodes$CurrentEmploymentType)
# names(em.cet) = email$nodes$label
# # remove the emails sent to all employees
# edges <- subset(email$edges, nrecipients < 54)
# # create network
# str(em)
# em.net <- edges[, c("From", "to") ]
# em.net <- network(em.net, directed = TRUE)
# 
