##########NBS operations
setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
source("functions/setup_globals.R")

#read nbs DISCONNECTIVITY results
library(R.matlab)
nbs.struct <- readMat("/Users/nth7/Documents/MATLAB/NBS1.2/NBS_t_3.0.mat")
nbs <- nbs.struct[["NBS"]]

nbsnetwork <- as.matrix(nbs[3,,]$con.mat[[1]][[1]])

###outputs binary NBS
outputdir <-  paste0("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/BNV_nodefiles/schaefer422_ridge/")
write.table(nbsnetwork, file = file.path(outputdir, "nbs.edge.txt"), row.names = FALSE, col.names = FALSE)


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



write.table(nbs.t.weighted, file = file.path(outputdir, "nbs.edge.low.conn.t.weighted.txt"), row.names = FALSE, col.names = FALSE)

nbs.t.weighted[lower.tri(nbs.t.weighted)] <- t(nbs.t.weighted)[lower.tri(nbs.t.weighted)]
#isSymmetric(nbs.t.weighted) 
nbs.t.weighted <- abs(nbs.t.weighted)



nbs.graph <- graph.adjacency(nbs.t.weighted, mode = "undirected", weighted = TRUE, diag = FALSE)
V(nbs.graph)$name <- atlas$name
nbs.graph <- tagGraph(nbs.graph, atlas)
nbs.graph.solo <- delete_vertices(nbs.graph, degree(nbs.graph) == 0)


data(email, package = 'geomnet')



degree.nbs <- degree(nbs.graph.solo, v = V(nbs.graph.solo))
sort(degree.nbs)

eigen.nbs <- sort(round(eigen_centrality(nbs.graph.solo)$vector, 2))
betweenness.nbs <- betweenness(nbs.graph.solo)

class(nbs.graph.solo)

a <- cluster_leading_eigen(nbs.graph.solo)
plot(a, layout = coords)
plot(nbs.graph.solo)

maximal.cliques(nbs.graph.solo)

ggplot(ggnetwork(nbs.graph.solo)) + 
  
##mess with ggnetwork example
data(email, package = 'geomnet')
em.cet <- as.character(
  email$nodes$CurrentEmploymentType)
names(em.cet) = email$nodes$label
# remove the emails sent to all employees
edges <- subset(email$edges, nrecipients < 54)
# create network
str(em)
em.net <- edges[, c("From", "to") ]
em.net <- network(em.net, directed = TRUE)

