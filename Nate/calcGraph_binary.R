####worker bee function which compiles nodal graph metrics including nodal metrics with reference to an input community structure
calcGraph_binary_nodal <- function(origgraph) {
  edgelist <- data.frame(get.edgelist(origgraph))
  origgraph <- tagGraph(origgraph, nodefile="~/Box Sync/RS_BPD_graph/power269_masterlookup_precull.csv", atlasname="power269")
  allmetrics <- list()
  
  allmetrics[["origgraph"]] <- origgraph
  allmetrics[["local.clustering"]] <- transitivity(origgraph, type = "local")
  allmetrics[["degree"]] <- degree(origgraph, v = V(origgraph))
  evcent<- eigen_centrality(origgraph, directed = FALSE)
  allmetrics[["eigen.cent"]] <- evcent$vector
  allmetrics[["closeness"]] <- closeness(origgraph, mode = "all", normalized=TRUE)#, weights=1/E(origgraph)$weight)
  allmetrics[["betweenness.node"]] <- betweenness(origgraph)#, weights=1/E(origgraph)$weight)
  allmetrics[["page.rank"]] <- page.rank(origgraph, algo="prpack")$vector

  allmetrics[["community.membership"]] <- V(origgraph)$community
 
  allmetrics[["within.module.deg.zscore"]] <- within_module_deg_z_score(origgraph, V(origgraph)$community, use.parallel=FALSE)
  allmetrics[["part.coeff"]] <- part.coeff(origgraph, V(origgraph)$community, use.parallel = FALSE)

  return(allmetrics)
}

####worker bee function which compiles global graph theoretical metrics
calcGraph_binary_global <- function(origgraph) {
  edgelist <- data.frame(get.edgelist(origgraph))
  origgraph <- tagGraph(origgraph, nodefile="~/Box Sync/RS_BPD_graph/power269_masterlookup_precull.csv", atlasname="power269")
  globmetrics <- list()

  globmetrics[["origgraph"]] <- origgraph
  #globmetrics[["shortest_path"]] <- shortest.paths(origgraph)
  globmetrics[["characteristic_path_length"]] <- mean_distance(origgraph)
  #globmetrics[["number_triangles"]] <- count_triangles(origgraph)
  globmetrics[["clustering_coefficient"]] <- transitivity(origgraph, type = "global")
  #globmetrics[["small_worldness"]] <- small.world(origgraph)
  
  return(globmetrics)
}
