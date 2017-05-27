####worker bee function which compiles nodal graph metrics including nodal metrics with reference to a $community attribute of the graph 
calcGraph_binary_nodal <- function(origgraph) {
  edgelist <- data.frame(get.edgelist(origgraph))
  origgraph <- tagGraph(origgraph, nodefile="~/Box Sync/RS_BPD_graph/power269_masterlookup_shift_nate.csv", atlasname="power269")
  allmetrics <- list()
  
  ##centrality
  allmetrics[["origgraph"]] <- origgraph
  allmetrics[["local.clustering"]] <- transitivity(origgraph, type = "local")
  allmetrics[["degree"]] <- degree(origgraph, v = V(origgraph))
  evcent<- eigen_centrality(origgraph, directed = FALSE)
  allmetrics[["eigen.cent"]] <- evcent$vector
  allmetrics[["closeness"]] <- closeness(origgraph, mode = "all", normalized=TRUE)#, weights=1/E(origgraph)$weight)
  allmetrics[["betweenness.node"]] <- betweenness(origgraph)#, weights=1/E(origgraph)$weight)
  allmetrics[["page.rank"]] <- page.rank(origgraph, algo="prpack")$vector

  allmetrics[["community.membership"]] <- V(origgraph)$community
 
  ##nodal within vs between modules
  wibw <- wibw_module_degree(origgraph)
  allmetrics[["within.module.deg.zscore"]] <- wibw$z_within
  allmetrics[["between.module.deg.zscore"]] <- wibw$z_between
  allmetrics[["part.coeff"]] <- part_coeff(origgraph, V(origgraph)$community)#, use.parallel = FALSE)
  allmetrics[["gateway.coeff.btw"]] <- gateway_coeff(origgraph, V(origgraph)$community, centr = "btwn.cent")
  allmetrics[["gateway.coeff.deg"]] <- gateway_coeff(origgraph, V(origgraph)$community, centr = "degree")

  return(allmetrics)
}

####worker bee function which compiles global graph theoretical metrics
calcGraph_binary_global <- function(origgraph) {
  edgelist <- data.frame(get.edgelist(origgraph))
  origgraph <- tagGraph(origgraph, nodefile="~/Box Sync/RS_BPD_graph/power269_masterlookup_shift_nate.csv", atlasname="power269")
  globmetrics <- list()

  globmetrics[["origgraph"]] <- origgraph
  globmetrics[["characteristic_path_length"]] <- mean_distance(origgraph)
  globmetrics[["clustering_coefficient"]] <- transitivity(origgraph, type = "global")
  globmetrics[["small_worldness"]] <- (mean_distance(origgraph)/mean_distance(rewire(origgraph, with = keeping_degseq(loops = FALSE, niter = 1000))))/(transitivity(origgraph, type = "global")/transitivity(rewire(origgraph, with = keeping_degseq(loops = FALSE, niter = 1000)), type = "global"))
  globmetrics[["modularity"]] <- modularity(origgraph, vertex_attr(origgraph)$community)
  
  return(globmetrics)
}
