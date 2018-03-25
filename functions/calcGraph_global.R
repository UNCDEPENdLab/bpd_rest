####worker bee function which compiles global graph theoretical metrics
calcGraph_global <- function(origgraph, community_attr="community", modular_weighted = NULL) {
  #edgelist <- data.frame(get.edgelist(origgraph)) #not used currently
  globmetrics <- list()
  
  globmetrics[["edge_density"]] <- edge_density(origgraph) #does not use edge weights
  globmetrics[["characteristic_path_length"]] <- mean_distance(origgraph) #does not use edge weights
  globmetrics[["clustering_coefficient"]] <- transitivity(origgraph, type = "global") #uses edge weights if the graph has an edge weight attribute
  #globmetrics[["small_worldness"]] <- (mean_distance(origgraph)/mean_distance(rewire(origgraph, with = keeping_degseq(loops = FALSE, niter = 1000))))/(transitivity(origgraph, type = "global")/transitivity(rewire(origgraph, with = keeping_degseq(loops = FALSE, niter = 1000)), type = "global")) 
  
  if (!is.null(vertex_attr(origgraph, community_attr))) {
    if (is.null(modular_weighted)) {
      globmetrics[["modularity"]] <- modularity(origgraph, vertex_attr(origgraph, community_attr))
    } else {
      globmetrics[["modularity"]] <- modularity(origgraph, vertex_attr(origgraph, community_attr), weights = E(origgraph)$weight)
    }
  }
  
  return(globmetrics)
}
