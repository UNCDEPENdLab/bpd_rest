####worker bee function which compiles nodal graph metrics including nodal metrics with reference to a $community attribute of the graph 
calcGraph_nodal <- function(origgraph, community_attr="community", weighted = FALSE) {
  #edgelist <- data.frame(get.edgelist(origgraph)) #not used currently
  allmetrics <- list()
  
  ####BINARY
  if(weighted == FALSE){
  ##centrality
  allmetrics$local.clustering <- transitivity(origgraph, type = "local")
  allmetrics$degree <- degree(origgraph, v = V(origgraph))
  evcent<- eigen_centrality(origgraph, directed = FALSE)
  allmetrics$eigen.cent <- evcent$vector
  allmetrics$closeness <- closeness(origgraph, mode = "all", normalized=TRUE)#, weights=1/E(origgraph)$weight)
  allmetrics$betweenness.node <- betweenness(origgraph)#, weights=1/E(origgraph)$weight)
  allmetrics$page.rank <- page.rank(origgraph, algo="prpack")$vector
  allmetrics$leverage.cent <- centr_lev(origgraph)

  #only compute community-related statistics if the expected community attribute is available
  if (!is.null(vertex_attr(origgraph, community_attr))) {
    allmetrics$community.membership <- vertex_attr(origgraph, community_attr)
    
    ##nodal within vs between modules
    wibw <- wibw_module_degree(origgraph, community_attr=community_attr)
    allmetrics[["within.module.deg.zscore"]] <- wibw$z_within
    allmetrics[["between.module.deg.zscore"]] <- wibw$z_between
    allmetrics[["within.module.deg"]] <- wibw$Ki
    allmetrics[["between.module.deg"]] <- wibw$Bi
    allmetrics[["part.coeff"]] <- part_coeff(origgraph, allmetrics$community.membership)#, use.parallel = FALSE)
    allmetrics[["gateway.coeff.btw"]] <- gateway_coeff(origgraph, allmetrics$community.membership, centr = "btwn.cent")
    allmetrics[["gateway.coeff.degree"]] <- gateway_coeff(origgraph, allmetrics$community.membership, centr = "degree")    
    }
  } else {
    #####WEIGHTED
    ##centrality
    allmetrics$local.clustering <- transitivity(origgraph, type = "local") #defaults to read edge weights when available
    # allmetrics$degree <- degree(origgraph, v = V(origgraph))
    allmetrics$strength <- strength(origgraph, v = V(origgraph))
    evcent<- eigen_centrality(origgraph, directed = FALSE) #defaults to read edge weights when available
    allmetrics$eigen.cent <- evcent$vector
    allmetrics$closeness <- closeness(origgraph, mode = "all", normalized=TRUE, weights=1/E(origgraph)$weight) 
    allmetrics$betweenness.node <- betweenness(origgraph, weights=1/E(origgraph)$weight)
    allmetrics$page.rank <- page.rank(origgraph, algo="prpack")$vector
    allmetrics$leverage.cent <- centr_lev(origgraph)
    
    #only compute community-related statistics if the expected community attribute is available
    if (!is.null(vertex_attr(origgraph, community_attr))) {
      allmetrics$community.membership <- vertex_attr(origgraph, community_attr)
      
      ##nodal within vs between modules
      wibw <- wibw_module_degree(origgraph, community_attr=community_attr)
      allmetrics[["within.module.deg.zscore"]] <- wibw$z_within
      allmetrics[["between.module.deg.zscore"]] <- wibw$z_between
      allmetrics[["within.module.deg"]] <- wibw$Ki
      allmetrics[["between.module.deg"]] <- wibw$Bi
      allmetrics[["part.coeff"]] <- part_coeff(origgraph, allmetrics$community.membership)#, use.parallel = FALSE)
      allmetrics[["gateway.coeff.btw"]] <- gateway_coeff_NH(origgraph, allmetrics$community.membership, centr = "btwn.cent") #does not inverse weights for betweenness in original code, so I quickly changed and added  as new function
      allmetrics[["gateway.coeff.degree"]] <- gateway_coeff(origgraph, allmetrics$community.membership, centr = "degree")    
    }
  }

  return(allmetrics)
}





