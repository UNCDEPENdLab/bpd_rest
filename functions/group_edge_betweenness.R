group_edge_betweenness <- function(allg, subj_info, allowCache = TRUE, weighted = TRUE){
  
  if (weighted){
    betweenness.subs <- list()
    for (sub in 1:length(allg)){
      eb <- edge_betweenness(allg[[sub]], e = E(allg[[sub]]), directed = FALSE)
      betweenness.subs[[sub]] <- eb
    }
  } else {
    betweenness.subs <- list()
    for (sub in 1:length(allg)){
      eb <- edge_betweenness(allg[[sub]], e = E(allg[[sub]]), directed = FALSE)
      betweenness.subs[[sub]] <- eb
    }
  }
  
}