generate_agg.g <- function(allmats, aggfun =mean, rm.neg = TRUE){
  agg.adjmat <- apply(allmats, c(2,3), aggfun, na.rm = TRUE)
  agg.g <- graph.adjacency(agg.adjmat, mode = "lower", weighted = TRUE, diag = FALSE)
  #remove negative weights
  if(rm.neg == TRUE) {agg.g <- delete.edges(agg.g, which(E(agg.g)$weight < 0))} 
  return(agg.g)
}