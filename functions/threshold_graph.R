#simple function to apply thresholding to a weighted graph, either based on density or strength threshold
threshold_graph <- function(g, w, method="density", rmweights=TRUE) {
  stopifnot(is_igraph(g))
  stopifnot(is.numeric(w))
  
  #Obtains desired density given graph diameter
  if (method=="density") {
    weights <- sort(E(g)$weight, decreasing=TRUE)
    threshold <- weights[length(V(g))*(length(V(g))-1)/2 * w]
  } else {
    threshold <- w #just use strength-based threshold
  }
  
  gthresh <- delete.edges(g, which(E(g)$weight < threshold))
  if (rmweights) {
    E(gthresh)$w_hidden <- E(gthresh)$weight #retain weight attribute for additional analyses (but remove "weight" to avoid igraph thinking it's a weighted graph)
    gthresh <- remove.edge.attribute(gthresh, "weight")
  }
  gthresh$wthresh <- threshold #copy threshold weight into object for tracking  
  if (method=="density") { gthresh$target_density <- w } #copy density into object for tracking  
  return(gthresh)  
}

threshold_glist <- function(glist, thresholds, method="density", ncores=4, ...) {
  allg_thresh <- lapply(thresholds, function(t, glist) {
    lapply(glist, function(group) {
      mclapply(group, function(g) {
        return(threshold_graph(g, t, method=method, ...))
      }, mc.cores=ncores)
    })
  }, glist)
  
  names(allg_thresh) <- as.character(paste(method, thresholds, sep="_"))
  return(allg_thresh)
}