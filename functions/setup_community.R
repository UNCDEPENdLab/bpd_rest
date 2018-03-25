run_community_detection_on_agg <- function(allmats, agg.g, algorithm="louvain", density=NULL,  ...) {
  #Run community detection on mean graph (across all subjects in 3d array passed in)
  #allmats_log <- apply(allmats, c(1,2,3), function(x) { log(x+.05)})
  #NB. If no density threshold is applied, cluster_* runs weighted community detection
  #The weights should scaled such that higher values represent stronger ties. This is opposite of
  #betweenness, where inverse weights are used to denote distance.  
  
  if (!is.null(density)) { agg.g <- density_threshold(agg.g, density) }
 
  if (algorithm=="louvain") {
    comm <- cluster_louvain(agg.g)    
  } else if (algorithm=="fast_greedy") {
    comm <- cluster_fast_greedy(agg.g)
  } else if (algorithm=="infomap") {
    comm <- runinfomap(agg.g, ...) #compiled C++ version of infomap 
  } else { stop ("Not supported yet") }
  
  return(comm)
  
}


assign_communities <- function(allg, comm, attribute="community") {
  stopifnot(is.list(allg))
  stopifnot("communities" %in% class(comm)) #ensure that it's a community object from igraph
  if (is_igraph(allg[[1]])) {
    #assume that we have a list of subjects where each element is a graph (currently used for weighted graphs)
    allg <- lapply(allg, function(subj) {
          subj <- set_vertex_attr(subj, attribute, value=comm$membership)
          return(subj)
        })
  } else if (is.list(allg[[1]]) && is_igraph(allg[[1]][[1]])) {
    #assuming we have a nested list structure: [subjects] x [density thresholds]
    
    allg <- lapply(allg, function(subj) {
          subj_graphs <- lapply(1:length(subj), function(d) {
                subj[[d]] <- set_vertex_attr(subj[[d]], attribute, value=comm$membership)
                return(subj[[d]])
              })
        })
  }
} 

yeo7_community <- function(agg.g){
  #stopifnot(is.list(allg)), if interested we can input a list of single subject matrices 
  membership.file <- file.path(getwd(), "data", paste0(parcellation, "_membership.yeo.csv"))
  stopifnot(file.exists(membership.file))
  membership.yeo <- as.numeric(as.matrix(read.csv(membership.file)))
  names(membership.yeo) <- atlas$name
   
  yeo7_community <- make_clusters(agg.g, membership = membership.yeo, algorithm = "Yeo_etal_2011_7Networks")
  
  return(yeo7_community)
}

# hist(E(mean.g)$weight)
# E(mean.g)
# dev.off()
#
#if (use.infomap == 1){
#  mean.g.infomap <- readRDS(paste0(basedir, "/cachedRfiles/infomap_communitylist.rds"))
#} else{
#  ##output will be 20 binary graphs AVERAGED ACROSS SUBJECTS (at 1-20% density) with community membership
#  mean.g.community <- lapply(densities_desired, function(d) {
#        weights <- sort(E(mean.g)$weight, decreasing=TRUE)
#        threshold <- weights[length(V(mean.g))*(length(V(mean.g))-1)/2 * d]
#        gthresh <- delete.edges(mean.g, which(E(mean.g)$weight < threshold))
#        gthresh <- remove.edge.attribute(gthresh, "weight")
#        louv.graph <- cluster_louvain(gthresh)
#        V(gthresh)$community <- louv.graph$membership
#        
#        return(gthresh)
#      })
#}
