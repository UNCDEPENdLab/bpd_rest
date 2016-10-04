
##############################################################################################
###############################################
####worker bee function which compiles graph metrics accross densities, after procedure looks good to go, throw a loop around it and run across participants
calcGraph <- function(origgraph, density) {
  edgelist <- data.frame(get.edgelist(origgraph) , weight=round(E(origgraph)$weight, 3 ))
  cut <- quantile(edgelist$weight, 1 - density)            ##sparsity

  cutgraph <- delete.edges(origgraph, E(origgraph)[weight<= cut])
  deg <- degree(cutgraph, mode = "all")
  deg0 <- which(deg<=1)
  cutgraph <- delete.vertices(cutgraph, deg0)
  allmetrics <- list()
  allmetrics[["nodes"]] <- V(cutgraph)
  #allmetrics[["nodes"]]$label <- node.file$V6 #if you need this, pass in as argument
  allmetrics[["edges"]] <- E(cutgraph)
  
  allmetrics[["diameter"]] <- get_diameter(cutgraph, directed = F, weights = NA)
  allmetrics[["efficiency"]] <- graph.efficiency(cutgraph, type = "global")
  eff<-1/(shortest.paths(cutgraph))
  eff[!is.finite(eff)]<-0
  allmetrics[["global.efficiency"]] <- mean(eff,na.rm=TRUE)
  allmetrics[["transitivity"]] <- transitivity(cutgraph, type = "global")
  mean.dist <- mean_distance(cutgraph, directed = FALSE)
  rand.sw <- sample_gnp(length(V(cutgraph)), 1/20)  ##generate random graph
  rand.meandist <- mean_distance(rand.sw, directed = FALSE)
  rand.transglobal <- transitivity(rand.sw, type = "globalundirected", isolates = "zero")
  allmetrics[["small.worldness"]]<- (allmetrics[["transitivity"]]/rand.transglobal)/(mean.dist/rand.meandist)
  allmetrics[["degree"]] <- deg
  allmetrics[["strength"]] <- strength(cutgraph, vids = V(cutgraph))
  allmetrics[["eigen.cent"]] <- centr_eigen(cutgraph, directed = FALSE)
  allmetrics[["closeness"]] <- centr_clo(cutgraph, mode = "all")
  allmetrics[["betweenness.node"]] <- centr_betw(cutgraph)
  allmetrics[["edge.betweenness"]] <- edge_betweenness(cutgraph, E(cutgraph))
  
  louv <- cluster_louvain(cutgraph, weights = E(cutgraph)$weights)
  cutgraph$louv.membership <- louv$membership
  allmetrics[["louv.membership"]] <- cutgraph$louv.membership
  allmetrics[["louv.intrad"]] <- rep(NA, length(louv))
  for (i in 1:length(louv)) {
    subg <- induced_subgraph(cutgraph, which(cutgraph$louv.membership == i))
    subden <- ecount(subg)/ecount(cutgraph)
    allmetrics[["louv.intrad"]] <- subden
  }
  allmetrics[["louv.unique"]] <- unique(allmetrics[["louv.membership"]])
  allmetrics[["part.coeff.louv"]] <- part.coeff(cutgraph, allmetrics[["louv.membership"]])
  
  lead.eig <- leading.eigenvector.community(cutgraph, weights = E(cutgraph)$weights)
  cutgraph$eigen.membership <- membership(lead.eig)
  allmetrics[["eigen.membership"]] <- cutgraph$eigen.membership
  allmetrics[["eigen.intrad"]] <- rep(NA, length(unique(allmetrics[["eigen.membership"]])))
  for (i in 1:length(unique(allmetrics[["eigen.membership"]]))) {
    subg <- induced_subgraph(cutgraph, which(cutgraph$eigen.membership == i))
    subden <- (ecount(subg)/ecount(cutgraph))
    allmetrics[["eigen.intrad"]] <- subden
  }
  allmetrics[["eigen.unique"]] <- unique(allmetrics[["eigen.membership"]])
  allmetrics[["part.coeff.eig"]] <- part.coeff(cutgraph, allmetrics[["eigen.unique"]])
  
  
  return(allmetrics)
}
