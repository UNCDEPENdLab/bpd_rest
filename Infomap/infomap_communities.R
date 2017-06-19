runinfomap <- function(agg.g, hierarchical=FALSE, verbose=FALSE, infomap_dir=file.path(basedir, "Infomap")) {
  library(igraph)
  cwd <- getwd()
  stopifnot(file.exists(infomap_dir))
  stopifnot(is_igraph(agg.g))
  setwd(infomap_dir)
  weighted <- !is.null(E(agg.g)$weight)
  
  ## Load Infomap library
  source("load-infomap.R")
  
  ### Zachary's karate club (for testing)
  ##g <- graph.famous("Zachary")
  
  prefix <- ifelse(verbose, "", "--silent")
  if (hierarchical) {
    conf <- init(paste(prefix, "-N 500 -u -M 50")) #multi-level partition
  } else {
    conf <- init(paste(prefix,"--two-level --silent -N 500 -u -M 50")) #two-level partition, allow 50 iterations of module assignment loop    
  }

  network <- Network(conf);
  
  ## Add links to Infomap network from igraph data
  edgelist <- as_edgelist(agg.g)
  
  #this assumes that nodes are labeled V<number>
  edgelist <- apply(edgelist, c(1,2), function(x) { as.numeric(sub("V", "", x, fixed=TRUE)) - 1 }) ##convert to 0-based numeric edgelist
  if (weighted) { edgelist <- cbind(edgelist, E(agg.g)$weight) } #tack on third column
  
  #add links between e1 and e2 (dummy does nothing -- just need $addLink to add onto network object)
  if (weighted) {
    dummy <- apply(edgelist, 1, function(e) network$addLink(e[1], e[2], e[3]))
  } else {
    dummy <- apply(edgelist, 1, function(e) network$addLink(e[1], e[2]))  
  }
  
  network$finalizeAndCheckNetwork(TRUE, vcount(agg.g))

  if (weighted) {
    if (abs(sum(E(agg.g)$weight) - network$totalLinkWeight()) > .001) { #make sure weights are added!
      stop("Weight check failed.")
    } else { 
      message("Weighted edge check complete!") 
    }
  }
  
  #network$totalLinkWeight() #has sums of weights on edges
  cat("Created network with", network$numNodes(), "nodes and", network$numLinks(), "links.\n")
  
  tree <- HierarchicalNetwork(conf)
  
  run(network, tree);
  
  cat("Partitioned network in", tree$numTopModules(), "modules with codelength", tree$codelength(), "bits:\n")
  
  clusterIndexLevel <- 1 ## 1, 2, ... or -1 for top, second, ... or lowest cluster level
  leafIt <- tree$leafIter(clusterIndexLevel)
  modules <- integer(length = network$numNodes())
  
  while (!leafIt$isEnd()) {
    modules[leafIt$originalLeafIndex + 1] = leafIt$moduleIndex() + 1
    leafIt$stepForward()
  }
  
  ## Create igraph community data
  comm <- create.communities(agg.g, modules, algorithm = 'Infomap')
  #print(comm)
  ## Plot communities and network
  #pdf(paste0("meang_infomap_", densities[i], ".pdf"), width=14, height=14)
  #plot(comm, agg.g)
  #dev.off()
  #table(commlist[[10]]$membership)
    
  return(comm)
  
  setwd(cwd)
}



