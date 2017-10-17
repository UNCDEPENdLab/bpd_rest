setup_graphs <- function(adjmats_all, aggfun = NULL, agg.rm.neg = TRUE, allowCache=TRUE, ncpus=4, dens_rm = .10) {
  #expects a subjects x ROIs x ROIs array OR list of density x ROIs x ROIs arrays (per subject) for dens.clime
  #creates igraph objects from adjacency matrix and label nodes V1, V2,...
  #looks up atlas from global environment
  #dens.rm correspoonds to the first density we will keep when calculating allg_noneg in the case of the dens.clime conn_method
  suppressMessages(require(igraph))
  suppressMessages(require(foreach))
  suppressMessages(require(doSNOW))
  
  # browser()
  
  if(conn_method == "dens.clime_partial"){stopifnot(length(atlas$name) == length(adjmats_all[[1]][1,1,]))} else{
  stopifnot(length(atlas$name) == dim(adjmats_all)[2]) #number of nodes must match between atlas and adjmats
  }
  
  
  
  if(conn_method == "dens.clime_partial"){
    #standard
    expectFile <- file.path(basedir, "cache", paste0("allg_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
    #Oct17 debug
    # expectFile <- file.path(basedir, "cache", paste0("allg_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
    if (file.exists(expectFile) && allowCache==TRUE) {
      message("Loading dens.clime weighted graph objects list from file: ", expectFile)
      load(expectFile)
    } else {
     
  allg <- lapply(adjmats_all, function(sub){
    sub_graphs <- apply(sub, 1, function(den){
      g <- graph.adjacency(den, mode = "undirected", weighted = TRUE, diag = FALSE)
      V(g)$name <- atlas$name
      
      g <- tagGraph(g, atlas) #populate all attributes from atlas to vertices
      return(g)
    })
    return(sub_graphs)
  })
  
  #label allg graphs with subject IDs 
  
  for (sub in 1:length(allg)) { 
    for (den in 1:length(densities_desired))
    allg[[sub]][[den]]$id <- names(adjmats_all)[sub]
    }
    
  
  save(file = expectFile, allg)
    }
    
    expectFile <- file.path(basedir, "cache", paste0("weightedgraphsnoneg_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
    if (file.exists(expectFile) && allowCache==TRUE) {
      message("Loading weighted non-negative graph objects from file: ", expectFile)
      load(expectFile)
    } else {
      
      ##remove negative correlations between nodes, if this is run on pearson, the number of edges remaining will be different across subjs 
      allg_noneg <- lapply(allg, function(sub) {
        sub_graphs <- lapply(sub, function(den){
          delete.edges(den, which(E(den)$weight < 0)) })
        })
      save(file=expectFile, allg_noneg) #save cache
    }
    
  #### Generate aggregate graph, defaults to mean of all adjmats with neg weights removed
  ridge.mean <- file.path(basedir, "cache", paste0("agg.g_", parcellation, "_", preproc_pipeline, "_ridge.net_partial.RData"))
  if (file.exists(ridge.mean)) {
    message("Loading aggregated ridge graph (DEFAULT MEAN + NO NEG) from file: ", ridge.mean)
    load(ridge.mean)
  } else {message("no agg.g to pull")}

  
  } else {
  ##create graph objects if inputted adjmats_all file is not a list of dens weighted graphs
  #### Setup basic weighted graphs if not dens.clime
  expectFile <- file.path(basedir, "cache", paste0("weightedgraphs_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading weighted graph objects from file: ", expectFile)
    load(expectFile)
  } else {
    
    allg <- apply(adjmats_all, 1, function(sub) {
      g <- graph.adjacency(sub, mode="undirected", weighted=TRUE, diag=FALSE)
      V(g)$name <- atlas$name

      g <- tagGraph(g, atlas) #populate all attributes from atlas to vertices
      return(g)
    })
    
    #add id to graphs
    for (i in 1:length(allg)) { allg[[i]]$id <- dimnames(adjmats_all)[["id"]][i] }
    save(file=expectFile, allg) #save cache
    
  }
  #### Setup weighted graphs, removing negative edges
  expectFile <- file.path(basedir, "cache", paste0("weightedgraphsnoneg_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading weighted non-negative graph objects from file: ", expectFile)
    load(expectFile)
  } else {
    ##remove negative correlations between nodes, if this is run on pearson, the number of edges remaining will be different across subjs 
    allg_noneg <- lapply(allg, function(g) { delete.edges(g, which(E(g)$weight < 0)) })
    save(file=expectFile, allg_noneg) #save cache
  }
  
  #### Generate aggregate graph, defaults to mean of all adjmats with neg weights removed
  expectFile <- file.path(basedir, "cache", paste0("agg.g_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading aggregated graph (DEFAULT MEAN + NO NEG) from file: ", expectFile)
    load(expectFile)
  } else {
    agg.g <- generate_agg.g(allmats, rm.neg = agg.rm.neg)
    save(file=expectFile, agg.g) #save cache
  }
  
  ### DENSITY THRESHOLDING: binarize and threshold graphs at densities ranging from 1-20%
  #uses densities_desired from setup_globals.R
  
  expectFile <- file.path(basedir, "cache", paste0("dthreshgraphs_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading density-thresholded graph objects from file: ", expectFile)
    load(expectFile)
  } else {
    setDefaultClusterOptions(master="localhost")
    clusterobj <- makeSOCKcluster(ncpus)
    registerDoSNOW(clusterobj)
    
    #make sure .inorder = TRUE to keep subject sorting the same in a parallel loop
    allg_density <- foreach(g=allg_noneg, .packages=c("igraph", "tidyr"), .export=c("densities_desired", "density_threshold"), .inorder = TRUE) %dopar% {
      #nnodes <- length(V(g))
      #maxedges <- (nnodes*(nnodes-1))/2
      
      dgraphs <- lapply(densities_desired, density_threshold, g=g) #implicitly passes the dth element of densities_desired to density_threshold function      
      names(dgraphs) <- paste0("d", densities_desired)
      return(dgraphs)
      }
    
    names(allg_density) <- names(allg_noneg) #propagate IDs at names to density-thresholded list
    stopCluster(clusterobj)
    
    #each element of allg_density is a list of 20 binary graphs for that subject at 1-20% density
    save(file=expectFile, allg_density)
    }
  return(list(allg=allg, allg_noneg=allg_noneg, allg_density=allg_density, agg.g = agg.g))
  }
  return(list(allg, allg_noneg, agg.g))
}
