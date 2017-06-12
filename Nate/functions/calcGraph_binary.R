####worker bee function which compiles nodal graph metrics including nodal metrics with reference to a $community attribute of the graph 
calcGraph_binary_nodal <- function(origgraph, community_attr="wcomm_louvain") {
  #edgelist <- data.frame(get.edgelist(origgraph)) #not used currently
  allmetrics <- list()

  ##centrality
  #allmetrics[["origgraph"]] <- origgraph #unclear why we should duplicate the graph object (big memory overhead)
  allmetrics$local.clustering <- transitivity(origgraph, type = "local")
  allmetrics$degree <- degree(origgraph, v = V(origgraph))
  evcent<- eigen_centrality(origgraph, directed = FALSE)
  allmetrics$eigen.cent <- evcent$vector
  allmetrics$closeness <- closeness(origgraph, mode = "all", normalized=TRUE)#, weights=1/E(origgraph)$weight)
  allmetrics$betweenness.node <- betweenness(origgraph)#, weights=1/E(origgraph)$weight)
  allmetrics$page.rank <- page.rank(origgraph, algo="prpack")$vector

  #only compute community-related statistics if the expected community attribute is available
  if (!is.null(vertex_attr(origgraph, community_attr))) {
    allmetrics$community.membership <- vertex_attr(origgraph, community_attr)
    
    ##nodal within vs between modules
    wibw <- wibw_module_degree(origgraph, community_attr=community_attr)
    allmetrics[["within.module.deg.zscore"]] <- wibw$z_within
    allmetrics[["between.module.deg.zscore"]] <- wibw$z_between
    allmetrics[["part.coeff"]] <- part_coeff(origgraph, allmetrics$community.membership)#, use.parallel = FALSE)
    allmetrics[["gateway.coeff.btw"]] <- gateway_coeff(origgraph, allmetrics$community.membership, centr = "btwn.cent")
    allmetrics[["gateway.coeff.deg"]] <- gateway_coeff(origgraph, allmetrics$community.membership, centr = "degree")    
  }

  return(allmetrics)
}

compute_nodal_metrics <- function(allg_density, ncpus=4, allowCache=TRUE) {
  require(foreach)
  require(doSNOW)
  
  expectFile <- file.path(basedir, "cache", paste0("dthreshnodemetrics_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading density-thresholded nodal statistics from file: ", expectFile)
    load(expectFile)
  } else {
    ###compute nodal metrics 
    setDefaultClusterOptions(master="localhost")
    clusterobj <- makeSOCKcluster(ncpus)
    registerDoSNOW(clusterobj)
    on.exit(try(stopCluster(clusterobj))) #shutdown cluster when function exits (either normally or crash)
    
    allmetrics.nodal <- foreach(subj=allg_density, .packages = c("igraph", "brainGraph"), .export=c("calcGraph_binary_nodal", "wibw_module_degree", "densities_desired")) %dopar% {
      #for (subj in allg_density) { #put here for more fine-grained debugging
      dl <- lapply(subj, function(dgraph) {
        glist <- calcGraph_binary_nodal(dgraph)
        glist$id <- dgraph$id #copy attributes for flattening to data.frame
        glist$density <- dgraph$density
        glist$node <- V(dgraph)$name
        return(glist)
      })
      names(dl) <- paste0("d", densities_desired)
      
      return(dl)
    }
    
    #flatten metrics down to a data.frame (assumes that each vector in the list (dens below) is of the same length
    #this should hold in general because these are nodal statistics and node number is constant
    #allmetrics.nodal is currently a [[subjects]][[densities]][[metrics]] list
    allmetrics.nodal.df <- do.call(rbind, lapply(allmetrics.nodal, function(subj) {
      do.call(rbind, lapply(subj, function(dens) {
        as.data.frame(dens) #should just be a list
      }))
    }))
    
    row.names(allmetrics.nodal.df) <- NULL #remove goofy d0.01 rownames
    save(allmetrics.nodal, allmetrics.nodal.df, file = expectFile)
  }
  
  return(list(allmetrics.nodal=allmetrics.nodal, allmetrics.nodal.df=allmetrics.nodal.df))
}

####worker bee function which compiles global graph theoretical metrics
calcGraph_binary_global <- function(origgraph, community_attr="wcomm_louvain") {
  #edgelist <- data.frame(get.edgelist(origgraph)) #not used currently
  globmetrics <- list()

  #globmetrics[["origgraph"]] <- origgraph #unclear why we should duplicate the graph object (big memory overhead)
  globmetrics[["characteristic_path_length"]] <- mean_distance(origgraph)
  globmetrics[["clustering_coefficient"]] <- transitivity(origgraph, type = "global")
  globmetrics[["small_worldness"]] <- (mean_distance(origgraph)/mean_distance(rewire(origgraph, with = keeping_degseq(loops = FALSE, niter = 1000))))/(transitivity(origgraph, type = "global")/transitivity(rewire(origgraph, with = keeping_degseq(loops = FALSE, niter = 1000)), type = "global"))
  
  if (!is.null(vertex_attr(origgraph, community_attr))) {
    globmetrics[["modularity"]] <- modularity(origgraph, vertex_attr(origgraph, community_attr))
  }
  
  return(globmetrics)
}


compute_global_metrics <- function(allg_density, ncpus=4, allowCache=TRUE) {
  require(foreach)
  require(doSNOW)
  
  expectFile <- file.path(basedir, "cache", paste0("dthreshglobmetrics_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading density-thresholded global statistics from file: ", expectFile)
    load(expectFile)
  } else {
    setDefaultClusterOptions(master="localhost")
    clusterobj <- makeSOCKcluster(ncpus)
    registerDoSNOW(clusterobj)
    on.exit(try(stopCluster(clusterobj))) #shutdown cluster when function exits (either normally or crash)
    
    ##compute global metrics 
    allmetrics.global <- foreach(subj=allg_density, .packages = c("igraph", "brainGraph"), .export=c("calcGraph_binary_global", "densities_desired")) %dopar% {
      #for (subj in allg_density) { #put here for more fine-grained debugging
      dl <- lapply(subj, function(dgraph) {
        glist <- calcGraph_binary_global(dgraph)
        glist$id <- dgraph$id #copy attributes for flattening to data.frame
        glist$density <- dgraph$density
        return(glist)
      })
      names(dl) <- paste0("d", densities_desired)
      dl
    }
    
    #flatten metrics down to a data.frame
    #allmetrics.global is currently a [[subjects]][[densities]][[metrics]] list
    allmetrics.global.df <- do.call(rbind, lapply(allmetrics.global, function(subj) {
      do.call(rbind, lapply(subj, function(dens) {
        as.data.frame(dens) #should just be a list
      }))
    }))
    
    row.names(allmetrics.global.df) <- NULL #remove goofy d0.01 rownames
    save(file=expectFile, allmetrics.global, allmetrics.global.df)
  }
  
  return(list(allmetrics.global=allmetrics.global, allmetrics.global.df=allmetrics.global.df))
}  
