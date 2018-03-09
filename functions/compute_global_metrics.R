compute_global_metrics <- function(graphs, ncpus=4, allowCache=TRUE, community_attr="community", weighted = 0) {
  require(foreach)
  require(doSNOW)
  if(weighted == 1){
    ####################################
    ##BEGIN WEIGHTED (no threshold)
    expectFile <- file.path(basedir, "cache", paste0("globmetrics_", file_tag_nothresh, ".RData"))
    if (file.exists(expectFile) && allowCache==TRUE) {
      message("Loading weighted global statistics from file: ", expectFile)
      load(expectFile)
    } else {
      setDefaultClusterOptions(master="localhost")
      clusterobj <- makeSOCKcluster(ncpus)
      registerDoSNOW(clusterobj)
      on.exit(try(stopCluster(clusterobj))) #shutdown cluster when function exits (either normally or crash)
      
      ##compute global metrics
      globmetrics_weighted <- foreach(subj=graphs, .packages = c("igraph", "brainGraph"), .export=c("calcGraph_global")) %dopar% {
        glist <- calcGraph_global(subj, community_attr = community_attr)
      }
      #flatten metrics down to a data.frame
      #allmetrics.global is currently a [[subjects]][[densities]][[metrics]] list
      allmetrics.global.df <- do.call(rbind, globmetrics_weighted)
      save(file=expectFile, globmetrics_weighted, allmetrics.global.df)
    }
    return(list(allmetrics.global=globmetrics_weighted, allmetrics.global.df=allmetrics.global.df))
  } else {
  ####################################
  ##BEGIN Thresholded global stats
  
  expectFile <- file.path(basedir, "cache", paste0("threshglobmetrics_", file_tag, ".RData"))
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading thresholded global statistics from file: ", expectFile)
    load(expectFile)
  } else {
    setDefaultClusterOptions(master="localhost")
    clusterobj <- makeSOCKcluster(ncpus)
    registerDoSNOW(clusterobj)
    on.exit(try(stopCluster(clusterobj))) #shutdown cluster when function exits (either normally or crash)
    # browser()
    ##compute global metrics 
    allmetrics.global <- foreach(subj=graphs, .packages = c("igraph", "brainGraph"), .export=c("calcGraph_global", "densities_desired", "rs_desired_log", "thresh")) %dopar% {
      #for (subj in graphs) { #put here for more fine-grained debugging
      dl <- lapply(subj, function(dgraph) {
        glist <- calcGraph_global(dgraph, community_attr=community_attr) #this will work for both weighted and unweighted graphs, right now modularity weighted components set to NULL but can change if desired. 
        glist$id <- dgraph$id #copy attributes for flattening to data.frame
        glist$density <- dgraph$density
        return(glist)
      })
      
      ifelse(thresh == "fc", names(dl) <- paste0("r", rs_desired_log), names(dl) <- paste0("d", densities_desired))
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
}  
