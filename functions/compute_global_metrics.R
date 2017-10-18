compute_global_metrics <- function(graphs, ncpus=4, allowCache=TRUE, community_attr="community") {
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
    allmetrics.global <- foreach(subj=graphs, .packages = c("igraph", "brainGraph"), .export=c("calcGraph_global", "densities_desired")) %dopar% {
      #for (subj in graphs) { #put here for more fine-grained debugging
      dl <- lapply(subj, function(dgraph) {
        glist <- calcGraph_global(dgraph, community_attr=community_attr) #this will work for both weighted and unweighted graphs, right now modularity weighted components set to NULL but can change if desired. 
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
