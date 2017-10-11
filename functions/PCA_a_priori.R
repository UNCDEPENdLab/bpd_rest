PCA_a_priori <- function(nodalmetrics_dthresh_df, pcametrics, allowCache = TRUE, den = .04){
  ##den refers to low densities, below which we want to remove from the PCA analysis. default is .04 
  
  require(wle)  
  require(psych)
  require(reshape2)
  require(tidyr)
  
  if(is.null(nodalmetrics_dthresh_df)){
    #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame
    nodalmetrics_dthresh_df <- load_nodal_metrics_df()
  }
  
  stopifnot(file.exists(file.path(basedir, "cache")))  
  expectFile.toanalyze <- file.path(basedir, "cache", paste0("toanalyze.pca.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method, ".csv"))
  expectFile.pclabs <- file.path(basedir, "cache", paste0("PClabs.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method, ".R"))
  
  
   if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading toanalyze from file: ", expectFile.toanalyze)
    toanalyze <- read.csv(expectFile.toanalyze)
    message("Loading PClabs from file:", expectFile.pclabs)
    load(expectFile.pclabs)
  } else {

  ####stack metrics in data.frame to have subjects, nodes, densities, *and* metrics on the rows
  metrics.raw <- nodalmetrics_dthresh_df %>% select(id, node, density, one_of(pcametrics)) %>%
    filter(density > den) %>% gather(key="variable", value="value", -id, -node, -density)
  
  #recast for PCA such that column names represent a metrics at a given density (e.g., 0.05_between.module.deg.zscore)
  metrics.raw_pca <- dcast(metrics.raw, id + node ~ density + variable, value.var = "value")
  
  toanalyze <- select(metrics.raw_pca, id, node)
  
  ################create a priori groupings based on past analyses using PCA_all
  if(conn_method == "pearson"){
    centrality.metrics <- c("degree", "eigen.cent", "page.rank")
    integration.metrics <- c("gateway.coeff.btw", "gateway.coeff.deg", "part.coeff")
    
    ###centrality
    centrality.metrics.scores <- metrics.raw[which(metrics.raw$variable %in% centrality.metrics),]  
    centrality.pca.obj <- dcast(centrality.metrics.scores, id+node~density+variable, value.var = "value" )
    pcaout <- pca(select(centrality.pca.obj, -id, -node), nfactors = 1, rotate = "varimax")
    print(pcaout$loadings)
    pcasolution.central <- data.frame(pcaout$scores)
    names(pcasolution.central) <- "central"
    pcasolution.central.scores <- cbind(select(centrality.pca.obj, id, node), pcasolution.central)
    
    toanalyze <- left_join(pcasolution.central.scores, toanalyze, by = c("id", "node"))
    
    ###integration
    integration.metrics.scores <- metrics.raw[which(metrics.raw$variable %in% integration.metrics),]  
    integration.pca.obj <- dcast(integration.metrics.scores, id+node~density+variable, value.var = "value" )
    pcaout <- pca(select(integration.pca.obj, -id, -node), nfactors = 1, rotate = "varimax")
    print(pcaout$loadings, cutoff = 0.3)
    pcasolution.integration <- data.frame(pcaout$scores)
    names(pcasolution.integration) <- "integration"
    pcasolution.integration.scores <- cbind(select(integration.pca.obj, id, node), pcasolution.integration)
    
    toanalyze <- left_join(toanalyze, pcasolution.integration.scores, by = c("id", "node"))
    
    ###within module connectivity
    within.mod.scores <- metrics.raw[which(metrics.raw$variable == "within.module.deg.zscore"),]  
    within.mod.pca.obj <- dcast(within.mod.scores, id+node~density+variable, value.var = "value" )
    pcaout <- pca(select(within.mod.pca.obj, -id, -node), nfactors = 1, rotate = "varimax")
    print(pcaout$loadings, cutoff = 0.3)
    pcasolution.within.mod <- data.frame(pcaout$scores)
    names(pcasolution.within.mod) <- "within.mod"
    pcasolution.within.mod.scores <- cbind(select(within.mod.pca.obj, id, node), pcasolution.within.mod)
    
    toanalyze <- left_join(toanalyze, pcasolution.within.mod.scores, by = c("id", "node"))
    
    ##closeness 
    closeness.scores <- metrics.raw[which(metrics.raw$variable == "closeness"),]  
    closeness.pca.obj <- dcast(closeness.scores, id+node~density+variable, value.var = "value" )
    pcaout <- pca(select(closeness.pca.obj, -id, -node), nfactors = 1, rotate = "varimax")
    print(pcaout$loadings, cutoff = 0.3)
    pcasolution.closeness <- data.frame(pcaout$scores)
    names(pcasolution.closeness) <- "closeness"
    pcasolution.closeness.scores <- cbind(select(closeness.pca.obj, id, node), pcasolution.closeness)
    
    toanalyze <- left_join(toanalyze, pcasolution.closeness.scores, by = c("id", "node"))
    
    ##betweenness centrality
    btw.scores <- metrics.raw[which(metrics.raw$variable == "betweenness.node"),]  
    btw.pca.obj <- dcast(btw.scores, id+node~density+variable, value.var = "value" )
    pcaout <- pca(select(btw.pca.obj, -id, -node), nfactors = 1, rotate = "varimax")
    print(pcaout$loadings, cutoff = 0.3)
    pcasolution.btw <- data.frame(pcaout$scores)
    names(pcasolution.btw) <- "betweenness"
    pcasolution.btw.scores <- cbind(select(btw.pca.obj, id, node), pcasolution.btw)
    
    toanalyze <- left_join(toanalyze, pcasolution.btw.scores, by = c("id", "node"))
    
    merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)
    
    toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")
    PClabs <- c("PC1: Centrality Aggregate", "PC2: Integration", "PC3: Within Module Connectivity", "PC4: Closeness Centrality", "PC5: Betweenness Centrality")
    
  }
  
  if(conn_method == "ridge.net_partial"){
    centrality.metrics <- c("degree", "eigen.cent", "betweenness.node", "page.rank")
    integration.metrics <- c("gateway.coeff.btw", "gateway.coeff.deg", "part.coeff")
    
    ###centrality
    centrality.metrics.scores <- metrics.raw[which(metrics.raw$variable %in% centrality.metrics),]  
    centrality.pca.obj <- dcast(centrality.metrics.scores, id+node~density+variable, value.var = "value" )
    pcaout <- pca(select(centrality.pca.obj, -id, -node), nfactors = 1, rotate = "varimax")
    print(pcaout$loadings)
    pcasolution.central <- data.frame(pcaout$scores)
    names(pcasolution.central) <- "central"
    pcasolution.central.scores <- cbind(select(centrality.pca.obj, id, node), pcasolution.central)
    
    toanalyze <- left_join(pcasolution.central.scores, toanalyze, by = c("id", "node"))
    
    ###integration
    integration.metrics.scores <- metrics.raw[which(metrics.raw$variable %in% integration.metrics),]  
    integration.pca.obj <- dcast(integration.metrics.scores, id+node~density+variable, value.var = "value" )
    pcaout <- pca(select(integration.pca.obj, -id, -node), nfactors = 1, rotate = "varimax")
    print(pcaout$loadings, cutoff = 0.3)
    pcasolution.integration <- data.frame(pcaout$scores)
    names(pcasolution.integration) <- "integration"
    pcasolution.integration.scores <- cbind(select(integration.pca.obj, id, node), pcasolution.integration)
    
    toanalyze <- left_join(toanalyze, pcasolution.integration.scores, by = c("id", "node"))
    
    ###within module connectivity
    within.mod.scores <- metrics.raw[which(metrics.raw$variable == "within.module.deg.zscore"),]  
    within.mod.pca.obj <- dcast(within.mod.scores, id+node~density+variable, value.var = "value" )
    pcaout <- pca(select(within.mod.pca.obj, -id, -node), nfactors = 1, rotate = "varimax")
    print(pcaout$loadings, cutoff = 0.3)
    pcasolution.within.mod <- data.frame(pcaout$scores)
    names(pcasolution.within.mod) <- "within.mod"
    pcasolution.within.mod.scores <- cbind(select(within.mod.pca.obj, id, node), pcasolution.within.mod)
    
    toanalyze <- left_join(toanalyze, pcasolution.within.mod.scores, by = c("id", "node"))
    
    ##closeness 
    closeness.scores <- metrics.raw[which(metrics.raw$variable == "closeness"),]  
    closeness.pca.obj <- dcast(closeness.scores, id+node~density+variable, value.var = "value" )
    pcaout <- pca(select(closeness.pca.obj, -id, -node), nfactors = 1, rotate = "varimax")
    print(pcaout$loadings, cutoff = 0.3)
    pcasolution.closeness <- data.frame(pcaout$scores)
    names(pcasolution.closeness) <- "closeness"
    pcasolution.closeness.scores <- cbind(select(closeness.pca.obj, id, node), pcasolution.closeness)
    
    toanalyze <- left_join(toanalyze, pcasolution.closeness.scores, by = c("id", "node"))
    
    merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)
    
    toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")
    PClabs <- c("PC1: Centrality Aggregate", "PC2: Integration", "PC3: Within Module Connectivity", "PC4: Closeness Centrality")
    
    }
  
  write.csv(toanalyze, file = paste0(basedir, "/cache/toanalyze.pca.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method, ".csv"), row.names = FALSE)
  save(PClabs, file = paste0(basedir, "/cache/PClabs.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method, ".R"))
  }
  
  metrics.pca <- select(toanalyze, -id, -node, -BPD, -Age)
  metrics.pca <- colnames(metrics.pca) 
  
  returnlist <- list(toanalyze, PClabs, metrics.pca)
   return(returnlist)
}
