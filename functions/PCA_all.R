PCA_all <- function(nodalmetrics_dthresh_df, pcametrics, allowCache = TRUE, den = .04){
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
  expectFile <- file.path(basedir, "cache", paste0("toanalyze.pca.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method, ".csv"))
  
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading toanalyze from file: ", expectFile)
    toanalyze <- read.csv(expectFile)
  } else {
    
    ####stack metrics in data.frame to have subjects, nodes, densities, *and* metrics on the rows
    metrics.raw <- nodalmetrics_dthresh_df %>% select(id, node, density, one_of(pcametrics)) %>%
      filter(density > den) %>% gather(key="variable", value="value", -id, -node, -density)
    
    #recast for PCA such that column names represent a metrics at a given density (e.g., 0.05_between.module.deg.zscore)
    metrics.raw_pca <- dcast(metrics.raw, id + node ~ density + variable, value.var = "value")
    
    toanalyze <- select(metrics.raw_pca, id, node)
    
    
    if(conn_method == "pearson"){
      pcaout5 <- pca(select(metrics.raw_pca, -id, -node), nfactors=5, rotate="varimax")
      print(pcaout5$loadings, cutoff = 0.3)
      pcasolution <- data.frame(pcaout5$scores) #use these, which are the rotated scores, not the raw eigenvectors from $x in prcomp
      #names(pcasolution) <- c("central", "between.node", "within.mod") #relevant for stats without gateway and participation
      names(pcasolution) <- c("integration", "central", "within.mod", "closeness", "betweenness") 
     
      
      toanalyze <- cbind(select(metrics.raw_pca, id, node), pcasolution)
      #subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz")
      merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)
      
      toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")
    }
    
    if(conn_method == "ridge.net_partial"){
      pcaout4 <- pca(select(metrics.raw_pca, -id, -node), nfactors=4, rotate="varimax")
      print(pcaout4$loadings, cutoff = 0.3)
      pcasolution <- data.frame(pcaout4$scores) #use these, which are the rotated scores, not the raw eigenvectors from $x in prcomp
      #names(pcasolution) <- c("central", "between.node", "within.mod") #relevant for stats without gateway and participation
      names(pcasolution) <- c("central", "integration", "within.mod", "closeness") 
      
      
      toanalyze <- cbind(select(metrics.raw_pca, id, node), pcasolution)
      #subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz")
      merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)
      
      toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")
    }
    
    write.csv(toanalyze, file = paste0(basedir, "/cache/toanalyze.pca.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method, ".csv"), row.names = FALSE)
  }
  return(toanalyze)
}
