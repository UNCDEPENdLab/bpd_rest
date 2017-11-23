reduce_centrality_fa <- function(nodalmetrics_df, reducemetrics, allowCache = TRUE, den = .05, weighted = FALSE, browse = FALSE){
  ##den refers to low densities, below which we want to remove from the analysis (due to instable component/factor loadings). default is .05
  
  suppressMessages(require(wle))
  suppressMessages(require(psych))
  suppressMessages(require(reshape2))
  suppressMessages(require(tidyr))
  suppressMessages(require(tidyverse))
  
  if(browse == TRUE){browser()}else{
#Binary  -----------------------------------------------------------------
  if(weighted == FALSE){
  
  stopifnot(file.exists(file.path(basedir, "cache")))
  expectFile <- file.path(basedir, "cache", paste0("toanalyze.fa.binary.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", data.reduce, ".csv"))

  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading toanalyze from file: ", expectFile)
    toanalyze <- read.csv(expectFile)
  } else {

    
    
    ####stack metrics in data.frame to have subjects, nodes, densities, *and* metrics on the rows
    metrics.raw <- nodalmetrics_df %>% dplyr::select(id, node, density, one_of(reducemetrics)) %>%
      filter(density > den) %>% gather(key="variable", value="value", -id, -node, -density)
    #recast for FA such that column names represent a metrics at a given density (e.g., 0.05_between.module.deg.zscore)
    metrics.raw_fa <- dcast(metrics.raw, id + node ~ density + variable, value.var = "value")
    toanalyze <- dplyr::select(metrics.raw_fa, id, node)

    # EFA if needed -----------------------------------------------------------
    
    # faout_3 <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 3,  missing = TRUE, rotate = "varimax")
    # print(faout_3$loadings)
    # 
    # faout_4 <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 4,  missing = TRUE, rotate = "varimax")
    # print(faout_4$loadings)
    # 
    # faout_6 <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 6,  missing = TRUE, rotate = "promax")
    # print(faout_6$loadings)
    # 
    # 
    # scree.dat <- data.frame(ev = faout_3$e.values, factor_number = seq(1:length(faout_3$e.values)))
    # ggplot(scree.dat, aes(x=factor_number, y=ev)) + geom_point() + geom_line() + theme_bw(base_size=15) +
    #   xlab("Factor number") + ylab("Eigenvalue")
    
    # Pearson -----------------------------------------------------------------
    
    
    if(conn_method == "pearson"){
      centrality.metrics <- c("degree", "eigen.cent", "page.rank")
      integration.metrics <- c("gateway.coeff.btw", "gateway.coeff.degree", "part.coeff")
      
      ###centrality
      centrality.metrics.scores <- metrics.raw[which(metrics.raw$variable %in% centrality.metrics),]  
      centrality.fa.obj <- dcast(centrality.metrics.scores, id+node~density+variable, value.var = "value" )
      faout <- fa(dplyr::select(centrality.fa.obj, -id, -node), nfactors = 1, fm = "ml", rotate = "varimax")
      print(faout$loadings)
      fasolution.central <- data.frame(faout$scores)
      names(fasolution.central) <- "central"
      fasolution.central.scores <- cbind(dplyr::select(centrality.fa.obj, id, node), fasolution.central)
      
      toanalyze <- left_join(fasolution.central.scores, toanalyze, by = c("id", "node"))
      
      ###integration
      integration.metrics.scores <- metrics.raw[which(metrics.raw$variable %in% integration.metrics),]  
      integration.fa.obj <- dcast(integration.metrics.scores, id+node~density+variable, value.var = "value" )
      
          #mean imputation
          for (i in 3:length(colnames(integration.fa.obj))){
            integration.fa.obj[which(is.nan(as.matrix(integration.fa.obj[,i]))), i] <- mean(na.omit(integration.fa.obj[,i]))
          }
          
          #any(is.nan(as.matrix(integration.fa.obj[,i]))) #check in MI worked
      
      faout <- fa(dplyr::select(integration.fa.obj, -id, -node), nfactors = 1,  fm = "ml", rotate = "varimax")
      print(faout$loadings, cutoff = 0.3)
      fasolution.integration <- data.frame(faout$scores)
      names(fasolution.integration) <- "integration"
      fasolution.integration.scores <- cbind(dplyr::select(integration.fa.obj, id, node), fasolution.integration)
      
      toanalyze <- left_join(toanalyze, fasolution.integration.scores, by = c("id", "node"))
      
      ###within module connectivity
      within.mod.scores <- metrics.raw[which(metrics.raw$variable == "within.module.deg"),]  
      within.mod.fa.obj <- dcast(within.mod.scores, id+node~density+variable, value.var = "value" )
      faout <- fa(dplyr::select(within.mod.fa.obj, -id, -node), nfactors = 1, fm = "ml", rotate = "varimax")
      print(faout$loadings, cutoff = 0.3)
      fasolution.within.mod <- data.frame(faout$scores)
      names(fasolution.within.mod) <- "within.mod"
      fasolution.within.mod.scores <- cbind(dplyr::select(within.mod.fa.obj, id, node), fasolution.within.mod)
      
      toanalyze <- left_join(toanalyze, fasolution.within.mod.scores, by = c("id", "node"))
      
      ##closeness 
      # closeness.scores <- metrics.raw[which(metrics.raw$variable == "closeness"),]  
      # closeness.fa.obj <- dcast(closeness.scores, id+node~density+variable, value.var = "value" )
      # faout <- fa(dplyr::select(closeness.fa.obj, -id, -node), nfactors = 1, fm = "ml", rotate = "varimax")
      # print(faout$loadings, cutoff = 0.3)
      # fasolution.closeness <- data.frame(faout$scores)
      # names(fasolution.closeness) <- "closeness"
      # fasolution.closeness.scores <- cbind(dplyr::select(closeness.fa.obj, id, node), fasolution.closeness)
      # 
      # toanalyze <- left_join(toanalyze, fasolution.closeness.scores, by = c("id", "node"))
      # 
      ##betweenness centrality
      btw.scores <- metrics.raw[which(metrics.raw$variable == "betweenness.node"),]  
      btw.fa.obj <- dcast(btw.scores, id+node~density+variable, value.var = "value" )
      faout <- fa(dplyr::select(btw.fa.obj, -id, -node), nfactors = 1, fm = "ml", rotate = "varimax")
      print(faout$loadings, cutoff = 0.3)
      fasolution.btw <- data.frame(faout$scores)
      names(fasolution.btw) <- "betweenness"
      fasolution.btw.scores <- cbind(dplyr::select(btw.fa.obj, id, node), fasolution.btw)
      
      toanalyze <- left_join(toanalyze, fasolution.btw.scores, by = c("id", "node"))
      
      merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)
      
      toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")
      #PClabs <- c("PC1: Centrality Aggregate", "PC2: Integration", "PC3: Within Module Connectivity", "PC4: Closeness Centrality", "PC5: Betweenness Centrality")
      
    }
    
    # Ridge-----------------------------------------------------------------
    
    if(conn_method == "ridge.net_partial"){
     
      centrality.metrics <- c("degree", "eigen.cent", "page.rank", "betweenness.node")
      integration.metrics <- c("gateway.coeff.btw", "gateway.coeff.degree", "part.coeff")
      # browser()
      ###centrality
      centrality.metrics.scores <- metrics.raw[which(metrics.raw$variable %in% centrality.metrics),]  
      centrality.fa.obj <- dcast(centrality.metrics.scores, id+node~density+variable, value.var = "value" )
      faout <- fa(dplyr::select(centrality.fa.obj, -id, -node), nfactors = 1, fm = "ml", rotate = "varimax")
      print(faout$loadings)
      fasolution.central <- data.frame(faout$scores)
      names(fasolution.central) <- "central"
      fasolution.central.scores <- cbind(dplyr::select(centrality.fa.obj, id, node), fasolution.central)
      
      toanalyze <- left_join(fasolution.central.scores, toanalyze, by = c("id", "node"))
      
      ###integration
      integration.metrics.scores <- metrics.raw[which(metrics.raw$variable %in% integration.metrics),]  
      integration.fa.obj <- dcast(integration.metrics.scores, id+node~density+variable, value.var = "value" )
      
          #mean imputation
          for (i in 3:length(colnames(integration.fa.obj))){
            integration.fa.obj[which(is.nan(as.matrix(integration.fa.obj[,i]))), i] <- mean(na.omit(integration.fa.obj[,i]))
          }
      
            #any(is.nan(as.matrix(integration.fa.obj[,i]))) #check in MI worked
      
      faout <- fa(dplyr::select(integration.fa.obj, -id, -node), nfactors = 1, fm = "ml", rotate = "varimax")
      print(faout$loadings, cutoff = 0.3)
      fasolution.integration <- data.frame(faout$scores)
      names(fasolution.integration) <- "integration"
      fasolution.integration.scores <- cbind(dplyr::select(integration.fa.obj, id, node), fasolution.integration)
      
      toanalyze <- left_join(toanalyze, fasolution.integration.scores, by = c("id", "node"))
      
      ###within module connectivity
      within.mod.scores <- metrics.raw[which(metrics.raw$variable == "within.module.deg"),]  
      within.mod.fa.obj <- dcast(within.mod.scores, id+node~density+variable, value.var = "value" )
      faout <- fa(dplyr::select(within.mod.fa.obj, -id, -node), nfactors = 1, fm = "ml", rotate = "varimax")
      print(faout$loadings, cutoff = 0.3)
      fasolution.within.mod <- data.frame(faout$scores)
      names(fasolution.within.mod) <- "within.mod"
      fasolution.within.mod.scores <- cbind(dplyr::select(within.mod.fa.obj, id, node), fasolution.within.mod)
      
      toanalyze <- left_join(toanalyze, fasolution.within.mod.scores, by = c("id", "node"))
      
      
      merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)
      
      toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")
      #PClabs <- c("PC1: Centrality Aggregate", "PC2: Integration", "PC3: Within Module Connectivity", "PC4: Closeness Centrality")
      }


    write.csv(toanalyze, file = paste0(basedir, "/cache/toanalyze.fa.binary.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", data.reduce, ".csv"), row.names = FALSE)
    }
  } 
# Weighted IN PROGRESS 11/11/17----------------------------------------------------------------
    else{
 
  message("WEIGHTED FA NOT CURRENTLY SUPPORTED")
  
  stopifnot(file.exists(file.path(basedir, "cache")))  
  expectFile <- file.path(basedir, "cache", paste0("toanalyze.fa.weighted.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", data.reduce, ".csv"))
  
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading toanalyze from file: ", expectFile)
    toanalyze <- read.csv(expectFile)
  } else {
    
    ####stack metrics in data.frame to have subjects, nodes *and* metrics on the rows
  
    metrics.raw <- nodalmetrics_df %>% dplyr::select(id, node, one_of(reducemetrics)) %>% gather(key="variable", value="value", -id, -node)#, -density)
    metrics.raw <- dcast(metrics.raw, id + node ~  variable, value.var = "value")
    
    toanalyze <- dplyr::select(metrics.raw, id, node)
    
    metrics.raw <- metrics.raw %>% dplyr::select(-id, -node)#, -closeness, -page.rank)
     
    # if(conn_method == "pearson"){
    #   pcaout5 <- pca(dplyr::select(metrics.raw_pca, -id, -node), nfactors=5, rotate="varimax")
    #   print(pcaout5$loadings, cutoff = 0.3)
    #   pcasolution <- data.frame(pcaout5$scores) #use these, which are the rotated scores, not the raw eigenvectors from $x in prcomp
    #   #names(pcasolution) <- c("central", "between.node", "within.mod") #relevant for stats without gateway and participation
    #   names(pcasolution) <- c("integration", "central", "within.mod", "closeness", "betweenness") 
    #   
    #   
    #   toanalyze <- cbind(dplyr::select(metrics.raw_pca, id, node), pcasolution)
    #   #subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz")
    #   merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)
    #   
    #   toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")
    # }
    
    if(conn_method == "ridge.net_partial"){
      
      faout3 <- fa(metrics.raw, nfactors = 3, fm = "ml", rotate = "promax")
      print(faout3$loadings)
      
      faout4 <- fa(metrics.raw, nfactors = 4, fm = "ml", rotate = "promax")
      print(faout4$loadings)
      
      faout5 <- fa(metrics.raw, nfactors = 5, fm = "ml", rotate = "varimax")
      print(faout5$loadings)
      
      pcaout4 <- pca(dplyr::select(metrics.raw_pca, -id, -node), nfactors=4, rotate="varimax")
      print(pcaout4$loadings, cutoff = 0.3)
      pcasolution <- data.frame(pcaout4$scores) #use these, which are the rotated scores, not the raw eigenvectors from $x in prcomp
      #names(pcasolution) <- c("central", "between.node", "within.mod") #relevant for stats without gateway and participation
      names(pcasolution) <- c("central", "integration", "within.mod", "closeness") 
      
      
      toanalyze <- cbind(dplyr::select(metrics.raw_pca, id, node), pcasolution)
      #subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz")
      merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)
      
      toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")
      }
    
    write.csv(toanalyze, file = paste0(basedir, "/cache/toanalyze.fa.weighted.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", data.reduce, ".csv"), row.names = FALSE)
      }
    }
  }
  return(toanalyze)
}



