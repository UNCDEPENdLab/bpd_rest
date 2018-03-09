reduce_centrality_fa <- function(nodalmetrics_df, reducemetrics, allowCache = TRUE, den = .05, browse = FALSE) {
  ##den refers to low densities, below which we want to remove from the analysis (due to instable component/factor loadings). default is .05
  
  suppressMessages(require(wle))
  suppressMessages(require(psych))
  suppressMessages(require(reshape2))
  suppressMessages(require(tidyr))
  suppressMessages(require(tidyverse))
  
  if(browse == TRUE){browser()}else{
    
    stopifnot(file.exists(file.path(basedir, "cache")))
    
    expectFile <- file.path(basedir, "cache", paste0("toanalyze.fa.thresh_", file_tag, ".RData"))
    if (file.exists(expectFile) && allowCache==TRUE) {
      message("Loading toanalyze from file: ", expectFile)
      if (thresh == "fc" && thresh_weighted =="binary"){
        output <- get(load(expectFile))
      } else {
       toanalyze <- read.csv(expectFile)
      }
    } else {
      
      
      reducemetrics <- c("eigen.cent", "betweenness.node", "leverage.cent", "degree", "page.rank", "part.coeff", "eigen.cent", "gateway.coeff.btw",  "within.module.deg")

      # stack metrics in data.frame to have subjects, nodes, densities,  --------
      #nodalmetrics_df <- get(load("/Users/nth7/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/cache/threshnodalmetrics_transformed_20rsschaefer422_nosmooth_aroma_bp_nonaggr_ridge.net_partial_fc_binary_all.RData"))
      
      metrics.raw <- nodalmetrics_df %>% dplyr::select(id, node, density, one_of(reducemetrics)) %>%
        filter(density > den) %>% gather(key="variable", value="value", -id, -node, -density)
      #recast for FA such that column names represent a metrics at a given density (e.g., 0.05_between.module.deg.zscore)
      metrics.raw_fa <- dcast(metrics.raw, id + node ~ density + variable, value.var = "value")
      toanalyze <- dplyr::select(metrics.raw_fa, id, node)
      
      
      # cor(nodalmetrics_df)
      
      
      fa.CFI<-function(x){
        nombre<-paste(x,"CFI",sep = ".")
        nombre<-
          ((x$null.chisq-x$null.dof)-(x$STATISTIC-x$dof))/(x$null.chisq-x$null.dof)
        return(nombre)
      }
      
      
      # faout_2 <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 2,  missing = TRUE, rotate = "varimax")
      # print(faout_2$loadings, cutoff = .3)
      # summary(faout_2)
      # fa.CFI(faout_2)
      # 
      faout_3 <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 3,  missing = TRUE, rotate = "promax")
      print(faout_3$loadings, cutoff = .3)
      summary(faout_3)
      fa.CFI(faout_3)

      faout_4 <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 4,  missing = TRUE, rotate = "promax")
      print(faout_4$loadings, cutoff = .3)
      summary(faout_4)
      fa.CFI(faout_4)
      # 
      # 
      faout_5 <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 5,  missing = TRUE, rotate = "promax")
      print(faout_5$loadings, cutoff = .3)
      summary(faout_5)
      fa.CFI(faout_5)

      
      #####retain for analysis, 4 factors shows imporvements in model fit with low factor loadings on factor 4
      #pcaout <- pca(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 3,  missing = TRUE, rotate = "promax")#, scores = "Bartlett")
      faout <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 3,  missing = TRUE, rotate = "promax")#, scores = "Bartlett")
      summary <- summary(faout)
      CFI <- fa.CFI(faout)
      fasolution <- data.frame(faout$scores)
      names(fasolution) <- c("central","integration", "within.mod")
      toanalyze_comb <- cbind(toanalyze, fasolution)
      
      head(toanalyze_comb)
      
      
      vss.out <- vss(dplyr::select(metrics.raw_fa, -id, -node), n = 3)
      
      CFI <- rbind(CFI, SRMR = vss.out[["vss.stats"]]$SRMR[3])
      
      
      
      if(!exists("subj_info")){subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)}
      merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)
      
      # browser()
      
      toanalyze <- dplyr::left_join(toanalyze_comb, merge.bpdage, by = "id")
      #PClab

  # what to retrn ----------------------------------------------------------


      
      if(thresh == "fc" && thresh_weighted =="binary"){
      output <- list()
      output[["toanalyze"]] <- toanalyze
      output[["faout"]] <- faout
      output[["CFI.fa"]] <- CFI
      
      save(output, file = expectFile)
      message("Saving FA-related objects to: ", expectFile)
      } else{
      write.csv(toanalyze, file = expectFile, row.names = FALSE)
      message("If you would like a list with the output of the factor analysis as well of the CFI see code under thresh == fc and thresh_weighted ==binary")
        }
    }
    
  } 
  if(thresh == "fc" && thresh_weighted =="binary"){
    return(output)} else{return(toanalyze)}
  
  } 
  
  



