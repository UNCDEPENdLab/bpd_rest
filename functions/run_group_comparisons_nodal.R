run_group_comparisons_nodal <- function(toanalyze, abst = 2.64, allowCache = TRUE, weighted = FALSE, browse = FALSE){
  if(browse == TRUE) {browser()}
  if (weighted == FALSE){

# Binary Graphs -----------------------------------------------------------

    
  stopifnot(file.exists(file.path(basedir, "results")))  
  expectFile.ttest <- file.path(basedir, "results", paste0(parcellation, "_", conn_method),   paste0("all.sigttest.nodal.binary_", data.reduce, ".csv"))
  expectFile.wle <- file.path(basedir, "results", paste0(parcellation, "_", conn_method), paste0("all.sigwlelm.nodal.binary_", data.reduce, ".csv"))
  expectFile.lm <- file.path(basedir, "results", paste0(parcellation, "_", conn_method), paste0("all.siglm.nodal.binary_", data.reduce, ".csv"))
  
  if (file.exists(expectFile.ttest) && file.exists(expectFile.lm) && file.exists(expectFile.wle) && allowCache==TRUE) {
    message("Loading binary nodal metric comparisons from: ", expectFile.ttest, ", ", expectFile.lm, ", ", expectFile.wle)
    all.sigwlelm.nodal <- read.csv(expectFile.wle)
    all.siglm.nodal <- read.csv(expectFile.lm) 
    all.sigttest.nodal <- read.csv(expectFile.ttest)
  } else {
  
  metrics.toanalyze <- names(dplyr::select(toanalyze, -id, -node, -BPD, -Age))
  ttest_bpd <- list()
  lm_eff <- list() #all effects from lm
  wle.lm_eff <- list() #all effects from weighted likelihood (robust estimation)
  
  for (metric in metrics.toanalyze) {
    for (node in levels(toanalyze$node)) {
      thismetric <- toanalyze[toanalyze$node == node, c(metric, "BPD", "Age")]
      thismetric$BPD <- factor(thismetric$BPD, levels = c(0,1), labels = c("control", "BPD"))
      colnames(thismetric) <- c("metric", "group", "age")
      
      node.test <- tryCatch(t.test(metric~group, thismetric, na.rm = TRUE), error = function(errorname) { print(errorname); return(NULL) })
      if (is.null(node.test)) { message("Error occurred for metric: ", metric, " and node: ", node) }
      
      if ((!is.null(node.test) && !is.nan(node.test$statistic) && abs(node.test$statistic) > abst)) {
        ttest_bpd[[length(ttest_bpd)+1]] <- data.frame(broom::tidy(node.test), nodename=as.character(atlas$anat_label[atlas$name == node]), metric=metric, nodenum=node)
      }
      
      
      age.test.wle <- tryCatch(wle.lm(metric~age*group, thismetric), error=function(e) { print(e); return(NULL) })
      if (!is.null(age.test.wle) && any(sigp <- summary(age.test.wle)[[1]]$coefficients[-1,"Pr(>|t|)"] < .005)) { #summary returns a top-level list. Use [[1]] to access coefficients
        df <- data.frame(summary(age.test.wle)[[1]]$coefficients[c(FALSE, sigp),,drop=FALSE]) #use sigp to just pull sig rows. FALSE specifies to skip intercept since this was omitted in check above. drop=FALSE forces it to stay a data.frame
        df$term <- row.names(df); row.names(df) <- NULL #manual tidying... a bit prone to explosion if columns change
        names(df) <- c("estimate", "std.error", "statistic", "p.value", "term")
        df <- df[,c("term", "estimate", "std.error", "statistic", "p.value")] #rearrange for consistency
        
        wle.lm_eff[[length(wle.lm_eff) + 1]] <- data.frame(df, nodenum=node, nodename=as.character(atlas$anat_label[atlas$name == node]), metric = metric)
      }
      
      age.test <- tryCatch(lm(metric~age*group, thismetric, na.action = "na.exclude"), error=function(e) { print(e); return(NULL) })
      if (is.null(age.test)) { pvec <- NULL
      } else { 
        agedf <- broom::tidy(age.test)
        pvec <- agedf$p.value[-1] #p-values of age, bpd, and age x bpd. -1 to drop off intercept
        names(pvec) <- agedf$term[-1]
      }
      
      if (!is.null(age.test) && !all(is.nan(pvec)) && pvec["groupBPD"] < .005) {
        lm_eff[[length(lm_eff) + 1]] <- data.frame(subset(agedf, term=="groupBPD"), nodenum=node, nodename=as.character(atlas$anat_label[atlas$name == node]), metric = metric)
      }
      
      if (!is.null(age.test) && !all(is.nan(pvec)) && pvec["age"] < .005) {
        lm_eff[[length(lm_eff) + 1]] <- data.frame(subset(agedf, term=="age"), nodenum=node, nodename=as.character(atlas$anat_label[atlas$name == node]), metric = metric)
      }
      
      if (!is.null(age.test) && !all(is.nan(pvec)) && pvec["age:groupBPD"] < .005) {
        lm_eff[[length(lm_eff) + 1]] <- data.frame(subset(agedf, term=="age:groupBPD"), nodenum=node, nodename=as.character(atlas$anat_label[atlas$name == node]), metric = metric)
      }
    }
  }
  
  #results.df <- do.call(rbind, results)
  
  ####BPD main effects
  results.ttestbpd <- dplyr::bind_rows(ttest_bpd) #efficient version of do.call(rbind, dflist) for list of data.frames
  results.lm <- dplyr::bind_rows(lm_eff)
  results.wle.lm <- dplyr::bind_rows(wle.lm_eff)
  
  
  if(use.yeo == 1){ 
    #message("Compiling communities based on Yeo 7 Networks")
    community <- yeo7_community(agg.g)
    
    membership <- data.frame(nodenum=names(community$membership), community = community$membership, 
                             community.name=mapvalues(community$membership, from = c("1","2","3","4","5","6","7"), to = c("VIS", "SOMMOT", "DORSATTN", "SALVENTATTN", "LIMBIC", "FPN", "DMN")),
                             stringsAsFactors = FALSE)
  } else {
    message("Compiling communities based on Louvain fast-unfolding at d_12 for N=83")
    community <- readRDS(paste0(getwd(), "/cache/power269/d12_louv_n83.rds"))
    
    membership <- data.frame(nodenum=names(community$membership), community = community$membership, 
                             community.name=mapvalues(community$membership, from = c("1","2","3","4","5"), to = c("SOMMOTOR", "FPN/DA", "OCC", "CO", "DMN")),
                             stringsAsFactors = FALSE)
  }
  
  
  all.siglm.nodal <- dplyr::left_join(results.lm, membership, by = "nodenum") %>% arrange(term, metric, nodenum)
  all.sigwlelm.nodal <- dplyr::left_join(results.wle.lm, membership, by = "nodenum") %>% arrange(term, metric, nodenum)
  all.sigttest.nodal <- dplyr::left_join(results.ttestbpd, membership, by = "nodenum") %>% arrange(metric, nodenum)
  
  write.csv(all.sigwlelm.nodal, file = paste0(basedir, "/results/", parcellation, "_", conn_method,  "/all.sigwlelm.nodal.binary_", data.reduce, ".csv"), row.names=FALSE)
  write.csv(all.siglm.nodal, file = paste0(basedir, "/results/", parcellation, "_", conn_method, "/all.siglm.nodal.binary_", data.reduce, ".csv"), row.names=FALSE)
  write.csv(all.sigttest.nodal, file = paste0(basedir, "/results/", parcellation, "_", conn_method,  "/all.sigttest.nodal.binary_", data.reduce,".csv"), row.names=FALSE)
  
  }
  
  return(list(all.sigwlelm.nodal=all.sigwlelm.nodal,
              all.siglm.nodal=all.siglm.nodal, 
              all.sigttest.nodal=all.sigttest.nodal))
  } else{
# Weighted graphs ---------------------------------------------------------

    stopifnot(file.exists(file.path(basedir, "results")))  
    expectFile.ttest <- file.path(basedir, "results", paste0(parcellation, "_", conn_method),   paste0("all.sigttest.nodal.weighted_", data.reduce, ".csv"))
    expectFile.wle <- file.path(basedir, "results", paste0(parcellation, "_", conn_method), paste0("all.sigwlelm.nodal.weighted_", data.reduce, ".csv"))
    expectFile.lm <- file.path(basedir, "results", paste0(parcellation, "_", conn_method), paste0("all.siglm.nodal.weighted_", data.reduce, ".csv"))
    
    if (file.exists(expectFile.ttest) && file.exists(expectFile.lm) && file.exists(expectFile.wle) && allowCache==TRUE) {
      message("Loading weighted nodal metric comparisons from: ", expectFile.ttest, ", ", expectFile.lm, ", ", expectFile.wle)
      all.sigwlelm.nodal <- read.csv(expectFile.wle)
      all.siglm.nodal <- read.csv(expectFile.lm) 
      all.sigttest.nodal <- read.csv(expectFile.ttest)
    } else {
      
      metrics.toanalyze <- names(dplyr::select(toanalyze, -id, -node, -BPD, -Age, -community.membership))
      ttest_bpd <- list()
      lm_eff <- list() #all effects from lm
      wle.lm_eff <- list() #all effects from weighted likelihood (robust estimation)
      
      # browser()
      for (metric in metrics.toanalyze) {
        for (node in levels(toanalyze$node)) {
          
          thismetric <- toanalyze[toanalyze$node == node, c(metric, "BPD", "Age")]
          thismetric$BPD <- factor(thismetric$BPD, levels = c(0,1), labels = c("control", "BPD"))
          colnames(thismetric) <- c("metric", "group", "age")
          
          node.test <- tryCatch(t.test(metric~group, thismetric, na.rm = TRUE), error = function(errorname) { print(errorname); return(NULL) })
          if (is.null(node.test)) { message("Error occurred for metric: ", metric, " and node: ", node) }
          
          if ((!is.null(node.test) && !is.nan(node.test$statistic) && abs(node.test$statistic) > abst)) {
            ttest_bpd[[length(ttest_bpd)+1]] <- data.frame(broom::tidy(node.test), nodename=as.character(atlas$anat_label[atlas$name == node]), metric=metric, nodenum=node)
          }
          
          
          age.test.wle <- tryCatch(wle.lm(metric~age*group, thismetric), error=function(e) { print(e); return(NULL) })
          if (!is.null(age.test.wle) &&!all(is.na(age.test.wle$coefficients)) && any(sigp <- summary(age.test.wle)[[1]]$coefficients[-1,"Pr(>|t|)"] < .005)) { #summary returns a top-level list. Use [[1]] to access coefficients
            
            df <- data.frame(summary(age.test.wle)[[1]]$coefficients[c(FALSE, sigp),,drop=FALSE]) #use sigp to just pull sig rows. FALSE specifies to skip intercept since this was omitted in check above. drop=FALSE forces it to stay a data.frame
            df$term <- row.names(df); row.names(df) <- NULL #manual tidying... a bit prone to explosion if columns change
            names(df) <- c("estimate", "std.error", "statistic", "p.value", "term")
            df <- df[,c("term", "estimate", "std.error", "statistic", "p.value")] #rearrange for consistency

            wle.lm_eff[[length(wle.lm_eff) + 1]] <- data.frame(df, nodenum=node, nodename=as.character(atlas$anat_label[atlas$name == node]), metric = metric)
          }

          age.test <- tryCatch(lm(metric~age*group, thismetric, na.action = "na.exclude"), error=function(e) { print(e); return(NULL) })
          if (is.null(age.test)) { pvec <- NULL
          } else { 
            agedf <- broom::tidy(age.test)
            pvec <- agedf$p.value[-1] #p-values of age, bpd, and age x bpd. -1 to drop off intercept
            names(pvec) <- agedf$term[-1]
          }
          
          if (!is.null(age.test) && !all(is.nan(pvec)) && pvec["groupBPD"] < .005) {
            lm_eff[[length(lm_eff) + 1]] <- data.frame(subset(agedf, term=="groupBPD"), nodenum=node, nodename=as.character(atlas$anat_label[atlas$name == node]), metric = metric)
          }
          
          if (!is.null(age.test) && !all(is.nan(pvec)) && pvec["age"] < .005) {
            lm_eff[[length(lm_eff) + 1]] <- data.frame(subset(agedf, term=="age"), nodenum=node, nodename=as.character(atlas$anat_label[atlas$name == node]), metric = metric)
          }
          
          if (!is.null(age.test) && !all(is.nan(pvec)) && pvec["age:groupBPD"] < .005) {
            lm_eff[[length(lm_eff) + 1]] <- data.frame(subset(agedf, term=="age:groupBPD"), nodenum=node, nodename=as.character(atlas$anat_label[atlas$name == node]), metric = metric)
          }
        }
      }
      
     
      #results.df <- do.call(rbind, results)
      
      ####BPD main effects
      results.ttestbpd <- dplyr::bind_rows(ttest_bpd) #efficient version of do.call(rbind, dflist) for list of data.frames
      results.lm <- dplyr::bind_rows(lm_eff)
      results.wle.lm <- dplyr::bind_rows(wle.lm_eff)
      
      
      if(use.yeo == 1){ 
        message("Compiling communities based on Yeo 7 Networks")
        community <- yeo7_community(agg.g)
        
        membership <- data.frame(nodenum=names(community$membership), community = community$membership, 
                                 community.name=mapvalues(community$membership, from = c("1","2","3","4","5","6","7"), to = c("VIS", "SOMMOT", "DORSATTN", "SALVENTATTN", "LIMBIC", "FPN", "DMN")),
                                 stringsAsFactors = FALSE)
      } else {
        message("Compiling communities based on Louvain fast-unfolding at d_12 for N=83")
        community <- readRDS(paste0(getwd(), "/cache/power269/d12_louv_n83.rds"))
        
        membership <- data.frame(nodenum=names(community$membership), community = community$membership, 
                                 community.name=mapvalues(community$membership, from = c("1","2","3","4","5"), to = c("SOMMOTOR", "FPN/DA", "OCC", "CO", "DMN")),
                                 stringsAsFactors = FALSE)
      }
      
      
      all.siglm.nodal <- dplyr::left_join(results.lm, membership, by = "nodenum") %>% arrange(term, metric, nodenum)
      all.sigwlelm.nodal <- dplyr::left_join(results.wle.lm, membership, by = "nodenum") %>% arrange(term, metric, nodenum)
      all.sigttest.nodal <- dplyr::left_join(results.ttestbpd, membership, by = "nodenum") %>% arrange(metric, nodenum)
      
      write.csv(all.sigwlelm.nodal, file = paste0(basedir, "/results/", parcellation, "_", conn_method,  "/all.sigwlelm.nodal.weighted_", data.reduce, ".csv"), row.names=FALSE)
      write.csv(all.siglm.nodal, file = paste0(basedir, "/results/", parcellation, "_", conn_method, "/all.siglm.nodal.weighted_", data.reduce, ".csv"), row.names=FALSE)
      write.csv(all.sigttest.nodal, file = paste0(basedir, "/results/", parcellation, "_", conn_method,  "/all.sigttest.weighted_", data.reduce, ".csv"), row.names=FALSE)
    
    }
    
    return(list(all.sigwlelm.nodal=all.sigwlelm.nodal,
                all.siglm.nodal=all.siglm.nodal, 
                all.sigttest.nodal=all.sigttest.nodal)) 
   }

  }
