run_group_comparisons_nodal <- function(toanalyze, abst = 2.64, allowCache = TRUE){
  stopifnot(file.exists(file.path(basedir, "cache")))  
  expectFile.ttest <- file.path(basedir, "output.files", paste0("all.sigttest.nodal.pca.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method, ".csv"))
  expectFile.wle <- file.path(basedir, "output.files", paste0("all.sigwlelm.nodal.pca.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method, ".csv"))
  expectFile.lm <- file.path(basedir, "output.files", paste0("all.siglm.nodal.pca.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method, ".csv"))
  
 
  if (file.exists(expectFile.ttest) && file.exists(expectFile.lm) && file.exists(expectFile.wle) && allowCache==TRUE) {
    message("Loading nodal metric results from: ", expectFile.ttest, ", ", expectFile.lm, ", ", expectFile.wle)
    all.sigwlelm.nodal.pca <- read.csv(expectFile.wle)
    all.siglm.nodal.pca <- read.csv(expectFile.lm) 
    all.sigttest.nodal.pca <- read.csv(expectFile.ttest)
  } else {
  
  metrics.toanalyze <- names(select(toanalyze, -id, -node, -BPD, -Age))
  ttest_bpd_eff <- list()
  lm_eff <- list() #all effects from lm
  wle.lm_eff <- list() #all effects from weighted likelihood (robust estimation)
  
  for (m in metrics.toanalyze) {
    for (n in levels(toanalyze$node)) {
      thismetric <- toanalyze[toanalyze$node == n, c(m, "BPD", "Age")]
      thismetric$BPD <- factor(thismetric$BPD, levels = c(0,1), labels = c("control", "BPD"))
      colnames(thismetric) <- c("metric", "group", "age")
      
      node.test <- tryCatch(t.test(metric~group, thismetric, na.rm = TRUE), error = function(errorname) { print(errorname); return(NULL) })
      if (is.null(node.test)) { message("Error occurred for metric: ", m, " and node: ", n) }
      
      if ((!is.null(node.test) && !is.nan(node.test$statistic) && abs(node.test$statistic) > abst)) {
        ttest_bpd_eff[[length(ttest_bpd_eff)+1]] <- data.frame(broom::tidy(node.test), nodename=as.character(atlas$anat_label[atlas$name == n]), metric=m, nodenum=n)
      }
      
      
      age.test.wle <- tryCatch(wle.lm(metric~age*group, thismetric), error=function(e) { print(e); return(NULL) })
      if (!is.null(age.test.wle) && any(sigp <- summary(age.test.wle)[[1]]$coefficients[-1,"Pr(>|t|)"] < .005)) { #summary returns a top-level list. Use [[1]] to access coefficients
        df <- data.frame(summary(age.test.wle)[[1]]$coefficients[c(FALSE, sigp),,drop=FALSE]) #use sigp to just pull sig rows. FALSE specifies to skip intercept since this was omitted in check above. drop=FALSE forces it to stay a data.frame
        df$term <- row.names(df); row.names(df) <- NULL #manual tidying... a bit prone to explosion if columns change
        names(df) <- c("estimate", "std.error", "statistic", "p.value", "term")
        df <- df[,c("term", "estimate", "std.error", "statistic", "p.value")] #rearrange for consistency
        
        wle.lm_eff[[length(wle.lm_eff) + 1]] <- data.frame(df, nodenum=n, nodename=as.character(atlas$anat_label[atlas$name == n]), metric = m)
      }
      
      age.test <- tryCatch(lm(metric~age*group, thismetric, na.action = "na.exclude"), error=function(e) { print(e); return(NULL) })
      if (is.null(age.test)) { pvec <- NULL
      } else { 
        agedf <- broom::tidy(age.test)
        pvec <- agedf$p.value[-1] #p-values of age, bpd, and age x bpd. -1 to drop off intercept
        names(pvec) <- agedf$term[-1]
      }
      
      if (!is.null(age.test) && !all(is.nan(pvec)) && pvec["groupBPD"] < .005) {
        lm_eff[[length(lm_eff) + 1]] <- data.frame(subset(agedf, term=="groupBPD"), nodenum=n, nodename=as.character(atlas$anat_label[atlas$name == n]), metric = m)
      }
      
      if (!is.null(age.test) && !all(is.nan(pvec)) && pvec["age"] < .005) {
        lm_eff[[length(lm_eff) + 1]] <- data.frame(subset(agedf, term=="age"), nodenum=n, nodename=as.character(atlas$anat_label[atlas$name == n]), metric = m)
      }
      
      if (!is.null(age.test) && !all(is.nan(pvec)) && pvec["age:groupBPD"] < .005) {
        lm_eff[[length(lm_eff) + 1]] <- data.frame(subset(agedf, term=="age:groupBPD"), nodenum=n, nodename=as.character(atlas$anat_label[atlas$name == n]), metric = m)
      }
    }
  }
  
  #results.df <- do.call(rbind, results)
  
  ####BPD main effects
  results.ttestbpd <- dplyr::bind_rows(ttest_bpd_eff) #efficient version of do.call(rbind, dflist) for list of data.frames
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
  
  
  all.siglm.nodal.pca <- dplyr::left_join(results.lm, membership, by = "nodenum") %>% arrange(term, metric, nodenum)
  all.sigwlelm.nodal.pca <- dplyr::left_join(results.wle.lm, membership, by = "nodenum") %>% arrange(term, metric, nodenum)
  all.sigttest.nodal.pca <- dplyr::left_join(results.ttestbpd, membership, by = "nodenum") %>% arrange(metric, nodenum)
  
  write.csv(all.sigwlelm.nodal.pca, file = paste0(basedir, "/output.files/all.sigwlelm.nodal.pca.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method, ".csv"), row.names=FALSE)
  write.csv(all.siglm.nodal.pca, file = paste0(basedir, "/output.files/all.siglm.nodal.pca.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method, ".csv"), row.names=FALSE)
  write.csv(all.sigttest.nodal.pca, file = paste0(basedir, "/output.files/all.sigttest.nodal.pca.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", PCA.method,".csv"), row.names=FALSE)
  
  }
  
  return(list(all.sigwlelm.nodal.pca=all.sigwlelm.nodal.pca,
              all.siglm.nodal.pca=all.siglm.nodal.pca, 
              all.sigttest.nodal.pca=all.sigttest.nodal.pca))
}
