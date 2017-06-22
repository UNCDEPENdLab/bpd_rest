#Analyses of nodal statistics using PCA to reduce across densities and metrics
#setup package dependencies and custom functions
#setwd("~/Box Sync/RS_BPD_graph")
# # setwd("/Users/mnh5174/Data_Analysis/bpd_rest")
# setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/")
# basedir <- getwd()

analyze_nodal_metrics_PCA <- function(nodalmetrics_dthresh_df = NULL, allowCache = TRUE){
  require(wle)  
  require(psych)
  
  stopifnot(file.exists(file.path(basedir, "cache")))  
  expectFile <- file.path(basedir, "output.files", paste0("all.siglm.nodal.pca.", preproc_pipeline, ".", conn_method, ".csv"))
  
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading significant results for PCA analysis from file: ", expectFile)
    #probably should write one .RData object for all three of these
    all.sigwlelm.nodal.pca <- read.csv(file = paste0(basedir, "/output.files/all.sigwlelm.nodal.pca.", preproc_pipeline, ".", conn_method, ".csv"))
    all.siglm.nodal.pca <- read.csv(file = paste0(basedir, "/output.files/all.siglm.nodal.pca.", preproc_pipeline, ".", conn_method, ".csv"))
    all.sigttest.nodal.pca <- read.csv(file = paste0(basedir, "/output.files/all.sigttest.nodal.pca.", preproc_pipeline, ".", conn_method, ".csv"))
    
    toanalyze <- read.csv(file = paste0(basedir, "/cache/toanalyze.pca.",preproc_pipeline, ".", conn_method, ".csv"))
  } else {
    #load nodal statistics for the current pipeline
    if(is.null(nodalmetrics_dthresh_df)){
      #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame
      nodalmetrics_dthresh_df <- load_nodal_metrics_df()
    } else {
      
      #not used at the moment
      #nodalmetrics_dthresh_df$density_fac <- factor(nodalmetrics_dthresh_df$density) #to overcome insanity of floating point imprecision if we use density==X approach 
      
      metrics_to_analyze <- c("degree", "eigen.cent", "betweenness.node", "within.module.deg.zscore", "between.module.deg.zscore")
      
      ######################################################################
      #PCA on latent degree/ eigenvector across densities
      
      #stack metrics in data.frame to have subjects, nodes, densities, *and* metrics on the rows
      #metrics.raw <- nodalmetrics_dthresh_df %>% select(id, node, density, degree, eigen.cent, betweenness.node, within.module.deg.zscore, between.module.deg.zscore) %>%
      #  filter(density > .04) %>% gather(key="variable", value="value", -id, -node, -density)
      
      #metrics including gateway and participation
      metrics.raw <- nodalmetrics_dthresh_df %>% select(id, node, density, degree, eigen.cent, betweenness.node, within.module.deg.zscore, between.module.deg.zscore,
                                                        gateway.coeff.btw, gateway.coeff.deg, part.coeff) %>%
        filter(density > .04) %>% gather(key="variable", value="value", -id, -node, -density)
      
      #na checks
      #natot <- function(x) { sum(is.na(x)) }
      #nodalmetrics_dthresh_df %>% filter(density > .04) %>% group_by(density) %>% summarize_at(vars(gateway.coeff.btw, gateway.coeff.deg, part.coeff), funs(natot))
      
      #recast for PCA such that column names represent a metrics at a given density (e.g., 0.05_between.module.deg.zscore)
      metrics.raw_pca <- dcast(metrics.raw, id + node ~ density + variable, value.var = "value")
      
      #does not handle NAs
      #metrics.pca <- prcomp(select(metrics.raw_pca, -id, -node), center = TRUE, scale. = TRUE)
      
      # pcaout1 <- pca(select(metrics.raw_pca, -id, -node), nfactors=1, rotate="varimax") #varimax (orthogonal) is default, just making that clear here
      # pcaout2 <- pca(select(metrics.raw_pca, -id, -node), nfactors=2, rotate="varimax")
      #this is the preferred solution for d > .04 with 8 nodal metrics (incl. gateway)
      pcaout3 <- pca(select(metrics.raw_pca, -id, -node), nfactors=3, rotate="varimax") #missing=TRUE, impute="median"
      #pcaout4 <- pca(select(metrics.raw_pca, -id, -node), nfactors=4, rotate="varimax")
      # pcaout5 <- pca(select(metrics.raw_pca, -id, -node), nfactors=5, rotate="varimax")
      # pcaout6 <- pca(select(metrics.raw_pca, -id, -node), nfactors=6, rotate="varimax")
      
      print(pcaout3$loadings, cutoff = 0.3) #if you need to look over PC loadings
      #print(pcaout4$loadings, cutoff = 0.3)
      
      pcasolution <- data.frame(pcaout3$scores) #use these, which are the rotated scores, not the raw eigenvectors from $x in prcomp
      #names(pcasolution) <- c("central", "between.node", "within.mod") #relevant for stats without gateway and participation
      names(pcasolution) <- c("central", "integration", "within.mod") #with participation
      
      #this was from Nate's analysis of Pearson. For the AROMA ridge, 3 factors covers 90% of variaence
      #pcasolution <- data.frame(metrics.pca$x[,1:4])
      #names(pcasolution) <- c("central", "between.node", "within.mod", "between.mod")
      
      toanalyze <- cbind(select(metrics.raw_pca, id, node), pcasolution)
      #subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz")
      merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)
      
      toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")
      write.csv(toanalyze, file = paste0(basedir, "/cache/toanalyze.pca.",preproc_pipeline, ".", conn_method, ".csv"))
      
      ##################################################################################################
      ##############
      ####run tests on nodal metrics collapsed across densities
      abst <- 2.64 ##about p < .005
      
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
          
          age.test.wle <- tryCatch(lm(metric~age*group, thismetric), error=function(e) { print(e); return(NULL) })
          if (!is.null(age.test.wle) && any(sigp <- summary(age.test.wle)$coefficients[-1,"Pr(>|t|)"] < .005)) {
            wle.lm_eff[[length(wle.lm_eff) + 1]] <- data.frame(summary(age.test.wle)$coefficients[c(FALSE, sigp),], nodenum=n, nodename=as.character(atlas$anat_label[atlas$name == n]), metric = m)
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

      community <- readRDS(paste0(getwd(), "/cache/d12_louv_n83.rds"))
      membership <- data.frame(nodenum=names(community$membership), community = community$membership, 
                               community.name=mapvalues(community$membership, from = c("1","2","3","4","5"), to = c("SOMMOTOR", "FPN/DA", "OCC", "CO", "DMN")),
                               stringsAsFactors = FALSE)

      all.siglm.nodal.pca <- dplyr::left_join(results.lm, membership, by = "nodenum")
      all.sigwlelm.nodal.pca <- dplyr::left_join(results.wle.lm, membership, by = "nodenum")
      all.sigttest.nodal.pca <- dplyr::left_join(results.ttestbpd, membership, by = "nodenum")
      
      write.csv(all.sigwlelm.nodal.pca, file = paste0(basedir, "/output.files/all.sigwlelm.nodal.pca.", preproc_pipeline, ".", conn_method, ".csv"))
      write.csv(all.siglm.nodal.pca, file = paste0(basedir, "/output.files/all.siglm.nodal.pca.", preproc_pipeline, ".", conn_method, ".csv"))
      write.csv(all.sigttest.nodal.pca, file = paste0(basedir, "/output.files/all.sigttest.nodal.pca.", preproc_pipeline, ".", conn_method, ".csv"))
    }
  }
  
  return(list(all.sigwlelm.nodal.pca=all.sigwlelm.nodal.pca,
              all.siglm.nodal.pca=all.siglm.nodal.pca, 
              all.sigttest.nodal.pca=all.sigttest.nodal.pca, 
              metrics.by.density=toanalyze))
}

