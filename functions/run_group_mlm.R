run_group_mlm <- function (toanalyze, allowCache = TRUE, browse = FALSE, weighted = FALSE){
  ##toanalyze must have be a df in long format with ID, node on the rows and a membership column to denote netwrk assignment
  #weighted is either TRUE or FALSE. Set to TRUE to run fully weighted analyses
  
  suppressMessages(require(tidyverse))
  suppressMessages(require(lme4))
  
  if (browse == TRUE){browser()} else{
    
    
    
    if(weighted == FALSE){
      
      stopifnot(file.exists(file.path(basedir, "cache")))  
      
      expectFile <- file.path(basedir, "cache", paste0("group_mlm_results_crossed_MH_update_", file_tag, ".RData"))
      
      if (file.exists(expectFile) && allowCache==TRUE) {
        message("Loading mlm results from file: ", expectFile)
        load(expectFile)
      } else {
        
        
        toanalyze_c <- toanalyze %>% mutate(membership=factor(membership), BPD = factor(BPD))#, levels=1:7, labels=c("VIsual")))
        metrics <- colnames(toanalyze_c)
        metrics <- toanalyze_c %>% dplyr::select(-id, -node, -BPD, -Age, -membership)
        #toanalyze
        
        
        metrics_results <- foreach(m=iter(colnames(metrics))) %dopar% {
          f <- as.formula(paste(m, "~ BPD*Age*membership + (1|id) + (1|node)"))
          group_mlm <- lmer(f, data = toanalyze_c)
          group_mlm #implicit return
        }
        
        metrics_results_crossed <- foreach(m=iter(colnames(metrics))) %dopar% {
          #test crossed effects: estimate subject-specific intercept that permits network-level variability
          f_crossed <- as.formula(paste(m, "~ BPD*Age*membership + (1 + membership|id) + (1|node)"))
          group_mlm_crossed <- lmer(f_crossed, data = toanalyze_c)
          group_mlm_crossed
        }
        
        results_for_comparison <- list(metrics_results, metrics_results_crossed)
        
        save(file = expectFile, results_for_comparison)
      }
      
    } else{
      
      stopifnot(file.exists(file.path(basedir, "cache")))  
      expectFile <- file.path(basedir, "cache", paste0("group_mlm_results_crossed_MH_update_", file_tag_nothresh, ".RData"))
      
      if (file.exists(expectFile) && allowCache==TRUE) {
        message("Loading mlm results from file: ", expectFile)
        load(expectFile)
      } else {
        
        
        toanalyze_c <- toanalyze %>% mutate(membership=factor(membership), BPD = factor(BPD))#, levels=1:7, labels=c("VIsual")))
        metrics <- colnames(toanalyze_c)
        metrics <- toanalyze_c %>% dplyr::select(-id, -node, -BPD, -Age, -membership)
        
        
        # metrics_results <- list()
        metrics_results <- foreach(m=iter(colnames(metrics))) %dopar% {
          f <- as.formula(paste(m, "~ BPD*Age*membership + (1|id) + (1|node)"))
          group_mlm <- lmer(f, data = toanalyze_c)
          group_mlm #implicit return
          
          #test crossed effects: estimate subject-specific intercept that permits network-level variability
          f_crossed <- as.formula(paste(m, "~ BPD*Age*membership + (1|id) + (membership|node)"))
          group_mlm_crossed <- lmer(f_crossed, data = toanalyze_c)
          group_mlm_crossed
        }
        # for(m in colnames(metrics)){
        #   f <- as.formula(paste(m, "~ BPD*Age*membership + (1|id) + (1|node)"))
        #   group_mlm <- lmer(f, data = toanalyze_c)
        #   metrics_results[[m]] <- group_mlm
        # }
        save(file = expectFile, metrics_results)
      }
      
    }
    
    if(weighted==FALSE){return(results_for_comparison)} else{return(metrics_results)}
  }
}

