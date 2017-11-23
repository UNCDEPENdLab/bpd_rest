run_group_mlm <- function (toanalyze, allowCache = TRUE, browse = FALSE){
  ##toanalyze must have be a df in long format with ID, node on the rows and a membership column to denote netwrk assignment
  
  require(tidyverse)
  require(lme4)
  
  if (browse == TRUE){browser()} else{
  
  stopifnot(file.exists(file.path(basedir, "cache")))  
  expectFile <- file.path(basedir, "cache", paste0("mlm_results.", parcellation, ".", preproc_pipeline, ".", conn_method, ".", data.reduce, ".RData"))
  
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading mlm results from file: ", expectFile)
    load(expectFile)
  } else {
    
 
  toanalyze_c <- toanalyze %>% mutate(membership=factor(membership), BPD = factor(BPD))#, levels=1:7, labels=c("VIsual")))
  metrics <- colnames(toanalyze_c)
  metrics <- toanalyze_c %>% dplyr::select(-id, -node, -BPD, -Age, -membership)
  
  
  metrics_results <- list()
  for(m in colnames(metrics)){
    f <- as.formula(paste(m, "~ BPD*Age*membership + (1|id) + (1|node)"))
    group_mlm <- lmer(f, data = toanalyze_c)
    metrics_results[[m]] <- group_mlm
      }
  save(file = expectFile, metrics_results)
    }
  return(metrics_results)
  }
}

