reduce_centrality_global <- function(global.df, global.metrics, allowCache = TRUE){
  
  suppressMessages(require(wle))
  suppressMessages(require(psych))
  suppressMessages(require(reshape2))
  suppressMessages(require(tidyr))
  suppressMessages(require(tidyverse))
  
  fa.CFI<-function(x){
    nombre<-paste(x,"CFI",sep = ".")
    nombre<-
      ((x$null.chisq-x$null.dof)-(x$STATISTIC-x$dof))/(x$null.chisq-x$null.dof)
    return(nombre)
  }
  
  stopifnot(file.exists(file.path(basedir, "cache")))
  
  expectFile <- file.path(basedir, "cache", paste0("toanalyze.fa.thresh_global", file_tag, ".RData"))
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading global toanalyze from file: ", expectFile)
    if (thresh == "fc" && thresh_weighted =="binary"){
      output <- get(load(expectFile))
    } else {
      message("only FC + binary supported")
    }
  } else {
  
    browser()
    head(global.df.raw)
    
    global.df.raw <- global.df %>% gather(key="variable", value="value", -id, -density)
    
    cor(global.df[,1:4])
    
    #split for multiple CFAs
    raw.path.length <- global.df.raw %>% dplyr::filter(variable == "characteristic_path_length")
    raw.clustering <- global.df.raw %>% dplyr::filter(variable == "clustering_coefficient")
    raw.sw <- global.df.raw %>% dplyr::filter(variable == "small_worldness")
    raw.mod <- global.df.raw %>% dplyr::filter(variable == "modularity")
    
    
    
    ####path length
    
    path.length.fa <- dcast(raw.path.length, id ~ density + variable, value.var = "value")
    
    faout_1 <- fa(dplyr::select(path.length.fa, -id), nfactors = 1,  missing = TRUE, rotate = "varimax")
    print(faout_1$loadings, cutoff = .3)
    summary(faout_1)
    fa.CFI(faout_1)
    faout_1$scores
    
    # faout_2 <- fa(dplyr::select(path.length.fa, -id), nfactors = 2,  missing = TRUE, rotate = "varimax")
    # print(faout_2$loadings, cutoff = .3)
    # summary(faout_2)
    # fa.CFI(faout_2)
    # 
    
    ####clustering coefficient
    
    clustering.fa <- dcast(raw.clustering, id ~ density + variable, value.var = "value")
    
    faout_1 <- fa(dplyr::select(clustering.fa, -id), nfactors = 1,  missing = TRUE, rotate = "varimax")
    print(faout_1$loadings, cutoff = .3)
    summary(faout_1)
    fa.CFI(faout_1)
    faout_1$scores
    
    ###small world-ness
    ####path length
    
    sw.fa <- dcast(raw.sw, id ~ density + variable, value.var = "value")
    
    faout_1 <- fa(dplyr::select(sw.fa, -id), nfactors = 1,  missing = TRUE, rotate = "varimax")
    print(faout_1$loadings, cutoff = .3)
    summary(faout_1)
    fa.CFI(faout_1)
    faout_1$scores
    
    ###modularity
    ####path length
    
    mod.fa <- dcast(raw.mod, id ~ density + variable, value.var = "value")
    
    faout_1 <- fa(dplyr::select(mod.fa, -id), nfactors = 1,  missing = TRUE, rotate = "varimax")
    print(faout_1$loadings, cutoff = .3)
    summary(faout_1)
    fa.CFI(faout_1)
    faout_1$scores
    
  
    
    
    
    ##global: try 4 factor solution on all
    global.all.fa <- dcast(global.df.raw, id ~ density + variable, value.var = "value")
    
    faout_4 <- fa(dplyr::select(global.all.fa, -id), nfactors = 4,  missing = TRUE, rotate = "promax")
    print(faout_4$loadings)
    summary(faout_4)
    fa.CFI(faout_4)
    
    x <- pca(dplyr::select(global.all.fa, -id), nfactors = 4,  missing = TRUE, rotate = "promax") 
    print(x$loadings)
    
    
    }
}