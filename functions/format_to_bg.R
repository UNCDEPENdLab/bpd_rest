format_to_bg <- function(allg_density_fc, covars, thresholds, allowCache = TRUE){
  
  cache_file <- paste0(basedir, "/cache/allg_density_bg_", file_tag, ".RData")
  
  
  if(allowCache == TRUE & file.exists(cache_file)){
    load(cache_file)
    message("loading brainGraph-ifized file from: ", cache_file)
  } else {
  
  #output: 3-deep list list[[group]][[threshold]][[subject]]
  allg_density_bg <- list()
  for(group in unique(covars$Group)){
    # browser()
    subdf <- covars[which(covars$Group == group),] %>% dplyr::select(Study.ID)
    subjs <- subdf$Study.ID
    group_by_thresh <- replicate(length(thresholds),list())
    for (tr in 1:length(thresholds)){
      for (sub in subjs){
        group_by_thresh[[tr]][[sub]] <- allg_density_fc[[sub]][[tr]]
      }
    }
    
    allg_density_fc_attrs <- lapply(group_by_thresh, function(thresh){
      lapply(thresh, function(g){
        ##cherry-picked from set_brainGraph_attr() to include the contents of calcGraph functions 
        # browser()
        ##avoid naming issues in glm 
        V(g)$name <- atlas$name
        g$name <- g$id
         
        ##global
        g$E.global <- efficiency(g, 'global', weights=NA)
        g$E.density <- edge_density(g) #does not use edge weights
        g$Lp <- mean_distance(g) #Characterisitc path length does not use edge weights
        g$CC <- transitivity(g, type = "global") #uses edge weights if the graph has an edge weight attribute
        g$modularity <- modularity(g, vertex_attr(g, "comm"))
        
        #nodal
        V(g)$degree <- degree(g)
        V(g)$PC <- part_coeff(g, V(g)$comm)
        V(g)$ev.cent <- centr_eigen(g)$vector
        V(g)$btwn.cent <- centr_betw(g)$res
        wibw <- wibw_module_degree(g, community_attr= "comm")
        V(g)$within.mod <- wibw$Ki
      
        
        #edge
        E(g)$btwn <- edge.betweenness(g)
        return(g)
      })
    })
    allg_density_bg[[group]] <- allg_density_fc_attrs
  }
  save(file = cache_file, allg_density_bg)
  }
  return(allg_density_bg)
}
