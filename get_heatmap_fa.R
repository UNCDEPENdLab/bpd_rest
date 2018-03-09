heatmap_FA <- function(faout) {
  loadings <- data.frame(unclass(faout$loadings))
  
  vis_quants <- c(seq(1,20,4), 20)
  rs <- unique(gsub( "_.*$", "", row.names(loadings)))
  vis_rs <- rs[vis_quants]
  
  reducemetrics <- sort(reducemetrics)
  x <- list()
  
  for (r in 1:length(vis_rs)){
    rn <- data.frame(paste0(vis_rs[r], "_", reducemetrics))
    x[[r]] <- rn
    
  }
  
  heat <- data.frame(do.call(rbind, x))
  colnames(heat) <- "rn"
  
  loadings <- data.frame(setDT(loadings, keep.rownames = TRUE))
  
  heat_plot <- left_join(heat, loadings, by = "rn")
  
  colnames(heat_plot) <- c("FC_Metric", "central", "integration", "within.mod")
  
  xx <- melt(heat_plot, id = "FC_Metric", variable.name = "Factor", value.name = "Loading")
  
  ggplot(xx, aes(Factor, FC_Metric, fill = Loading)) + geom_tile() +scale_fill_viridis(na.value = "transparent")
  
  
  stat_density(aes(fill = ..density..), geom = "raster", position = "identity")
  
  ggplot(df,aes(x = Var1,y = Var2,fill = z)) + 
    geom_tile() + 
    scale_fill_gradient(low = "green", high = "red")
  (p <- ggplot(nba.m, aes(variable, Name)) + geom_tile(aes(fill = rescale),
                                                       +     colour = "white") + scale_fill_gradient(low = "white",
                                                                                                     +     high = "steelblue"))
  
  # rs <- round(as.numeric(unique(gsub( "_.*$", "", row.names(loadings)))), 4)
  # row.names(loadings)
  # 
  # 
  # 
  # 
  # test <- data.frame(rn =  melt(do.call(cbind, x))$value) 
  # 
  # # paste0(rep(reducemetrics, length(reducemetrics)), "_", rep(rs, length(rs)))
  # # 
  # # apply(expand.grid(rs, reducemetrics), 2, paste0, collapse = "_")
  # # 
  # # quant_rs <- as.numeric(quantile(rs))
  # # 
  # # rs_desired_log
  # # 
  # 
  # 
  # quant_rs <- c( .01, )
  # 
  # 
  # heat <- data.frame(rn = apply(expand.grid(quant_rs, reducemetrics), 1, paste0, collapse = "_"))
  # 
  # 
  # left_join(heat, loadings, by = "rn")
  # 
  # 
  # 
  
}
