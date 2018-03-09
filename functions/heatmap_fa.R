heatmap_FA <- function(faout) {
  library(viridis)
  
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
  
  colnames(heat_plot) <- c("FC_Metric", "centrality", "integration",  "within.mod")
  
  xx <- melt(heat_plot, id = "FC_Metric", variable.name = "Factor", value.name = "Loading")
  
    plot <-  ggplot(xx, aes(Factor, FC_Metric, fill = Loading)) + geom_tile() +scale_fill_viridis(na.value = "transparent")#, option = "plasma")
  
  
  pdf(plot, file = "Factor_loading_heatplot.pdf", width = 8, height = 11)
  
  print(plot)
  
  dev.off()
  
  message("outputting to: Factor_loading_heatplot.pdf")
}
