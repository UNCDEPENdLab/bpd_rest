plot_significant_groupeffects <- function(bpd.main){
  # bpd.main <- subset(all.sig.nodal.pca, is.na(term))
  
  # bpd.main <- signod.bpd
  bpd.main$ordernode <- as.numeric(substring(bpd.main$nodenum, 2))
  bpd.main <- bpd.main[order(bpd.main$metric, bpd.main$ordernode),]
  bpd.main <- select(bpd.main, -ordernode)
  row.names(bpd.main) <- seq(1, length(bpd.main$nodenum), 1)
  bpd.main.all <- data.frame()
  
  pdf(file = paste0(basedir, "/Figures/bpd_main_plot_", parcellation, "_", preproc_pipeline, "_", conn_method, "_", PCA.method, ".pdf"), width=10, height=7)
  for(res in 1:length(bpd.main[,1])){
    ss.bpd <- data.frame(subj_info[,c("BPD", "AgeAtScan")], metric = bpd.main[res,]$metric)
    toanalyze$node <- as.character(toanalyze$node)
    bpd.main$nodenum <- as.character(bpd.main$nodenum)
    
    ss.bpd$metric.val <- select_(toanalyze[which(toanalyze$node == bpd.main[res,]$nodenum),], as.character(ss.bpd$metric[1])) 
    
    ss.bpd$metric.val <- as.numeric(as.matrix(ss.bpd$metric.val))
    ss.bpd$BPD <- factor(ss.bpd$BPD, levels = c(0,1), labels = c("Control", "BPD"))
    ss.bpd$roiname <- paste0(bpd.main[res,]$nodenum, ": ", bpd.main[res,]$nodename)
    ss.bpd$roinum <- bpd.main[res,]$nodenum
    
    for(m in 1:length(metrics.pca)){
      if(ss.bpd$metric[1] == metrics.pca[m]){
        ss.bpd$metric.label <- PClabs[m]
      }
    }
    # 
    # if(ss.bpd$metric[1] == metrics.pca[1]){
    #   ss.bpd$metric.label <- PClabs[1]
    # } else if(ss.bpd$metric[1] == metrics.pca[2]){
    #   ss.bpd$metric.label <- "PC2: Betweenness Centrality"
    # } else if(ss.bpd$metric[1] == "within.mod"){
    #   ss.bpd$metric.label <- "PC3: Within Module Connectivity"
    # } else {
    #   ss.bpd$metric.label <- "PC4: Between Module Connectivity"
    # }
    # 
    colnames(ss.bpd) <- c("BPD", "AgeAtScan", "metric", "metric.val", "roiname", "roinum", "metric.label")
    # g <- ggplot(ss.bpd, aes(x = AgeAtScan, color=BPD, y = metric.val)) + geom_point() + stat_smooth(method="lm", se=FALSE)+ 
    #   labs(x = "Age", y = as.character(ss.bpd$metric.label[1]), title = as.character(ss.bpd$roiname[1])) +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))
    # plot(g)
    # 
    bpd.main.all <- rbind(bpd.main.all, ss.bpd)
  }
  
  for (m in metrics.pca){
    #bc some metrics may see no effects
    if(m %in% bpd.main.all$metric){
    this.metric <- bpd.main.all[which(bpd.main.all$metric ==m), ]
    this.metric$roiname <- factor(this.metric$roiname, levels = unique(this.metric$roiname))
    g <- ggplot(this.metric, aes(x = factor(BPD), y = this.metric[,"metric.val"], color = BPD)) + stat_summary(fun.data="mean_cl_boot", size=1.5, fatten=1.5) + theme(legend.position="none") +   #fill = factor(BPD)
      labs(x = "", y = this.metric$metric.label[1], title = "") + scale_x_discrete(breaks=c(0,1), labels=c("Control", "BPD"))  + facet_wrap(~roiname, scales="free_y")  +
      theme_grey(base_size = 10) + theme(legend.title=element_blank(), strip.text = element_text(size=6)) + theme(legend.position="bottom")
    plot(g)
    }
  }
  message("exporting significant group effects to: ", paste0(basedir, "/Figures/bpd_main_plot_", parcellation, "_", preproc_pipeline, "_", conn_method, "_", PCA.method,".pdf"))
  dev.off()
  return(bpd.main.all)
}