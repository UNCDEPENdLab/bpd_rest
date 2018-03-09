plot_significant_groupeffects <- function(bpd.main, reducemetrics){
  # bpd.main <- subset(all.sig.nodal.pca, is.na(term))
  # browser()
  # bpd.main <- signod.bpd
  bpd.main$ordernode <- as.numeric(substring(bpd.main$nodenum, 2))
  bpd.main <- bpd.main[order(bpd.main$metric, bpd.main$ordernode),]
  bpd.main <- dplyr::select(bpd.main, -ordernode)
  row.names(bpd.main) <- seq(1, length(bpd.main$nodenum), 1)
  #set colors of rectangles to correspond to Yeo network colors
  bpd.main$shade <- ifelse(bpd.main$community.name == "SOMMOT", "blue",
                          ifelse(bpd.main$community.name == "FPN", "orange",
                                 ifelse(bpd.main$community.name == "DMN", "red",
                                        ifelse(bpd.main$community.name == "DORSATTN", "green",
                                               ifelse(bpd.main$community.name == "SALVENTATTN", "violet",
                                                      ifelse(bpd.main$community.name == "VIS", "purple", "papayawhip"))))))
  
  bpd.main.all <- data.frame()
  
  
  
  for(res in 1:length(bpd.main[,1])){
    ss.bpd <- data.frame(subj_info[,c("BPD", "AgeAtScan")], metric = bpd.main[res,]$metric)
    toanalyze$node <- as.character(toanalyze$node)
    bpd.main$nodenum <- as.character(bpd.main$nodenum)
    
    ss.bpd$metric.val <- select_(toanalyze[which(toanalyze$node == bpd.main[res,]$nodenum),], as.character(ss.bpd$metric[1])) 
    
    ss.bpd$metric.val <- as.numeric(as.matrix(ss.bpd$metric.val))
    ss.bpd$BPD <- factor(ss.bpd$BPD, levels = c(0,1), labels = c("Control", "BPD"))
    ss.bpd$roiname <- paste0(bpd.main[res,]$nodenum, ": ", bpd.main[res,]$nodename)
    ss.bpd$roinum <- bpd.main[res,]$nodenum
    
    ss.bpd$metric <- levels(droplevels(ss.bpd$metric))
    ss.bpd$network <- bpd.main[res,]$community.name
    ss.bpd$shade <- bpd.main[res,]$shade
    
    
    for(m in 1:length(fa.metrics)){
      if(ss.bpd$metric[1] == fa.metrics[m]){
        ss.bpd$metric.label <- fa.metrics[m]
      }
    }

    # if(ss.bpd$metric[1] == reducemetrics[1]){
    #   ss.bpd$metric.label <- PClabs[1]
    # } else if(ss.bpd$metric[1] == reducemetrics[2]){
    #   ss.bpd$metric.label <- "PC2: Betweenness Centrality"
    # } else if(ss.bpd$metric[1] == "within.mod"){
    #   ss.bpd$metric.label <- "PC3: Within Module Connectivity"
    # } else {
    #   ss.bpd$metric.label <- "PC4: Between Module Connectivity"
    # }
    # 
    #colnames(ss.bpd) <- c("BPD", "AgeAtScan", "metric", "metric.val", "roiname", "roinum", "metric.label")
    # g <- ggplot(ss.bpd, aes(x = AgeAtScan, color=BPD, y = metric.val)) + geom_point() + stat_smooth(method="lm", se=FALSE)+ 
    #   labs(x = "Age", y = as.character(ss.bpd$metric.label[1]), title = as.character(ss.bpd$roiname[1])) +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))
    # plot(g)
    # 
    bpd.main.all <- rbind(bpd.main.all, ss.bpd)
  }
  
  # m <- "central"
  
  pdf(file = paste0(basedir, "/Figures/bpd_main_plot_", parcellation, "_", preproc_pipeline, "_", conn_method, "_", data.reduce, ".pdf"), width=10, height=7)
  for (m in fa.metrics){
    #bc some metrics may see no effects
    if(m %in% bpd.main.all$metric){
    this.metric <- bpd.main.all[which(bpd.main.all$metric ==m), ]
    this.metric$roiname <- factor(this.metric$roiname, levels = unique(this.metric$roiname))
    
    g <- ggplot(this.metric, aes(x = factor(BPD), y = this.metric[,"metric.val"], color = BPD)) + stat_summary(fun.data="mean_cl_boot", size=1.5, fatten=1.5) + theme(legend.position="none") +   #fill = factor(BPD)
      labs(x = "", y = this.metric$metric.label[1], title = "") + scale_x_discrete(breaks=c(0,1), labels=c("Control", "BPD"))  + facet_wrap(~roiname, scales="free_y")  +
      theme_bw() + theme(legend.title=element_blank(), strip.text = element_text(size=6)) + theme(legend.position="bottom") + scale_color_brewer("BPD", palette="Set1") +
      theme(strip.background=element_blank())
    plot(g)
   

# trying to get rects to sort by network ----------------------------------

    
     # dummy <- g
    # dummy$layers <- NULL
    # dummy <- dummy +geom_rect(data = this.metric,xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, aes(fill=shade))
    # 
    # gtable_select <- function (x, ...) 
    # {
    #   matches <- c(...)
    #   x$layout <- x$layout[matches, , drop = FALSE]
    #   x$grobs <- x$grobs[matches]
    #   x
    # }
    # 
    # require(gtable); require(grid)
    # 
    # g1 <- ggplotGrob(g)
    # g2 <- ggplotGrob(dummy)
    #  
    # panels <- grepl(pattern="panel", g2$layout$name)
    # strips <- grepl(pattern="strip-right", g2$layout$name)
    # g2$grobs[strips] <- replicate(sum(strips), nullGrob(), simplify = FALSE)
    # g2$layout$l[panels] <- g2$layout$l[panels] + 1
    # g2$layout$r[panels] <- g2$layout$r[panels] + 2
    # 
    # new_strips <- gtable_select(g2, panels | strips)
    # grid.newpage()
    # grid.draw(new_strips)
    
    # dummy <- ggplot(data = this.metric, aes(y = this.metric[,"metric.val"]))+ facet_wrap(~roiname, scales="free_y") + 
    #   geom_rect( aes(fill = network),xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha = .6) +
    #   theme_minimal() #+ scale_fill_manual(values = shade)
    # 
     # plot(dummy)
    
   
    
    }
  }
  message("exporting significant group effects to: ", paste0(basedir, "/Figures/bpd_main_plot_", parcellation, "_", preproc_pipeline, "_", conn_method, "_", data.reduce,".pdf"))
  dev.off()
  return(bpd.main.all)
}