plot_significant_ixn <- function(results){
  
  age.bpd.ixn <- subset(results, term == "age:groupBPD")
  age.bpd.ixn.all <- data.frame()
  
  pdf(file = paste0(basedir, "/Figures/age_bpd_ixn_plot_", parcellation, "_", preproc_pipeline, "_", conn_method, "_", PCA.method,".pdf"), width=10, height=7)
  for(res in 1:length(age.bpd.ixn[,1])){
   ss.age <- data.frame(subj_info[,c("BPD", "AgeAtScan")], metric = age.bpd.ixn[res,]$metric)
    toanalyze$node <- as.character(toanalyze$node)
    age.bpd.ixn$nodenum <- as.character(age.bpd.ixn$nodenum) 
    ss.age$metric.val <- select_(toanalyze[which(toanalyze$node == age.bpd.ixn[res,]$nodenum),], as.character(ss.age$metric[1])) 
    ss.age$metric.val <- as.numeric(as.matrix(ss.age$metric.val))
    ss.age$BPD <- factor(ss.age$BPD, levels = c(0,1), labels = c("Control", "BPD"))
    ss.age$roiname <- paste0(age.bpd.ixn[res,]$nodenum, ": ", age.bpd.ixn[res,]$nodename)
    ss.age$roinum <- age.bpd.ixn[res,]$nodenum
    for(m in 1:length(metrics.pca)){
      if(ss.age$metric[1] == metrics.pca[m]){
        ss.age$metric.label <- PClabs[m]
      }
    }
    colnames(ss.age) <- c("BPD", "AgeAtScan", "metric", "metric.val", "roiname", "roinum", "metric.label") 
    g <- ggplot(ss.age, aes(x = AgeAtScan, color = BPD, y = metric.val)) + geom_point() + stat_smooth(method="lm", se=FALSE)+ 
      labs(x = "Age", y = as.character(ss.age$metric.label[1]), title = as.character(ss.age$roiname[1])) +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))
    plot(g)
    # browser()
    age.bpd.ixn.all <- rbind(age.bpd.ixn.all, ss.age)
  }
  message("exporting significant age X BPD interactions to: ", paste0(basedir, "/Figures/age_bpd_ixn_plot_", parcellation, "_", preproc_pipeline, "_", conn_method, "_", PCA.method,".pdf"))
  dev.off()
}

ss.age <- data.frame(subj_info[,c("BPD", "AgeAtScan")], metric = "within.mod")
toanalyze$node <- as.character(toanalyze$node)
age.bpd.ixn$nodenum <- as.character(age.bpd.ixn$nodenum) 
ss.age$metric.val <- select_(toanalyze[which(toanalyze$node == age.bpd.ixn[res,]$nodenum),], as.character(ss.age$metric[1])) 
ss.age$metric.val <- as.numeric(as.matrix(ss.age$metric.val))
ss.age$BPD <- factor(ss.age$BPD, levels = c(0,1), labels = c("Control", "BPD"))
ss.age$roiname <- paste0(age.bpd.ixn[res,]$nodenum, ": ", age.bpd.ixn[res,]$nodename)
ss.age$roinum <- age.bpd.ixn[res,]$nodenum
for(m in 1:length(metrics.pca)){
  if(ss.age$metric[1] == metrics.pca[m]){
    ss.age$metric.label <- PClabs[m]
  }
}
colnames(ss.age) <- c("BPD", "AgeAtScan", "metric", "metric.val", "roiname", "roinum", "metric.label")

ss.rand <- toanalyze[which(toanalyze$node == "V392"),]
ss.rand$BPD <- factor(ss.rand$BPD, levels = c(0,1), labels = c("Control", "BPD"))
ss.rand <- ss.rand[2:8]


g <- ggplot(ss.rand, aes(x = Age, color = BPD, y = within.mod)) + geom_point() +stat_smooth(method = "lm")
plot(g)
