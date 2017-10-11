plot_significant_ageeffects <- function(results){
  
age.main <- subset(results, term == "age")
age.main.all <- data.frame()

pdf(file = paste0(basedir, "/Figures/age_main_plot_", parcellation, "_", preproc_pipeline, "_", conn_method, "_", PCA.method,".pdf"), width=10, height=7)
for(res in 1:length(age.main[,1])){
  ss.age <- data.frame(subj_info[,c("BPD", "AgeAtScan")], metric = age.main[res,]$metric)
  toanalyze$node <- as.character(toanalyze$node)
  age.main$nodenum <- as.character(age.main$nodenum) 
  ss.age$metric.val <- select_(toanalyze[which(toanalyze$node == age.main[res,]$nodenum),], as.character(ss.age$metric[1])) 
  ss.age$metric.val <- as.numeric(as.matrix(ss.age$metric.val))
  ss.age$BPD <- factor(ss.age$BPD, levels = c(0,1), labels = c("Control", "BPD"))
  ss.age$roiname <- paste0("V_", age.main[res,]$nodenum, ": ", age.main[res,]$nodename)
  ss.age$roinum <- age.main[res,]$nodenum
  for(m in 1:length(metrics.pca)){
    if(ss.age$metric[1] == metrics.pca[m]){
      ss.age$metric.label <- PClabs[m]
    }
  }
   colnames(ss.age) <- c("BPD", "AgeAtScan", "metric", "metric.val", "roiname", "roinum", "metric.label")
  g <- ggplot(ss.age, aes(x = AgeAtScan, y = metric.val)) + geom_point() + stat_smooth(method="lm", se=FALSE)+ 
    labs(x = "Age", y = as.character(ss.age$metric.label[1]), title = as.character(ss.age$roiname[1])) +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))
  plot(g)
  # browser()
  age.main.all <- rbind(age.main.all, ss.age)
}
message("exporting significant age effects to: ", paste0(basedir, "/Figures/age_main_plot_", parcellation, "_", preproc_pipeline, "_", conn_method, "_", PCA.method,".pdf"))
dev.off()
}