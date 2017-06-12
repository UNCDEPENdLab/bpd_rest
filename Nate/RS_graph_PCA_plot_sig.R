#RS_graph_PCA_plot_sig findings
######plots for sig effects
##main effect of BPD: 
########
setwd("~/Box Sync/RS_BPD_graph")
basedir <- getwd()
outputdir <- paste0(getwd(), "/Viz_files")
source("Graph_util_redux.R")

atlas <- read.csv("power269_masterlookup_shift_nate.csv", header = TRUE)
all.sig.nodal.pca <- read.csv(paste0(getwd(), "/output.files/all.sig.nodal.pca.pearson.csv"))
toanalyze <- read.csv(paste0(getwd(), "/output.files/nodal.pca.toanalyze.pearson.csv"))
SPECC_rest <- read.csv("SPECC_rest.csv")
nnodes <- 269
metrics.pca <- c("central", "between.node", "within.mod", "between.mod")
library(ggplot2)
library(dplyr)

bpd.main <- subset(all.sig.nodal.pca, is.na(term))
bpd.main <- bpd.main[order(bpd.main$metric, bpd.main$nodenum),]
row.names(bpd.main) <- seq(1, length(bpd.main$X), 1)
bpd.main.all <- data.frame()

pdf(file = paste0(basedir, "/Figures/bpd_main_plot_PCA.pdf"), width=10, height=7)
for(res in 1:length(bpd.main[,1])){
  ss.bpd <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric = bpd.main[res,]$metric)
  ss.bpd$metric.val <- select_(toanalyze[which(toanalyze$node == bpd.main[res,]$nodenum),], as.character(ss.bpd$metric[1])) 
  ss.bpd$metric.val <- as.numeric(as.matrix(ss.bpd$metric.val))
  ss.bpd$BPD <- factor(ss.bpd$BPD, levels = c(0,1), labels = c("Control", "BPD"))
  ss.bpd$roiname <- paste0("V_", bpd.main[res,]$nodenum, ": ", bpd.main[res,]$nodename)
  ss.bpd$roinum <- bpd.main[res,]$nodenum
  if(ss.bpd$metric[1] == "central"){
    ss.bpd$metric.label <- "PC1: Eigenvector and Degree Centrality"
  } else if(ss.bpd$metric[1] == "between.node"){
    ss.bpd$metric.label <- "PC2: Betweenness Centrality"
  } else if(ss.bpd$metric[1] == "within.mod"){
    ss.bpd$metric.label <- "PC3: Within Module Connectivity"
  } else {
    ss.bpd$metric.label <- "PC4: Between Module Connectivity"
  }
  
  colnames(ss.bpd) <- c("BPD", "AgeAtScan", "metric", "metric.val", "roiname", "roinum", "metric.label")
  g <- ggplot(ss.bpd, aes(x = AgeAtScan, color=BPD, y = metric.val)) + geom_point() + stat_smooth(method="lm", se=FALSE)+ 
    labs(x = "Age", y = as.character(ss.bpd$metric.label[1]), title = as.character(ss.bpd$roiname[1])) +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))
  plot(g)
  
  bpd.main.all <- rbind(bpd.main.all, ss.bpd)
}

##would like to have each of these before the line plots but, need to calculate bpd.main.all first 
for (m in metrics.pca){
  this.metric <- bpd.main.all[which(bpd.main.all$metric ==m), ]
  this.metric$roiname <- factor(this.metric$roiname, levels = unique(this.metric$roiname))
  g <- ggplot(this.metric, aes(x = factor(BPD), y = this.metric[,"metric.val"], color = BPD)) + stat_summary(fun.data="mean_cl_boot", size=1.5, fatten=1.5) + theme(legend.position="none") +   #fill = factor(BPD)
    labs(x = "", y = this.metric$metric.label[1], title = "") + scale_x_discrete(breaks=c(0,1), labels=c("Control", "BPD"))  + facet_wrap(~roiname, scales="free_y")  +
    theme_grey(base_size = 12) + theme(legend.title=element_blank(), strip.text = element_text(size=10)) + theme(legend.position="bottom")
  plot(g)
}
dev.off()

####Age Main effects
age.main <- subset(all.sig.nodal.pca, term == "age")
age.main.all <- data.frame()

pdf(file = paste0(basedir, "/Figures/age_main_plot_PCA.pdf"), width=10, height=7)
for(res in 1:length(age.main[,1])){
  ss.age <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric = age.main[res,]$metric)
  ss.age$metric.val <- select_(toanalyze[which(toanalyze$node == age.main[res,]$nodenum),], as.character(ss.age$metric[1])) 
  ss.age$metric.val <- as.numeric(as.matrix(ss.age$metric.val))
  ss.age$BPD <- factor(ss.age$BPD, levels = c(0,1), labels = c("Control", "BPD"))
  ss.age$roiname <- paste0("V_", age.main[res,]$nodenum, ": ", age.main[res,]$nodename)
  ss.age$roinum <- age.main[res,]$nodenum
  if(ss.age$metric[1] == "central"){
    ss.age$metric <- "PC1: Eigenvector and Degree Centrality"
  } else if(ss.age$metric[1] == "between.node"){
    ss.age$metric <- "PC2: Betweenness Centrality"} else if(ss.age$metric[1] == "within.mod"){
      ss.age$metric <- "PC3: Within Module Connectivity"} else {ss.age$metric <- "PC4: Between Module Connectivity"}
  g <- ggplot(ss.age, aes(x = AgeAtScan, color=BPD, y = metric.val)) + geom_point() + stat_smooth(method="lm", se=FALSE)+ 
    labs(x = "Age", y = as.character(ss.age$metric[1]), title = as.character(ss.age$roiname[1])) +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))
  plot(g)
  # browser()
  age.main.all <- rbind(age.main.all, ss.age)
}
dev.off()

####Age x BPD interactions
age.bpd.ixn <- subset(all.sig.nodal.pca, term == "age:groupBPD")
age.bpd.ixn.all <- data.frame()

pdf(file = paste0(basedir, "/Figures/agexbpd_ixns_plot_PCA.pdf"), width=10, height=7)
for(res in 1:length(age.bpd.ixn[,1])){
  ss.ixn <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric = age.bpd.ixn[res,]$metric)
  ss.ixn$metric.val <- select_(toanalyze[which(toanalyze$node == age.bpd.ixn[res,]$nodenum),], as.character(ss.ixn$metric[1])) 
  ss.ixn$metric.val <- as.numeric(as.matrix(ss.ixn$metric.val))
  ss.ixn$BPD <- factor(ss.ixn$BPD, levels = c(0,1), labels = c("Control", "BPD"))
  ss.ixn$roiname <- paste0("V_", age.bpd.ixn[res,]$nodenum, ": ", age.bpd.ixn[res,]$nodename)
  ss.ixn$roinum <- age.bpd.ixn[res,]$nodenum
  if(ss.ixn$metric[1] == "central"){
    ss.ixn$metric <- "PC1: Eigenvector and Degree Centrality"
  } else if(ss.ixn$metric[1] == "between"){
    ss.ixn$metric <- "PC2: Betweenness Centrality"} else if(ss.ixn$metric[1] == "within.mod"){
    ss.ixn$metric <- "PC3: Within Module Connectivity"} else {ss.ixn$metric <- "PC4: Between Module Connectivity"}
  g <- ggplot(ss.ixn, aes(x = AgeAtScan, color=BPD, y = metric.val)) + geom_point() + stat_smooth(method="lm", se=FALSE)+ 
    labs(x = "Age", y = as.character(ss.ixn$metric[1]), title = as.character(ss.ixn$roiname[1])) +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))
  plot(g)
  g.loess <- ggplot(ss.ixn, aes(x = AgeAtScan, color=BPD, y = metric.val)) + geom_point() + stat_smooth()+#method="lm", se=FALSE)+ 
    labs(x = "Age", y = as.character(ss.ixn$metric[1]), title = as.character(ss.ixn$roiname[1])) +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))
  plot(g.loess)
  # browser()
  age.bpd.ixn.all <- rbind(age.bpd.ixn.all, ss.ixn)
}
dev.off()

###export to brainNet Viewer
####to be exported: BPD main effects (central, betweenness, between.mod, within.mod), age main effects, ixns
###export BPD main effects central
mean.g.infomap <- readRDS(paste0(basedir, "/cachedRfiles/infomap_communitylist.rds"))
df <- as.data.frame(table(mean.g.infomap[[10]]$membership))
badcomm <- df$Var1[df$Freq < 4]
goodcomm <- as.numeric(df$Var1[df$Freq >= 4])
mean.g.infomap[[10]]$membership[mean.g.infomap[[10]]$membership %in% badcomm] <- length(goodcomm) + 1   #where highest community number is extraneous nodes
community <- mean.g.infomap[[10]]$membership


compare.list <- list()

###following computes a compare list coded to indicate the direction of significant effects for all nodes 
###compare.string: 0 = no diff, 1 = BPD > Control, 2 = BPD < Control
for (m in metrics.pca){
  this.metric <- bpd.main.all[which(bpd.main.all$metric == m), ]
  # this.metric$compare <- rep(NA, nnodes)
  compare.string <- rep(NA, nnodes +2) #we dropped 249 and 250 due to unstable TS data
  for(roi in unique(this.metric$roinum)){
    this.roi <- this.metric[which(this.metric$roinum == roi),]
    mean.Control <- mean(this.roi[which(this.roi$BPD == "Control"),"metric.val"])
    mean.BPD <- mean(this.roi[which(this.roi$BPD == "BPD"),"metric.val"])
    
    ##assumes proper ordering
    if(mean.Control > mean.BPD){compare.string[roi] <- 1} else {compare.string[roi] <- 2}
  }
  compare.string[which(is.na(compare.string))] <- 0
  compare.string <- compare.string[c(-249, -250)]
  compare.list[[m]] <- compare.string
}  

####export node files to be plotted 
node.file.list <- list()
for(m in metrics.pca){
  this.metric <- bpd.main.all[which(bpd.main.all$metric == m),]
  nodestp <- unique(this.metric$roinum)
  Node.File <- NodeFile(atlas = atlas, 
                       community = community, 
                       nodestp = nodestp, 
                       nodevals = compare.list[[m]], ##can make nodal statistics etc.
                       nnodes = nnodes, 
                       labels = 0, 
                       filename = m, 
                       outputdir = outputdir)
  node.file.list[[m]] <- Node.File 
}

