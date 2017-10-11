##export node files and BNV jpegs for significant effects
plot_sig_results_nodal <- function(signod.bpd, signod.lm, metrics.by.density, metrics.pca){ 
  outputdir <- paste0(basedir, "/BNV_nodefiles/OHBM_files/")
  stopifnot(file.exists(file.path(basedir, "BNV_nodefiles")))
  stopifnot(file.exists(file.path(basedir, "cache")))  
  
  ##main effect of BPD: 
  ########
  bpd.main <- subset(signod.bpd, select = c(nodename, metric, nodenum, community, community.name))
  row.names(bpd.main) <- seq(1, length(bpd.main[,1]), 1)
  bpd.main.all <- data.frame()
  
  pdf(file = paste0(outputdir, "bpd_main_plot_PCA.pdf"), width=10, height=7)
  for(res in 1:length(bpd.main[,1])){
    ss.bpd <- data.frame(subj_info[,c("SPECC_ID", "BPD", "AgeAtScan")], metric = bpd.main[res,]$metric)
    metrics.by.density$node <- as.character(metrics.by.density$node)
    metric.val <- subset(metrics.by.density, node == as.character(bpd.main[res,]$nodenum), select = c("id", as.character(ss.bpd$metric[1])))
    metric.val$id <- as.character(metric.val$id)
    colnames(metric.val) <- c("SPECC_ID", "metric.val")
    ss.bpd <- left_join(ss.bpd, metric.val, by = "SPECC_ID")
    ss.bpd <- ss.bpd %>% select(-SPECC_ID)
    #ss.bpd$metric.val <- select_(metrics.by.density[which(metrics.by.density$node == bpd.main[res,]$nodenum),], as.character(ss.bpd$metric[1])) 
    #ss.bpd$metric.val <- as.numeric(as.matrix(ss.bpd$metric.val))
    
    
    ss.bpd$BPD <- factor(ss.bpd$BPD, levels = c(0,1), labels = c("Control", "BPD"))
    ss.bpd$roiname <- paste0(bpd.main[res,]$nodenum, ": ", bpd.main[res,]$nodename)
    ss.bpd$roinum <- bpd.main[res,]$nodenum
    if(ss.bpd$metric[1] == "central"){
      ss.bpd$metric.label <- "PC1: Eigenvector and Degree Centrality"
    } else if(ss.bpd$metric[1] == "within.mod"){
      ss.bpd$metric.label <- "PC2: Within Module Connectivity"
    } else if (ss.bpd$metric[1] == "between.node") {
      ss.bpd$metric.label <- "PC3: Betweenness Centrality"
    } else if (ss.bpd$metric[1] == "integration") {
      ss.bpd$metric.label <- "PC3: Integration"
    } else if (ss.bpd$metric[1] == "closeness") {
      ss.bpd$metric.label <- "PC4: Closeness Centrality"
    } else { stop("what is: ", ss.bpd$metric[1]) }
    
    colnames(ss.bpd) <- c("BPD", "AgeAtScan", "metric", "metric.val", "roiname", "roinum", "metric.label")
    g <- ggplot(ss.bpd, aes(x = AgeAtScan, color=BPD, y = metric.val)) + geom_point() + stat_smooth(method="lm", se=FALSE)+ 
      labs(x = "Age", y = as.character(ss.bpd$metric.label[1]), title = as.character(ss.bpd$roiname[1])) +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))
    plot(g)
    
    bpd.main.all <- rbind(bpd.main.all, ss.bpd)
  }
  ##would like to have each of these before the line plots but, need to calculate bpd.main.all first 
 
  pdf(file = paste0(outputdir, "bpd_main_central_interesting.pdf"), width=10, height=7)
   for (m in metrics.pca){
    this.metric <- bpd.main.all[which(bpd.main.all$metric ==m), ]
    this.metric$roiname <- factor(this.metric$roiname, levels = unique(this.metric$roiname))
    g <- ggplot(this.metric, aes(x = factor(BPD), y = this.metric[,"metric.val"], color = BPD)) + stat_summary(fun.data="mean_cl_boot", size=1.5, fatten=1.5) + theme(legend.position="none") +   #fill = factor(BPD)
      labs(x = "", y = this.metric$metric.label[1], title = "") + scale_x_discrete(breaks=c(0,1), labels=c("Control", "BPD"))  + facet_wrap(~roiname, scales="free_y")  +
      theme_grey(base_size = 10) + theme(legend.title=element_blank(), strip.text = element_text(size=10)) + theme(legend.position="bottom")
    plot(g)
  }
  dev.off()

####Age Main effects
age.main <- subset(signod.wle, term == "age")
age.main <- age.main[order(age.main$metric, age.main$nodenum),]
row.names(age.main) <- seq(1, length(age.main[,1]), 1)
age.main.all <- data.frame()

pdf(file = paste0(outputdir, "age_main_plot_PCA.pdf"), width=10, height=7)
for(res in 1:length(age.main[,1])){
 
  #ss.bpd <- data.frame(subj_info[,c("SPECC_ID", "BPD", "AgeAtScan")], metric = bpd.main[res,]$metric)
  # 
  # metric.val <- subset(metrics.by.density, node == as.character(bpd.main[res,]$nodenum), select = c("id", as.character(ss.bpd$metric[1])))
  # metric.val$id <- as.character(metric.val$id)
  # colnames(metric.val) <- c("SPECC_ID", "metric.val")
  # ss.bpd <- left_join(ss.bpd, metric.val, by = "SPECC_ID")
  # ss.bpd <- ss.bpd %>% select(-SPECC_ID)
  # 
  ss.age <- data.frame(subj_info[,c("SPECC_ID", "BPD", "AgeAtScan")], metric = age.main[res,]$metric)
  metric.val <- subset(metrics.by.density, node == as.character(age.main[res,]$nodenum), select = c("id", as.character(ss.age$metric[1])))
  metric.val$id <- as.character(metric.val$id)
  colnames(metric.val) <- c("SPECC_ID", "metric.val")
  ss.age <- left_join(ss.age, metric.val, by = "SPECC_ID")
  ss.age <- ss.age %>% select(-SPECC_ID)
  # ss.age$metric.val <- select_(metrics.by.density[which(metrics.by.density$node == age.main[res,]$nodenum),], as.character(ss.age$metric[1])) 
  # ss.age$metric.val <- as.numeric(as.matrix(ss.age$metric.val))
  ss.age$BPD <- factor(ss.age$BPD, levels = c(0,1), labels = c("Control", "BPD"))
  ss.age$roiname <- paste0(age.main[res,]$nodenum, ": ", age.main[res,]$nodename)
  ss.age$roinum <- age.main[res,]$nodenum
  if(ss.age$metric[1] == "central"){
    ss.age$metric <- "PC1: Eigenvector and Degree Centrality"
  } else if(ss.age$metric[1] == "within.mod"){
    ss.age$metric <- "PC2: Within Module Connectivity"
  } else if (ss.age$metric[1] == "between.node") {
    ss.bpd$metric.label <- "PC3: Betweenness Centrality"
  } else if (ss.age$metric[1] == "integration") {
    ss.bpd$metric.label <- "PC3: Integration"
  } else if (ss.age$metric[1] == "closeness") {
    ss.bpd$metric.label <- "PC4: Closeness Centrality"
  } else { stop("what is: ", ss.bpd$metric[res]) }
  
 
  
  g <- ggplot(ss.age, aes(x = AgeAtScan, color=BPD, y = metric.val)) + geom_point() + stat_smooth(method="lm", se=FALSE)+ 
    labs(x = "Age", y = as.character(ss.age$metric[1]), title = as.character(ss.age$roiname[1])) +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))
  plot(g)
  # browser()
  age.main.all <- rbind(age.main.all, ss.age)
}
dev.off()

####Age x BPD interactions
age.bpd.ixn <- subset(signod.wle, term == "age:groupBPD")
age.bpd.ixn <- age.bpd.ixn[order(age.bpd.ixn$metric, age.bpd.ixn$nodenum),]
row.names(age.bpd.ixn) <- seq(1, length(age.bpd.ixn[,1]),1)
age.bpd.ixn.all <- data.frame()

pdf(file = paste0(outputdir, "agexbpd_ixns_plot_PCA.pdf"), width=10, height=7)
for(res in 1:length(age.bpd.ixn[,1])){
  # ss.ixn <- data.frame(subj_info[,c("BPD", "AgeAtScan")], metric = age.bpd.ixn[res,]$metric)
  # ss.ixn$metric.val <- select_(metrics.by.density[which(metrics.by.density$node == age.bpd.ixn[res,]$nodenum),], as.character(ss.ixn$metric[1])) 
  # ss.ixn$metric.val <- as.numeric(as.matrix(ss.ixn$metric.val))
  
  ss.ixn <- data.frame(subj_info[,c("SPECC_ID", "BPD", "AgeAtScan")], metric = age.bpd.ixn[res,]$metric)
  metric.val <- subset(metrics.by.density, node == as.character(age.bpd.ixn[res,]$nodenum), select = c("id", as.character(ss.ixn$metric[1])))
  metric.val$id <- as.character(metric.val$id)
  colnames(metric.val) <- c("SPECC_ID", "metric.val")
  ss.ixn <- left_join(ss.ixn, metric.val, by = "SPECC_ID")
  ss.age <- ss.ixn %>% select(-SPECC_ID)
  
  ss.ixn$BPD <- factor(ss.ixn$BPD, levels = c(0,1), labels = c("Control", "BPD"))
  ss.ixn$roiname <- paste0(age.bpd.ixn[res,]$nodenum, ": ", age.bpd.ixn[res,]$nodename)
  ss.ixn$roinum <- age.bpd.ixn[res,]$nodenum
  if(ss.ixn$metric[1] == "central"){
    ss.ixn$metric <- "PC1: Eigenvector and Degree Centrality"
  } else if(ss.age$metric[1] == "within.mod"){
    ss.age$metric <- "PC2: Within Module Connectivity"
  } else if (ss.age$metric[1] == "between.node") {
    ss.bpd$metric.label <- "PC3: Betweenness Centrality"
  } else if (ss.age$metric[1] == "integration") {
    ss.bpd$metric.label <- "PC3: Integration"
  } else if (ss.age$metric[1] == "closeness") {
    ss.bpd$metric.label <- "PC4: Closeness Centrality"
  } else { stop("what is: ", ss.bpd$metric[res]) }
  
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


# mean.g.infomap <- readRDS(paste0(basedir, "/cachedRfiles/infomap_communitylist.rds"))
# # basedir <- "~/Box Sync/DEPENd/Projects/RS_BPD_graph"
# # mean.g <- apply(allmats, c(2,3), aggfun, na.rm = TRUE)
# df <- as.data.frame(table(mean.g.infomap[[10]]$membership))
# badcomm <- df$Var1[df$Freq < 4]
# goodcomm <- as.numeric(df$Var1[df$Freq >= 4])
# mean.g.infomap[[10]]$membership[mean.g.infomap[[10]]$membership %in% badcomm] <- length(goodcomm) + 1   #where highest community number is extraneous nodes
#community <- mean.g.infomap[[10]]$membership



compare.list <- list()

metrics.bpd.main <- "integration"

###following computes a compare list coded to indicate the direction of significant effects for all nodes 
###compare.string: 0 = no diff, 1 = BPD > Control, 2 = BPD < Control
for (m in metrics.pca){
  this.metric <- bpd.main.all[which(bpd.main.all$metric == m), ]
  # this.metric$compare <- rep(NA, nnodes)
  compare.string <- rep(NA, nnodes +2) #we dropped 249 and 250 due to unstable TS data
  for(roi in unique(this.metric$roinum)){
    # roi <- as.numeric(gsub("V", "", roi))
    this.roi <- this.metric[which(this.metric$roinum == roi),]
    mean.Control <- mean(this.roi[which(this.roi$BPD == "Control"),"metric.val"], na.rm = TRUE)
    mean.BPD <- mean(this.roi[which(this.roi$BPD == "BPD"),"metric.val"], na.rm = TRUE)
    
    ##assumes proper ordering
    if(mean.Control > mean.BPD){compare.string[as.numeric(gsub("V", "", roi))] <- 1} else {compare.string[as.numeric(gsub("V", "", roi))] <- 2}
  }
  compare.string[which(is.na(compare.string))] <- 0
  compare.string <- compare.string[c(-249, -250)]
  compare.list[[m]] <- compare.string
}  

####export node files to be plotted 
# node.file.list <- list()
# for(m in metrics.pca){
#   this.metric <- bpd.main.all[which(bpd.main.all$metric == m),]
#   nodestp <- unique(this.metric$roinum)
#   Node.File <- NodeFile(atlas = atlas, 
#                         community = community, 
#                         nodestp = nodestp, 
#                         nodevals = compare.list[[m]], ##can make nodal statistics etc.
#                         nnodes = nnodes, 
#                         labels = 0, 
#                         filename = m, 
#                         outputdir = outputdir)
#   node.file.list[[m]] <- Node.File 
# }
# 
# }

# compare.string.allmetrics <- compare.list[[1]] + compare.list[[2]] + compare.list[[3]]
# names(compare.string.allmetrics) <- names(community$membership)
names(compare.string) <- names(community$membership)
bpd.main.all.integration <- subset(bpd.main.all, metric == "integration")
nodestp.bpd.integration <- unique(bpd.main.all.integration$roinum)
bpd.main.all.wit <- subset(bpd.main.all, metric == "within.mod")
nodestp.bpd.wit <- unique(bpd.main.all.wit$roinum)


membership <- community$membership
names(membership) <- NULL

nf.comm <- NodeFile(atlas = atlas,
                        community = membership,
                        nnodes = nnodes,
                        labels = 0,
                        filename = "OHBM_bpd_all_communities",
                        outputdir = outputdir)




nodestp.ixn <- unique(age.bpd.ixn.all$roinum)

nf.ixn <- NodeFile(atlas = atlas,
                   community = membership,
                   nodestp = nodestp.ixn,
                   labels = 0,
                   filename = "OHBM_ixn",
                   outputdir = outputdir
                    )
}
