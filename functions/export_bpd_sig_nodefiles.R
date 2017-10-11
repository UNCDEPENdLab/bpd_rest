export_bpd_sig_nodefiles <- function(bpd.main.all, outputdir){
compare.list <- list()

###following computes a compare list coded to indicate the direction of significant effects for all nodes 
###compare.string: 0 = no diff, 1 = BPD > Control, 2 = BPD < Control
for (m in metrics.pca){
  this.metric <- bpd.main.all[which(bpd.main.all$metric == m), ]
  # this.metric$compare <- rep(NA, nnodes)
  # compare.string <- rep(NA, nnodes +2) 
  compare.string <- rep(NA, nnodes) 
  for(roi in unique(this.metric$roinum)){
    this.roi <- this.metric[which(this.metric$roinum == roi),]
    mean.Control <- mean(this.roi[which(this.roi$BPD == "Control"),"metric.val"])
    mean.BPD <- mean(this.roi[which(this.roi$BPD == "BPD"),"metric.val"])
    
    if(is.na(mean.Control) | is.na(mean.BPD)){compare.string[roi] <- NA; next} #if we have an 
    ##assumes proper ordering
    if(mean.Control > mean.BPD){compare.string[roi] <- 1} else if (mean.Control < mean.BPD) {compare.string[roi] <- 2} else{compare.string[roi] <- 0}
  }
  compare.string[which(is.na(compare.string))] <- 0
  #compare.string <- compare.string[c(-249, -250)]
  compare.list[[m]] <- compare.string
}  


##Read in Yeo 7 networks
membership.file <- file.path(getwd(), "data", "membership.yeo.csv")
stopifnot(file.exists(membership.file))
community <- as.numeric(as.matrix(read.csv(membership.file)))
#names(community) <- paste0("V", seq(1,422,1))

#community <- make_clusters(agg.g, membership = community, algorithm = "Yeo_etal_2011_7Networks")

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
}


