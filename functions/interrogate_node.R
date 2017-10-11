interrogate_node <- function(allmats, subj_info, node, quant=NULL, t.stat=NULL, outputdir){
  if(!is.null(quant) & !is.null(t.stat)){
    message("Choose either quantile or t statistic to view")
  }

  edge.vals.vec <- allmats[,node,]
  
  controls <- which(subj_info$BPD == 0)
  bpd <- which(subj_info$BPD == 1)
 
  controlcorrs <- edge.vals.vec[controls, ]
  bpdcorrs <- edge.vals.vec[bpd,]
  
  edgepvals <- matrix(NA, nnodes)
  edgetstats <- matrix(NA, nnodes)
  for(n in 1:nnodes){
  res <- t.test(atanh(controlcorrs[,n]), atanh(bpdcorrs[,n]))
  # edgepvals[n,] <- res$p.value
  edgetstats[n,] <- res$statistic 
  
  colnames(edgetstats) <- "t.stat"
  ##returns a few NaNs, substitute p-value of 1 and t statisitic of 0
  if(is.na(edgepvals[n,])){
    edgepvals[n,] <- 1
  }
  if(is.na(edgetstats[n,])){
    edgetstats[n,] <- 0 
  }
  }
 
  edgetstats <- data.frame(edgetstats)
  ggplot(edgetstats, aes(t.stat)) + geom_histogram() 
  
  ##export edge files
  if(!is.null(t.stat)){
    to.bnv <- data.frame()
    for (n in 1:length(edgetstats[,1])){
      this.edge <- data.frame(edgetstats[n,])
      rownames(this.edge) <- row.names(edgetstats)[n]
      colnames(this.edge) <- "t.stat"
      if(abs(as.numeric(this.edge)) > t.stat){
        to.bnv <- rbind(to.bnv, this.edge)
      }
    }
    
    edge.file <- matrix(0, nrow = nnodes, ncol = nnodes)
    for(e in 1:length(to.bnv[,1])){
      edge.file[node, as.numeric(rownames(to.bnv)[e])] <- to.bnv[e,1]
      edge.file[as.numeric(rownames(to.bnv)[e]), node] <- to.bnv[e,1]
    }
  }
  
  write.table(edge.file, file = file.path(outputdir, paste0("V_", node, "_edge_diffs.edge.txt")), row.names = FALSE, col.names = FALSE)
}
