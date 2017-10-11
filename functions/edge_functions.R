run_edge_comparisons <- function(subj_info){

controls <- which(subj_info$BPD == 0)
bpd <- which(subj_info$BPD == 1)

edgepvals <- matrix(NA, nrow=nnodes, ncol=nnodes)
abst <- 3.2 #about p = .0001 (two-tailed)
edge.diffs <- data.frame()
for (i in 1:nnodes) {
  for (j in 1:nnodes) {
    if (i > j) {
      #lower triangle
      controlcorrs <- allmats[controls,i,j]
      bpdcorrs <- allmats[bpd, i, j]
      res <- t.test(atanh(controlcorrs), atanh(bpdcorrs))
      edgepvals[i,j] <- res$p.value
      edgepvals[j,i] <- res$statistic #t stat on upper triangle
      if (!is.nan(res$statistic) && abs(res$statistic) > abst) {
        cat("Significant edge difference: ", i, j, ", t=", round(res$statistic, 3), "p =", round(res$p.value, 3), "means are:", round(res$estimate, 5), "\n")
        sig.edge <- data.frame(edge1 = i, edge2 = j, tval = round(res$statistic, 3), p = round(res$p.value, 3), mean.control = round(res$estimate[1], 5), mean.bpd = round(res$estimate[2], 5))
        edge.diffs <- rbind(edge.diffs, sig.edge)
        
      }
    }
  }
}
for (i in 1:length(edge.diffs[,1])){
  if(edge.diffs$mean.control[i] > edge.diffs$mean.bpd[i]){
    edge.diffs$control.bpd[i] <- .5
  } else {edge.diffs$control.bpd[i] <- 1}
}

row.names(edge.diffs) <- NULL
return(edge.diffs)
}



# ####IN PROGRESS: FDR correction for multiple comparisons
# pvec <- edgepvals[lower.tri(edgepvals)]
# #tvec <- edgepvals[upper.tri(edgepvals)]
# lookup <- which(lower.tri(edgepvals) == TRUE, arr.ind=TRUE)
# padj <- p.adjust(pvec, "fdr")
# lookup[which(pvec < .0001),]
# padj[which(padj < .05)]
# lookup[which(padj < .05),]


edge_bnv_files <- function(edge.comp, outputdir){
  # browser()
  nodes1 <- as.matrix(unique(subset(edge.comp, edge1 %in% seq(1,269,1), select = edge1)))
  nodes2 <- as.matrix(unique(subset(edge.comp, edge2 %in% seq(1,269,1), select = edge2)))
  # nodes.involved <- paste0("V", sort(unique(rbind(nodes1, nodes2))))
  # nnodes.involved <- length(nodes.involved)

  nodes.involved <- sort(unique(rbind(nodes1, nodes2)))

  nodevals <- rep(0, nnodes)
  nodevals[nodes.involved] <- 1
 
  
  edge.diffs <- NodeFile(
    atlas = atlas,
    nnodes = nnodes,
    # nodestp = nodes.involved,
    nodevals = nodevals,
    labels = 0,
    filename = "edge.diffs",
    outputdir = outputdir
  )

  # edge.file <- matrix(0, nrow=nnodes.involved, ncol=nnodes.involved)
  edge.file <- matrix(0, nrow=nnodes, ncol=nnodes)
  # browser()
  for (sig.e in 1:length(edge.comp[,1])){
    e1 <- edge.comp$edge1[sig.e]
    e2 <- edge.comp$edge2[sig.e]
  edge.file[e1,e2] <- edge.comp$control.bpd[sig.e]
  edge.file[e2,e1] <- edge.comp$control.bpd[sig.e]
  }
  
  # edge.file <- matrix(0, nrow=nnodes, ncol=nnodes) 
  # 
  # # for (i in 1:nnodes) {
  # #   for (j in 1:nnodes) {
  # #     #populate lower triangle only for speed, then transpose onto upper triangle
  # #     if (i %in% nodes.involved) { 
  # #       # roiDist[i,j] <- sqrt((roiMat[i,1] - roiMat[j,1])^2 + (roiMat[i,2] - roiMat[j,2])^2 + (roiMat[i,3] - roiMat[j,3])^2) 
  # #       
  # #       }
  # #   }
  # # }
  # # 
  # #copy lower triangle to upper triangle
  # roiDist[upper.tri(roiDist)] <- t(roiDist)[upper.tri(roiDist)]
  # diag(roiDist) <- 0
  # 
  # edge.file <- data.frame(edge.file)
  write.table(edge.file, file = file.path(outputdir, "edge.diffs.edge.txt"), row.names = FALSE, col.names = FALSE)

  outputlist <- list()
  outputlist[["nodes"]] <- nodes.involved
  outputlist[["node.file"]] <- edge.diffs
  outputlist[["edge.list"]] <- edge.file
  return(outputlist)
  }
