# Euclidean Distance estimation -------------------------------------------

#estimate ROI-ROI Euclidean distance (use cache if available since the double for loop takes ~4 seconds)
expectFile <- file.path(basedir, "cache", paste0("distmat_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
if (file.exists(expectFile)) {
  load(expectFile)
} else {
  roiDist <- matrix(NA, nrow=nnodes, ncol=nnodes) 
  
  system.time(for (i in 1:nnodes) {
    for (j in 1:nnodes) {
      #populate lower triangle only for speed, then transpose onto upper triangle
      if (i > j) { roiDist[i,j] <- sqrt((roiMat[i,1] - roiMat[j,1])^2 + (roiMat[i,2] - roiMat[j,2])^2 + (roiMat[i,3] - roiMat[j,3])^2) }
    }
  })
  
  #copy lower triangle to upper triangle
  roiDist[upper.tri(roiDist)] <- t(roiDist)[upper.tri(roiDist)]
  diag(roiDist) <- 0
  
  ####quick QA
  # hist(roiDist)
  # vecDist <- roiDist[lower.tri(roiDist)]
  # sum(vecDist < 20)/length(vecDist)
  
  #roi.dist can be changed at the front end of the script
  rmShort <- roiDist > roi.dist
  #creates binary matrix, in which 0 denotes a short distanced connection that is to be removed. 
  rmShort <- apply(rmShort, c(1,2), as.numeric)
  save(file=expectFile, rmShort, roiDist)
}
