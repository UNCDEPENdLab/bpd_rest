
# Final Pipeline comparison----------------------------------------------------------

setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()

pipelines <- list.files(pattern="schaefer.*", path="adjmats", full.names=FALSE)

massive_df <- list()

for (p in pipelines) {
  plist <- list()
  adjmats <- list.files(path=file.path("adjmats", p), pattern=".*\\.txt(\\.gz)*")
  if (length(adjmats) == 0L) { next }
  ids <- sub("(.*)_schaefer422.*", "\\1", adjmats, perl=TRUE)
  for (a in 1:length(adjmats)) {
    df <- read.table(file.path("adjmats", p, adjmats[a]))
    plist[[ ids[a] ]] <- as.matrix(df)
  }
  
  massive_df[[ p ]] <- plist
}

#pairwise combinations of pipelines
pcomb <- combn(length(massive_df), 2)

#pnames <- abbreviate(names(massive_df), minlength=6, strict=FALSE)
pnames <- sub("schaefer422_", "", names(massive_df), fixed=TRUE)

dumsum <- summary(c(NA, rnorm(10)))
m_mat <- array(NA, dim=c(length(massive_df), length(massive_df), length(dumsum)),
                dimnames=list(p1=pnames, p2=pnames, stat=names(dumsum)))

library(robust)

for(p in 1:ncol(pcomb)) {
  p1 <- massive_df[[ pcomb[1,p] ]]
  p2 <- massive_df[[ pcomb[2,p] ]]
  cvec <- rep(NA, length(p1))
  for (i in 1:length(p1)) {
    id_target <- names(p1)[i]
    adj1 <- p1[[id_target]]
    adj2 <- p2[[id_target]]
    if (!(is.null(adj1) || is.null(adj2))) {
      cvec[i] <- cor(adj1[lower.tri(adj1)], adj2[lower.tri(adj2)])
      #cvec[i] <- covRob(cbind(adj1[lower.tri(adj1)], adj2[lower.tri(adj2)]), corr=TRUE, na.action=na.exclude)$cov[1,2] #robust to outliers
    }
  }
  sv <- summary(cvec)
  m_mat[ pcomb[1,p], pcomb[2,p], 1:length(sv)] <- sv
}

#reflect to lower triangle
m_mat2 <- plyr::aaply(m_mat, 3, function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
})

m_mat2 <- aperm(m_mat2, c(2,3,1))

sink(file="pipeline_corr_pearson.txt")
options(width=400)
for (i in 1:dim(m_mat2)[3]) {
  cat("\n\nMETRIC: ", dimnames(m_mat2)[[3]][i])
  print(round(m_mat2[,,i], 3))  
}
options(width=120)
sink()


#dens clime is out
#aroma_pcor.shrink_partial: BP + AROMA using unsmoothed data
#aroma_pcor.shrink_partial_lambda.01: BP + AROMA prob with undo smooth
#aroma_pearson: pearson BP + AROMA using unsmoothed
#aroma_ridge.net_partial: HP + AROMA using unsmoothed
#aroma_ridge.net_partial_bp: BP + AROMA using unsmoothed
#aroma_nosmooth_aroma_bp_nonaggr_pcor.shrink_partial: EMPTY nonaggressive AROMA
#aroma_nosmooth_aroma_bp_nonaggr_ridge.net_partial: BP + nonaggressive AROMA using undosmooth
#aroma_nosmooth_aroma_bp_pcor.shrink_partial: BP + aggressive AROMA with undosmooth
#nosmooth_aroma_bp_ridge.net_partial: AUG2017 PAUSE: BP + aggressive AROMA with undosmooth
