pcaout <- pca(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 3)#, scores = "Bartlett")
cor(dplyr::select(metrics.raw_fa, -id, -node))
vv <- cor(dplyr::select(metrics.raw_fa, -id, -node))

head(vv)
vv[lower.tri(vv)]
vv[upper.tri(vv)]

lapply(dplyr::select(metrics.raw_fa, -id, -node), is.na)
lapply(dplyr::select(metrics.raw_fa, -id, -node), function(x) { sum(is.na(x)) })
str(dplyr::select(metrics.raw_fa, -id, -node))
names(metrics.raw_fa)
misspos <- lapply(dplyr::select(metrics.raw_fa, -id, -node), function(x) { which(is.na(x)) })
misspos
metrics.raw_fa$id[misspos$`0.02_part.coeff`]
length(unique(metrics.raw_fa$id[misspos$`0.02_part.coeff`]))
vv <- cor(dplyr::select(metrics.raw_fa, -id, -node), use="pairwise")
vv[lower.tri(vv)]
range(vv[lower.tri(vv)])
corrvec <- vv[lower.tri(vv)]
hugcorr <- corrvec[corrvec > .96]
hugcorr
badness = which(vv > .96, arr.ind=TRUE)
badness
badness = which(vv > .98, arr.ind=TRUE)
str(badness)
badness
str(badness)
which(vv == 1, arr.ind=T)
badpos <- which(vv == 1, arr.ind=T)
badpos <- badpos[badpos[,1] != badpos[,2]]
badpos
str(badpos)
badpos <- which(vv == 1, arr.ind=T)
str(vv)
badpos <- which(vv == 1, arr.ind=T)
str(badpos)
badpos
badpos <- badpos[badpos[,1] != badpos[,2]]
badpos
str(badpos)
badpos <- which(vv == 1, arr.ind=T)
badpos
str(badpos)
badpos <- as.data.frame(badpos)
str(badpos)
filter(badpos, row!=col)
filter(badpos, row!=col)
str(vv)
filter(badpos, row!=col)
dimnames(vv)[[1]][filter(badpos, row!=col)]
dimnames(vv)
dimnames(vv)[[1]]
dimnames(vv)[[1]][as.matrix(filter(badpos, row!=col))]
badpos
badpos <- which(vv == 1, arr.ind=T)
