####diagnose SCOTMI
##grab 2 SCOTMI files and correlate lower triangles
setwd("~/Box Sync/RS_BPD_graph")
basedir <- getwd()
adjbase.scotmi <- file.path(paste0(basedir,"/adjmats_269_scotmi/"))
adjbase.pearson <- file.path(paste0(basedir,"/adjmats_269_pearson/"))

#s1.s <- as.matrix(read.csv(as.character(paste0(adjbase.scotmi,"/005ai_06Nov2013_S-raw.txt"))))
#s1.p <- as.matrix(read.table(as.character(paste0(adjbase.pearson,"/005ai_06Nov2013_adjmat269_pearson_fisherz.txt.gz"))))
s1.s <- as.matrix(read.csv(as.character(paste0(adjbase.scotmi,"/001ra_07dec2013_S-raw.txt"))))
s1.p <- as.matrix(read.table(as.character(paste0(adjbase.pearson,"/001ra_07dec2013_adjmat269_pearson_fisherz.txt.gz"))))

s1.s <- rbind(array(NaN, ncol(s1.s)), s1.s)
s1.s[upper.tri(s1.s)] <- t(s1.s)[upper.tri(s1.s)]

svec <- s1.s[lower.tri(s1.s)]
pvec <- s1.p[lower.tri(s1.p)]

par(mfrow=c(1,2))
#hist(log(svec + .05), main="SCoTMI"); hist(pvec, main="Pearson")
hist(sqrt(svec), main="SCoTMI"); hist(pvec, main="Pearson")

cor.subs <- cor(s1.s[lower.tri(s1.s)], s1.p[lower.tri(s1.p)])


cor(sqrt(svec), pvec)
cor(log(svec + .05), pvec)

#check for mismatch in patterns
#maxcorr <- -Inf
cmatch <- matrix(NA, nrow=nrow(s1.s), ncol=ncol(s1.s))
minervamatch <- matrix(NA, nrow=nrow(s1.s), ncol=ncol(s1.s))
plotgrid <- vector('list', 269^2)
dim(plotgrid) <- c(269, 269)
library(ggplot2)
for (i in 1:10) { #nrow(s1.s)) {
  for (j in 1:10) { #ncol(s1.s)) {
    library(robust)
    # for (i in 1:nrow(s1.s)) {
    # for (j in 1:ncol(s1.s)) {
    #browser()
    #plot(s1.s[i,-i], s1.p[j,-j])
    v1 <- s1.s[i,] #scotmi row
    v2 <- abs(s1.p[j,]) #pearson row
    v1[c(i,j)] <- v2[c(i,j)] <- NA #remove any values in the nodes with themselves
    
    df <- data.frame(scotmi = v1, pearson=v2, i=i, j=j)
    g <- ggplot(df, aes(x=scotmi, y=pearson)) + geom_point() + theme_bw() + ggtitle(paste("SCoTMI", i, "with Pearson", j))
    plotgrid[[i,j]] <- g
    
    
    cmatch[i,j] <- cor(v1, v2, use="pairwise.complete.obs", method="spearman") #SCoTMI with absolute value of correlation
    #cmatch[i,j] <- cor(log(v1) + 0.05), abs(v2))
    
    #cmatch[i,j] <- robust::covRob(cbind(v1, v2), corr=TRUE, estim="donostah", na.action=na.exclude)$cov[1,2]
    #minervamatch[i,j] <- minerva::mine(s1.s[i,-i], s1.p[j,-j])$MIC
  }
}

par(mfrow = c(1,2))
hist(apply(cmatch[1:10,1:10], 1, max))

hist(diag(cmatch))


# hist(apply(minervamatch, 1, max))

library(cowplot)
plotgrid <- plotgrid[1:10,1:10]
pdf("10nodes_scotmivspearson.pdf", width=30, height=30)
do.call(plot_grid, plotgrid)
dev.off()

?cowplot::plot_grid
