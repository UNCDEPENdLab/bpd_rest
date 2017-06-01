nodemeasures <- list.files("/Volumes/Serena/SPECC/Neil/bpd_rest/neil/stats_output_v2", pattern='^[01]',full.names=TRUE)
setwd("~/Data_Analysis/bpd_rest")
library(ggplot2)
library(reshape2)
library(lattice)


allsubjs <- c()
for (f in nodemeasures) {
  d <- read.csv(f, header=TRUE)
  d$id <- basename(f)
  allsubjs <- rbind(allsubjs, d)
}

bysubj <- split(allsubjs, f=allsubjs$id)
byroi <- split(allsubjs, f=allsubjs$ROI)

graph_measures <- names(allsubjs[!names(allsubjs) %in% c("id", "ROI")])

pdf("per_subject_histograms.pdf", width=11, height=8)
lapply(bysubj, function(subj) {
      mdata <- melt(subj[,graph_measures])
      g <- ggplot(mdata, aes(x=value)) + geom_histogram() + facet_wrap(~variable, scales="free") + ggtitle(subj$id[1L])
      plot(g)
      return(NULL)
    })
dev.off()

pdf("per_roi_histograms.pdf", width=11, height=8)
lapply(byroi, function(roi) {
      mdata <- melt(roi[,graph_measures])
      g <- ggplot(mdata, aes(x=value)) + geom_histogram() + facet_wrap(~variable, scales="free") + ggtitle(roi$ROI[1L])
      plot(g)
      return(invisible(NULL))
    })
dev.off()

pdf("per_subject_splom.pdf", width=11, height=8)
sink(file="corr_table.txt")
for (subj in bysubj) {
  cat("id: ", subj$id[1L], "\n\n")
  print(round(cor(subj[,graph_measures]), 2))
  splom(subj[,graph_measures], auto.key = list(title = subj$id[1L]))
  cat("\n\n")
}
sink()
dev.off()

str(allsubjs)

pr <- prcomp(allsubjs[,graph_measures], scale.=TRUE)
sum(pr$sdev)
cumvariance <- (cumsum((pr$sdev)^2) / sum(pr$sdev^2))
pr$rotation #varimax or promax rotation if you need to interpret this


first2 <- pr$x[,1:2]
