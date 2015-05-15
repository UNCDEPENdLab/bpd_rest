nodemeasures <- list.files("/Volumes/Serena/SPECC/Neil/bpd_rest/neil/stats_output_v2/", pattern='^[01]',full.names=TRUE)
library(ggplot2)
library(reshape2)
library(lattice)
library(reshape2)

nodemeasures = nodemeasures[1:2]

allsubjs <- c()
for (f in nodemeasures) {
  d <- read.csv(f, header=TRUE)
  d_wide <- dcast(d,roi ~ edge_def+parameter+statistic, value.var='value')
  d_wide$id <- basename(f)
  allsubjs <- rbind(allsubjs, d_wide)
}

bysubj <- split(allsubjs, f=allsubjs$id)
byroi <- split(allsubjs, f=allsubjs$roi)

graph_measures <- grep('HARD_0.*local_eff',names(allsubjs),value=TRUE)
