
# generate degree distributions --------------------------------------------



##plot single subject degree distributions across thresholds

pdf("single_subjs_degree_dist_binary.pdf", width = 8, height = 11)
for(i in 1:length(allg_density_fc)){
  subj <- allg_density_fc[i]
all_degree <- array(NA, c(nnodes, length(rs_desired_log)))
colnames(all_degree) <- rs_desired_log
row.names(all_degree) <- paste0("V", seq(1,422,1))
  for(r in 1:length(rs_desired_log)){
    this_r_deg <- degree(subj[[1]][[r]]) 
    all_degree[,r] <- this_r_deg
  }

all_degree.melt <- data.frame(melt(all_degree));colnames(all_degree.melt) <- c("node", "r_val", "degree")
x <- ggplot(data = all_degree.melt, aes(x = degree)) + geom_histogram(bins = 12) + facet_wrap(~r_val)
print(x)
}

dev.off()