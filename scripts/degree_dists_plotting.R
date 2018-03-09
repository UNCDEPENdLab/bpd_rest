# looking at degree distributions across thresholds between groups-----------------------

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

agg.g.bpd_fcthresh <- list()
for(i in 1:length(rs_desired_log)){
  agg.g.bpd_fcthresh[[i]] <- threshold_graph(agg.g.bpd, rs_desired_log[i], method="strength")#, ncores=1)
}

agg.g.controls_fcthresh <- list()
for(i in 1:length(rs_desired_log)){
  agg.g.controls_fcthresh[[i]] <- threshold_graph(agg.g.controls, rs_desired_log[i], method="strength")#, ncores=1)
}

agg.g.groups <- list(bpd = agg.g.bpd_fcthresh, controls = agg.g.controls_fcthresh)

agg.g.groups["bpd"]
names(agg.g.groups)
pdf("degree_dists_groups_.pdf", width = 15, height = 10)
plot_file <- data.frame()
descriptives_distributions <- data.frame()

for(i in names(agg.g.groups)){
  group <- agg.g.groups[i]
  
  group_degree <- array(NA, c(nnodes, length(rs_desired_log)))
  colnames(group_degree) <- rs_desired_log
  row.names(group_degree) <- paste0("V", seq(1,422,1))
  
  group_descriptives <- array(NA,c(length(rs_desired_log), 8))
  group_descriptives <- data.frame(group_descriptives)
  colnames(group_descriptives) <- c("mean", "SD", "mode", "median", "range_low", "range_high", "FC_thresh", "group")
  
  for(r in 1: length(rs_desired_log)){
    group_degree[,r] <- degree(group[[1]][[r]])
    group_descriptives$mean[r] <- mean(group_degree[,r]) 
    group_descriptives$SD[r] <- sd(group_degree[,r]) 
    group_descriptives$mode[r] <- Mode(group_degree[,r]) 
    group_descriptives$median[r] <- median(group_degree[,r]) 
    group_descriptives$range_low[r] <- range(group_degree[,r])[1]
    group_descriptives$range_high[r] <- range(group_degree[,r])[2]
    group_descriptives$FC_thresh[r] <- rs_desired_log[r]
    group_descriptives$group <- i
  }
  
  
  
  descriptives_distributions <- rbind(descriptives_distributions,group_descriptives)
  
  group_melt <- data.frame(melt(group_degree)) %>% dplyr::select(-Var1); colnames(group_melt) <- c("r_val", "degree")
  group_melt$group <- i
  
  plot_file <- rbind(plot_file, group_melt)
}


degree_hists <- ggplot(data = plot_file, aes(x = degree)) + 
  geom_histogram(data=subset(plot_file, group == 'bpd'), aes(fill=group), alpha = .4) + 
  geom_histogram(data=subset(plot_file, group == 'controls'), aes(fill=group), alpha = .4) + 
  facet_wrap(~r_val) + scale_fill_manual(name = "group", values = c("red", "blue"), labels = c("BPD", "Controls")) +
  theme_bw() +
  labs(title = "Degree distributions across FC thresholds based on partial correlations estimated via ridge regression", y = "Density", x = "Degree") +
  theme(plot.title = element_text(hjust = 0.5))
print(degree_hists)

descriptives_distributions_plot <- descriptives_distributions %>% dplyr::select(-range_low, -range_high)


melted_descriptives <- melt(descriptives_distributions_plot[,c(1,3,4,6)])
melted_descriptives$FC_thresh <- rep(rs_desired_log, 6)
melted_descriptives$SD <- rep(descriptives_distributions_plot$SD)

describe_lineplot <- ggplot(melted_descriptives, aes(x = FC_thresh, color = group, y = value)) +
  geom_smooth() +facet_grid(~variable) +theme_bw() +labs(title= "Descriptive Statisitcs ")

print(describe_lineplot)
rev(rs_desired_log)
