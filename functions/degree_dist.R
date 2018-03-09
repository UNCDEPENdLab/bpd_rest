##degree(strength) distributions script 

##BPD get single subject graphs and aggregated graph
bpd_vec <- which(subj_info$BPD == 1)
allg_noneg_bpd <- allg_noneg[bpd_vec]
allmats.bpd <- allmats[bpd_vec,,]
agg.g.bpd <- generate_agg.g(allmats.bpd)

##controls get single subject graphs and aggregated graph
control_vec <- which(subj_info$BPD == 0)
allg_noneg_control <- allg_noneg[control_vec]
allmats.control <- allmats[control_vec,,]
agg.g.control <- generate_agg.g(allmats.control)

aggregate_all <- qplot(degree_distribution(agg.g), geom = "histogram", bins = 12)
aggregate_bpd <- qplot(degree_distribution(agg.g.bpd), geom = "histogram", bins = 12)
aggregate_control <- qplot(degree_distribution(agg.g.control), geom = "histogram", bins = 12)

aggregate_all.df <- data.frame(strength = degree_distribution(agg.g))
aggregate_all <- ggplot(data = aggregate_all.df, aes(x = strength)) + geom_histogram(bins = 12)

print(aggregate_all)


pdf("figures/degree_vs_degreedistribution.pdf", width =8, height = 11)
qplot(data.frame(strength = degree_distribution(agg.g)), geom = "histogram", bins =12)
qplot(data.frame(strength = strength(agg.g)), geom= "histogram", bins = 12)
qplot(data.frame(degree = degree(agg.g)), geom= "histogram", bins = 12)
dev.off()


count.degree.distribution <- function (graph, cumulative = FALSE, ...) 
{
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  cs <- degree(graph)
  hi <- hist(cs, -1:max(cs), plot = FALSE)$count
  if (!cumulative) {
    res <- hi
  }
  else {
    res <- rev(cumsum(rev(hi)))
  }
  res
}

g <- agg.g 
degree.distribution(g)
# g.ridge <-g

qplot(strength(g))
x <- 0:max(degree(g))
qplot(x,degree.distribution(g))
qplot(degree_distribution(g.ridge), geom = "histogram", bins = 1)
qplot(degree(g))


plot_degree_distribution = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
       col = 1, main = "Degree Distribution")
}

plot_degree_distribution(agg.g)

fit_power_law = function(graph) {
  # calculate degree
  d = degree(graph, mode = "all")
  dd = degree.distribution(graph, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
       col = 1, main = "Degree Distribution")
  curve(power.law.fit, col = "red", add = T, n = length(d))
}

fit_power_law(agg.g)
####aggregate strength distributions
pdf("figures/strength.distributions.bygroup.pdf", width = 8, height = 11)
print(aggregate_all); print(aggregate_bpd); print(aggregate_control)



pdf("figures/bpd.strength.distributions.pdf", width = 8, height = 11)
for (i in 1:length(allg_noneg_bpd)){
  plot <- qplot(degree_distribution(allg_noneg_bpd[[i]]), geom = "histogram", bins = 12)
  print(plot)
} 
dev.off()


#######CONTROL strength distributions
pdf("Figures/bpd.strength.distributions.pdf", width = 8, height = 11)
for (i in 1:length(allg_noneg_bpd)){
  plot <- qplot(degree_distribution(allg_noneg_bpd[[i]]), geom = "histogram", bins = 12)
  print(plot)
} 

dev.off()

