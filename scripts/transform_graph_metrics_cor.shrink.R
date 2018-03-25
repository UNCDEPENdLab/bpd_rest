# # Transform nodal metrics for input to FA ---------------------------------

trans.file <- paste0(basedir, "/cache/threshnodalmetrics_transformed_MARCH_light_trans_", file_tag, ".RData")

library("car")
library("MASS")

# metrics.file <- paste0(basedir, "/cache/threshnodalmetrics_", file_tag, ".RData")
# message("Loading thresholded nodal statistics from file: ", metrics.file)
# load(metrics.file)

#bc for some reason, densities are not stored
#allmetrics.nodal.df$density <- rep(rep(rs_desired_log, each = 422),length(as.character(unique(allmetrics.nodal.df$id)))) 

#qplot(allmetrics.nodal.df.trans$closeness)
# qplot(allmetrics.nodal.df$part.coeff, binwidth = .01)
# BoxCox_extract((allmetrics.nodal.df$within.module.deg + 1), lambdas = seq(-10,10,.1))

allmetrics.nodal.df.trans <- allmetrics.nodal.df

allmetrics.nodal.df.trans$degree <- log(allmetrics.nodal.df$degree)
#allmetrics.nodal.df.trans$eigen.cent <- Winsorize(allmetrics.nodal.df.trans$eigen.cent, minval = quantile(allmetrics.nodal.df.trans$eigen.cent, .001), maxval = quantile(allmetrics.nodal.df.trans$degree, .999))
allmetrics.nodal.df.trans$closeness <- log(allmetrics.nodal.df$closeness)
allmetrics.nodal.df.trans$closeness <- Winsorize(allmetrics.nodal.df.trans$closeness, minval = quantile(allmetrics.nodal.df.trans$closeness, .05), maxval = quantile(allmetrics.nodal.df.trans$closeness, .99))
allmetrics.nodal.df.trans$betweenness.node <- log((allmetrics.nodal.df$betweenness.node))
#allmetrics.nodal.df.trans$page.rank <- BoxCox_extract(allmetrics.nodal.df$page.rank, lambdas = seq(-10,10,.1))
allmetrics.nodal.df.trans$within.module.deg <- sqrt(allmetrics.nodal.df$within.module.deg)
# allmetrics.nodal.df.trans$part.coeff <- sqrt(allmetrics.nodal.df$part.coeff + 1)
# allmetrics.nodal.df.trans$gateway.coeff.btw <- sqrt(allmetrics.nodal.df$gateway.coeff.btw + 1)
# allmetrics.nodal.df.trans$gateway.coeff.degree <- sqrt(allmetrics.nodal.df$gateway.coeff.degree + 1)


##just pull the metrics we are interested in and plot distributions
allmetrics_raw <- allmetrics.nodal.df.trans %>% dplyr::select(degree, eigen.cent, closeness, betweenness.node, within.module.deg,page.rank,
                                                              part.coeff, gateway.coeff.btw, gateway.coeff.degree, leverage.cent, local.clustering)

for(metric in 1:length(colnames(allmetrics_raw))){
  this.metric <- data.frame(allmetrics_raw[,metric])
  colnames(this.metric) <- colnames(allmetrics_raw)[metric]
  a <- qplot(this.metric, main = colnames(this.metric))
  print(a)
}


save(allmetrics.nodal.df.trans, file = trans.file)
