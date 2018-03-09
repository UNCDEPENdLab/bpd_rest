# # Transform nodal metrics for input to FA ---------------------------------
# 
# ##to initialize globals:
# ##########
# setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
# #initialize the graph analysis pipeline (includes sourcing pipeline inputs, creating graphs, thresholding, binarization, community assignment, and calculation and reduction of global and nodal graph metrics)
# 
# source("scripts/setup_globals.R")
# #RIDGE
# inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, preproc_pipeline = "nosmooth_aroma_bp_nonaggr") #leave blank for defaults, consult with function for deviations from default
# #PEARSON
# #inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, conn_method = "pearson", rs_desired_log = logspace(log10(.2), log10(.4), 20)) #leave blank for defaults, consult with function for deviations from default
# for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements
# ##########
# 

# trans.file <- paste0(basedir, "/cache/threshnodalmetrics_transformed_", file_tag, ".RData")
# if(file.exists(trans.file)){
#   load(trans.file)
#   message("Loaded transformed graph metrics from: ", trans.file )
# } else{
library("car")
library("MASS")

metrics.file <- paste0(basedir, "/cache/threshnodalmetrics_", file_tag, ".RData")
message("Loading thresholded nodal statistics from file: ", metrics.file)
load(metrics.file)

#bc for some reason, densities are not stored
allmetrics.nodal.df$density <- rep(rep(rs_desired_log, each = 422),length(as.character(unique(allmetrics.nodal.df$id)))) 



allmetrics.nodal.df.trans <- allmetrics.nodal.df

allmetrics.nodal.df.trans$degree <- Winsorize(allmetrics.nodal.df.trans$degree, minval = quantile(allmetrics.nodal.df.trans$degree, .001), maxval = quantile(allmetrics.nodal.df.trans$degree, .999))
allmetrics.nodal.df.trans$eigen.cent <- Winsorize(allmetrics.nodal.df.trans$eigen.cent, minval = quantile(allmetrics.nodal.df.trans$eigen.cent, .001), maxval = quantile(allmetrics.nodal.df.trans$degree, .999))
allmetrics.nodal.df.trans$betweenness.node <- Winsorize(log(allmetrics.nodal.df.trans$betweenness.node), minval = quantile(log(allmetrics.nodal.df.trans$betweenness.node), .01), maxval = quantile(log(allmetrics.nodal.df.trans$betweenness.node), .99))
allmetrics.nodal.df.trans$page.rank <- Winsorize(allmetrics.nodal.df.trans$page.rank, minval = quantile(allmetrics.nodal.df.trans$page.rank, .001), maxval = quantile(allmetrics.nodal.df.trans$page.rank, .999))
allmetrics.nodal.df.trans$leverage.cent <- Winsorize(allmetrics.nodal.df.trans$leverage.cent, minval = quantile(allmetrics.nodal.df.trans$leverage.cent, .005, na.rm = TRUE), maxval = quantile(allmetrics.nodal.df.trans$leverage.cent, .999, na.rm = TRUE))
#within-mod is okay for now


###Use box-cox transform on the integration metrics

##tiny function to test a range of lambdas, extract the lambda leading to the highest log-likelihood and transform the data
BoxCox_extract <- function(data, lambdas, plot = TRUE){
  Box <- boxcox(data ~ 1, lambda = lambdas)
  Cox = data.frame(Box$x, Box$y)            # Create a data frame with the results
  
  Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] # Order the new data frame by decreasing y
  Cox2[1,]    #display lambda with highest log-likelihood
  lambda = Cox2[1, "Box.x"] #extract lambda
  
  transformed <- (data ^lambda -1)/lambda
  
  if(plot) {
  before <- qplot(data)
  print(before)
  
  x <- qplot(transformed)
  print(x)
  }
  return(transformed)
}

part.coeff <- Winsorize(allmetrics.nodal.df.trans$part.coeff, minval = quantile(allmetrics.nodal.df.trans$part.coeff, .005, na.rm = TRUE), maxval = quantile(allmetrics.nodal.df.trans$part.coeff, .999, na.rm = TRUE))

allmetrics.nodal.df.trans$part.coeff <- BoxCox_extract(part.coeff, seq(5,20,.1))

##gateway coefficient
gateway.btw <- Winsorize(allmetrics.nodal.df.trans$gateway.coeff.btw, minval = quantile(allmetrics.nodal.df.trans$gateway.coeff.btw, .005, na.rm = TRUE), maxval = quantile(allmetrics.nodal.df.trans$gateway.coeff.btw, .999, na.rm = TRUE))

allmetrics.nodal.df.trans$gateway.coeff.btw <- BoxCox_extract(gateway.btw, seq(-3,25,.1))

#based on degree
gateway.degree <- Winsorize(allmetrics.nodal.df.trans$gateway.coeff.degree, minval = quantile(allmetrics.nodal.df.trans$gateway.coeff.degree, .005, na.rm = TRUE), maxval = quantile(allmetrics.nodal.df.trans$gateway.coeff.degree, .999, na.rm = TRUE))

allmetrics.nodal.df.trans$gateway.coeff.degree <- BoxCox_extract(gateway.degree, seq(-3,25,.1))

##just pull the metrics we are interested in and plot distributions
allmetrics_raw <- allmetrics.nodal.df.trans %>% dplyr::select(degree, eigen.cent, closeness, betweenness.node, page.rank, leverage.cent, within.module.deg,
                                                        part.coeff, gateway.coeff.btw, gateway.coeff.degree)

for(metric in 1:length(colnames(allmetrics_raw))){
  this.metric <- data.frame(allmetrics_raw[,metric])
  colnames(this.metric) <- colnames(allmetrics_raw)[metric]
  a <- qplot(this.metric, main = colnames(this.metric))
  print(a)
}


save(allmetrics.nodal.df.trans, file = trans.file)
# 
# message("Saved transformed graph metrics to: ", trans.file )
# }