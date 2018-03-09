pdf(file = paste0(basedir, "/figures/check_mean_FC_cor.shrink.pdf"), height = 8, width = 11)
get_meanfc_outliers(allg_density_fc)

testing <- lapply(allg_density_fc, function(subj){
  this.subj <- lapply(subj, function(density){
    edge_density(density)
    })
  message("Subj ", subj, ": ", this.subj)
})

density <- allg_density_fc[[1]][[1]]

length(allg_density_fc[[1]])

allg_density_fc <- threshold_glist(allg_noneg, rs_desired_log, method="strength", ncores=1)
for(den in 1:20){
  x <- c()
  for(subj in 1:83){
    x <- c(x, edge_density(allg_density_fc[[subj]][[den]]))
  }
  message("For threshold ", den, ", the average density is: ", mean(x))
}

subj_pearson <- allmats_pearson[1,,]
subj_shrink <- allmats
par(mfrow=c(2,1))
hist(subj_shrink)
hist(subj_pearson)
