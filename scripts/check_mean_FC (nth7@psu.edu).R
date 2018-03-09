
# Check Mean FC between groups --------------------------------------------
source()


mean.density <- lapply(allg_noneg, function(sub){
  mean(E(sub)$weight)
})

mean.density <- data.frame(do.call(rbind, mean.density)) 
mean.density$SPECC_ID <- rownames(mean.density) 
colnames(mean.density) <- c("mean.den", "SPECC_ID")

compare.mean.density <- subj_info %>% select(SPECC_ID, BPD) 
compare.mean.density <- left_join(compare.mean.density, mean.density, by = "SPECC_ID")
t.test(mean.den~BPD, compare.mean.density)

ggplot(compare.mean.density, aes(x = BPD, y = mean.den, group = BPD)) + geom_boxplot()


mean.fc <- mean(compare.mean.density$mean.den)
sd.fc <- sd(compare.mean.density$mean.den)

nooutliers <- compare.mean.density[which(compare.mean.density$mean.den < mean.fc + 2*sd.fc),] #drops 2 outliers with especially high FC (in ridge at least): one in BPD group, one in HC group
#nooutliers <- nooutliers[which(nooutliers$mean.den > mean.fc - sd.fc),] 

ggplot(nooutliers, aes(x = BPD, y = mean.den, group = BPD)) + geom_boxplot()
t.test(mean.den~BPD, nooutliers)

subjs_outliers <- which(compare.mean.density$mean.den > mean.fc + 2*sd.fc)

