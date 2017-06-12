#####################################
setwd("~/Box Sync/RS_BPD_graph")

allmetrics <- readRDS("allmetrics.pearson.FDremove.rds")
length(allmetrics)

library(ggplot2)
SPECC_rest
ss.metric <- array(NA, length(allmetrics), dimnames = list("BPD" = SPECC_rest[,"BPD"]))
#####################

##misc code from plotting main effects of BPD and bpd X age interactions


###[[i]][[12]] for all subjects i at density 12.
##here you need to specify what metric you'd like to compute and for what node to compute

####plot main effect of BPD differences and degree

for (i in 1:length(allmetrics)){
  ss.metric[i] <- allmetrics[[i]][[11]]$degree["V64"]
}

ss.metric.1 <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric=ss.metric)
ss.metric.1$BPD <- factor(ss.metric.1$BPD, levels=c(0,1), labels=c("Control", "BPD"))
ss.metric.1$roiname <- "L Rolandic Operculum" 

for (i in 1:length(allmetrics)){
  ss.metric[i] <- allmetrics[[i]][[12]]$degree["V202"]
}

ss.metric.2 <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric=ss.metric)
ss.metric.2$BPD <- factor(ss.metric.2$BPD, levels=c(0,1), labels=c("Control", "BPD"))
ss.metric.2$roiname <- "L Superior Medial Gyrus" 

for (i in 1:length(allmetrics)){
  ss.metric[i] <- allmetrics[[i]][[12]]$degree["V233"]
}

ss.metric.3 <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric=ss.metric)
ss.metric.3$BPD <- factor(ss.metric.3$BPD, levels=c(0,1), labels=c("Control", "BPD"))
ss.metric.3$roiname <- "R Caudate"

for (i in 1:length(allmetrics)){
  ss.metric[i] <- allmetrics[[i]][[12]]$degree["V216"]
}

ss.metric.4 <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric=ss.metric)
ss.metric.4$BPD <- factor(ss.metric.4$BPD, levels=c(0,1), labels=c("Control", "BPD"))
ss.metric.4$roiname <- "R Middle Cingulate"

for (i in 1:length(allmetrics)){
  ss.metric[i] <- allmetrics[[i]][[12]]$degree["V193"]
}

ss.metric.5 <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric=ss.metric)
ss.metric.5$BPD <- factor(ss.metric.5$BPD, levels=c(0,1), labels=c("Control", "BPD"))
ss.metric.5$roiname <- "R Middle Frontal Gyrus"

for (i in 1:length(allmetrics)){
  ss.metric[i] <- allmetrics[[i]][[12]]$degree["V231"]
}

ss.metric.6 <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric=ss.metric)
ss.metric.6$BPD <- factor(ss.metric.6$BPD, levels=c(0,1), labels=c("Control", "BPD"))
ss.metric.6$roiname <- "R Putamen" 

for (i in 1:length(allmetrics)){
  ss.metric[i] <- allmetrics[[i]][[12]]$degree["V238"]
}

ss.metric.7 <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric=ss.metric)
ss.metric.7$BPD <- factor(ss.metric.7$BPD, levels=c(0,1), labels=c("Control", "BPD"))
ss.metric.7$roiname <- "R Superior Temporal Gyrus" 

ss.all <- rbind(ss.metric.2, ss.metric.3, ss.metric.4, ss.metric.5, ss.metric.6, ss.metric.7)



dev.off()

pdf("degree_plot.pdf", width=10, height=7)
g <- ggplot(ss.all, aes(x = factor(BPD), y = metric, color = BPD)) + stat_summary(fun.data="mean_cl_boot", size=1.5, fatten=1.5) + theme(legend.position="none") +   #fill = factor(BPD)
  labs(x = "", y = "degree centrality", title = "") + scale_x_discrete(breaks=c(0,1), labels=c("Control", "BPD"))  + facet_wrap(~roiname, scales="free_y")  +
  theme_grey(base_size = 16) + theme(legend.title=element_blank(), strip.text = element_text(size=16)) + theme(legend.position="bottom")
plot(g)
dev.off()
SPECC_rest
# The data structure should look like this:
# age   BPD metric roiname
# dsfdf 0   16     mfg
# sdf   1          mfg
# df
# d
# df
# dsfdf 0   16     amygdala
# sdf   1          amygdala
# df
# d
# df



ggplot(ss.metric.df, aes(x = factor(BPD), y = metric, fill = factor(BPD))) +geom_boxplot() + theme(legend.position="none") + labs(x = "", y = "degree centrality", title = "")

ggplot(ss.metric.df, aes(x = factor(BPD), y = metric)) + stat_summary(fun.data="mean_cl_boot", size=1.5, fatten=2) + theme(legend.position="none") +  # fill = factor(BPD))
  labs(x = "", y = "degree centrality", title = "") + scale_x_discrete(breaks=c(0,1), labels=c("Control", "BPD")) + theme_bw(base_size=14)

####################################################################
#####interactions of interest
## V109  L Middle Orbital Gyrus
ss.metric.ixn1 <- array(NA, length(allmetrics), dimnames = list("BPD" = SPECC_rest[,"BPD"]))
for (i in 1:length(allmetrics)){
  ss.metric[i] <- allmetrics[[i]][[10]]$part.coeff[109]
}

ss.metric.ixn1 <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric=ss.metric)
ss.metric.ixn1$BPD <- factor(ss.metric.ixn1$BPD, levels=c(0,1), labels=c("Control", "BPD"))
ss.metric.ixn1$roiname <- "L Middle Orbital Gyrus" 

## V 
ss.metric <- array(NA, length(allmetrics), dimnames = list("BPD" = SPECC_rest[,"BPD"]))
for (i in 1:length(allmetrics)){
  ss.metric[i] <- allmetrics[[i]][[10]]$part.coeff[230]
}

ss.metric.2 <- data.frame(SPECC_rest[,c("BPD", "AgeAtScan")], metric=ss.metric)
ss.metric.2$BPD <- factor(ss.metric.2$BPD, levels=c(0,1), labels=c("Control", "BPD"))
ss.metric.2$roiname <- "L Middle Orbital Gyrus" 


#interaction
ggplot(ss.metric.2, aes(x = AgeAtScan, color=BPD, y = metric)) + geom_point() + stat_smooth(method="lm", se=FALSE)+ 
  labs(x = "Age", y = "participation coefficient", title = "R Putamen") +scale_color_brewer("Group", palette="Set1") + theme(plot.title = element_text(hjust = 0.5))

+ theme(legend.position="none")

ggplot(ss.metric.df, aes(x = AgeAtScan, color=BPD, y = metric)) + geom_point() + stat_smooth(method="lm", se=FALSE) + labs(x = "Age", y = "", title = "") +
  scale_color_brewer("Group", palette="Set1") + theme(legend.position="none")



library(cowplot)

cowplot::plot_grid(g1, g2)
