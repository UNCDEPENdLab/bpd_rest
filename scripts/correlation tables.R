setwd("~/Box Sync/RS_BPD_graph")
basedir <- getwd()

source("Graph_util_redux.R")
library(igraph)
library(RGtk2)
library(brainGraph)
library(ggplot2)
library(dplyr)
library(gdata)
library(beepr)
library(parallel)
library(data.table)
library(psych)

#computed nodal metrics
if(file.exists(paste0(basedir, "/cachedRfiles/allmetrics.binary.pearson.RData")) == TRUE) {
  allmetrics.pearson  <- get(load(paste0(basedir, "/cachedRfiles/allmetrics.binary.pearson.RData")))
} 
if(file.exists(paste0(basedir, "/cachedRfiles/allmetrics.binary.scotmi.RData")) == TRUE) {
  allmetrics.scotmi  <- get(load(paste0(basedir, "/cachedRfiles/allmetrics.binary.scotmi.RData")))
} 
if(file.exists(paste0(basedir, "/cachedRfiles/allg.weighted.pearson.RData")) == TRUE) {
  allg.weighted.pearson  <- get(load(paste0(basedir, "/cachedRfiles/allg.weighted.pearson.RData")))
} 
if(file.exists(paste0(basedir, "/cachedRfiles/allg.weighted.pearson.RData")) == TRUE) {
  allg.weighted.scotmi  <- get(load(paste0(basedir, "/cachedRfiles/allg.weighted.scotmi.RData")))
} 
rm(allg, allmetrics)


atlas <- read.csv("power269_masterlookup_shift_nate.csv", header = TRUE)
SPECC_rest <- read.csv("SPECC_info_trimmed.csv", header = TRUE)
SPECC_rest$SPECC_ID <- as.character(SPECC_rest$SPECC_ID)
nnodes <- 269

################################################################################
#####filter subjects with over .15 brain volumes displaces .5mm or more
table(SPECC_rest[,c(3,5)])
SPECC_rest <- filter(SPECC_rest, pr_over5mm <= .15)
table(SPECC_rest[,c(3,5)])

describe(SPECC_rest[,c(1:6, 8)])

################################################################################
##create table of nodal metrics by subject and by network across densities
node.table.pearson <- data.frame()

names(allmetrics.pearson) <- SPECC_rest$SPECC_ID #assumes elements of allmetrics in same order as this DF
attr(allmetrics.pearson, "bpd") <- SPECC_rest[,"BPD"]
attr(allmetrics.pearson, "age") <- SPECC_rest[,"AgeAtScan"]

node.table.pearson <- do.call(rbind, lapply(1:length(allmetrics.pearson), function (s) {
  subj <- allmetrics.pearson[[s]]
  
  ddf <- lapply(1:length(subj), function(d) {
    density <- subj[[d]]
    
    metrics_to_grab <- names(density)[!names(density) %in% c("origgraph")]
    nvec <- lapply(1:length(metrics_to_grab), function(m) { 
      density[[metrics_to_grab[m]]]
    })
    
    names(nvec) <- metrics_to_grab
    ndf <- as.data.frame(do.call(cbind, nvec))
    ndf$density <- names(subj)[d]
    ndf$node <- atlas$vname
    ndf$anat_label <- atlas$anat_label
    ndf$scotmi <- 0
    return(ndf)
  })
  
  #names(ddf) <- names(subj)
  
  
  ddf <- do.call(rbind, ddf)
  ddf$subj <- names(allmetrics.pearson)[s]
  ddf$bpd <- attr(allmetrics.pearson, "bpd")[s]
  ddf$age <- attr(allmetrics.pearson, "age")[s]
  
  return(ddf)
  
}))

metrics <- names(node.table.pearson[, c(2:6, 8:12)])
densities_desired <- seq(.01, .2, .01) 
nodes <- 1:269
nsubj <- 79



###########output from node.metrics.pearson is a list of 20 densities with 269 nodal t-tests (including group means and t values) by all 10 nodal metrics (excluding local.clsutering)
if(file.exists(paste0(basedir, "/cachedRfiles/node.metrics.pearson.corr.RData")) == TRUE){
  node.metrics.pearson  <- get(load(paste0(basedir, "/cachedRfiles/node.metrics.pearson.corr.RData")))
} else{
  
library(foreach)
library(doSNOW)
setDefaultClusterOptions(master="localhost") #move away from 10187 to avoid collisions
clusterobj <- makeSOCKcluster(4)
registerDoSNOW(clusterobj)

node.metrics.pearson <- foreach(d=1:length(densities_desired), .packages="broom") %dopar% {
  results <- data.frame()
  results.df <- data.frame()
  
  for (m in metrics) {
     
      for (n in nodes) {
        thismetric <- do.call(c, lapply(allmetrics.pearson, function(subj) {
          subj[[d]][[m]][n] #NB!!! If you inadvertently pass in 0.1, 0.2 etc., R rounds up to 1. Therefore, use d (index), not densities_desired[d].
        }))
        
        names(thismetric) <- names(allmetrics.pearson)
        
        #combine these values with identifying variables -- age and BPD
        thisdf <- data.frame(group=factor(attr(allmetrics.pearson, "bpd"), levels=c(0,1), labels=c("control", "bpd")), 
                             age=attributes(allmetrics.pearson)$age, metric=thismetric)
        
        node.test <- tryCatch(t.test(metric~group, thisdf, na.rm = TRUE), error = function(errorname) { print(errorname); return(NULL) })
        browser()
        if (is.null(node.test)) { message("Error occurred for density: ", densities_desired[d], ", metric: ", m, " and node: ", n)} 
        
        if(!is.null(node.test)){
          
          results <- data.frame(subj = SPECC_rest$SPECC_ID, nodenum=n, BPD = SPECC_rest[,"BPD"], Age = SPECC_rest[,"AgeAtScan"], metric.name = m, 
                                density=densities_desired[d], mean.control = node.test$estimate[1], mean.bpd = node.test$estimate[2], group.t=node.test$statistic) 
          
          single.subs <- t(data.frame(thismetric))
          raw <- as.numeric(single.subs)
          browser()
  
          
          z_single.subs <- as.numeric(scale(raw))
    
          results$metric.val <- raw
          results$metric.z.val <- z_single.subs
    
          row.names(results) <- NULL
          results.df <- rbind(results.df, results)
        }
      }
    }
   return(results.df)
  }
save(node.metrics.pearson, file = paste0(basedir,"/cachedRfiles/node.metrics.pearson.corr.RData"))
}

##merge large list into dataframe
node.metrics.pearson.df <- do.call(rbind, node.metrics.pearson)
node.metrics.pearson.df$density <- factor(node.metrics.pearson.df$density)
# levels(node.metrics.pearson.df$density)
densities.desired <- factor(densities.desired)

#### heatmaps across densities
if(!file.exists("den.cors.pdf")){
pdf("den.cors.pdf", width = 18, height = 18)
for(m in metrics){
metric <- subset(node.metrics.pearson.df, metric.name == m)

#####raw metric. ignore z-scores
metric.agg <- aggregate(x = metric$metric.val, by = list(metric$nodenum, metric$density), mean)
colnames(metric.agg) <- c("nodenum", "density", "metric.val")

densities_desired <- as.character(densities_desired)
cormat.den <- matrix(NA, nrow = length(densities_desired), ncol = length(densities_desired))
dimnames(cormat.den) <- list(d1=densities_desired, d2=densities_desired)
for(d1 in densities_desired){
  for(d2 in densities_desired){
    # if(d1==.07){browser()}
    cormat.den[d1,d2] <- cor(metric.agg[which(metric.agg$density ==d1),]$metric.val, metric.agg[which(metric.agg$density == d2),]$metric.val)
    }
  }

melt.cormat.den <- melt(cormat.den)
melt.cormat.den$d1 <- ordered(melt.cormat.den$d1, levels = c(seq(.01,.20,.01))) 
melt.cormat.den$d2 <- ordered(melt.cormat.den$d2, levels = rev(c(seq(.01,.20,.01))))
colnames(melt.cormat.den) <- c("Density_1", "Density_2", "metric.score")
g <- ggplot(data = melt.cormat.den, aes(x = Density_1, y = Density_2, fill = metric.score)) + geom_tile() + ggtitle(paste0(m, "_metric"))
print(g)
  }
dev.off()
}
######################################################################
#PCA on latent degree/ eigenvector across densities

###all metrics (for some reason)
metrics.raw <- node.table.pearson[,c(metrics, "subj", "node", "density")]
#metrics.raw.noNA <- metrics.raw[,c(1:7)]
metrics.raw$density_num <- as.numeric(sub("d(.*)", "\\1", metrics.raw$density, perl=TRUE))
metrics.raw <- metrics.raw %>% select(subj, node, density, density_num, degree, eigen.cent, betweenness.node) %>%
  filter(density_num > .03) %>% select(-density_num) %>% gather(key="variable", value="value", degree, eigen.cent, betweenness.node)

library(data.table)
#dfm <- melt(metrics.raw, id.vars = c("subj", "node", "density"))
  
metrics.raw_pca <- dcast(metrics.raw, subj + node ~ density + variable, value.var = "value")

metrics.pca <- prcomp(select(metrics.raw_pca, -subj, -node), center = TRUE, scale. = TRUE)

library(psych)
pcaout1 <- pca(select(metrics.raw_pca, -subj, -node), nfactors=1)
pcaout2 <- pca(select(metrics.raw_pca, -subj, -node), nfactors=2)
pcaout3 <- pca(select(metrics.raw_pca, -subj, -node), nfactors=3)
pcaout4 <- pca(select(metrics.raw_pca, -subj, -node), nfactors=4)

str(pcaout2$scores)

print(pcaout4$loadings, cutoff = 0.3)

pcasolution <- data.frame(metrics.pca$x[,1:2])
names(pcasolution) <- c("central", "between")

toanalyze <- cbind(select(metrics.raw_pca, subj, node), pcasolution)

####PCA outputs
print(metrics.pca)
plot(metrics.pca, type = "l")
summary(metrics.pca)


library(devtools)
library(ggbiplot)

pdf("PCA_all_metrics.pdf")
g <- ggbiplot(metrics.pca, obs.scale = 1, var.scale = 1,  ellipse = TRUE, circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)
dev.off()

####degree PCA
degree.raw <- metrics.raw.noNA[,1]
degree.pca <- prcomp(degree.raw, center = TRUE, scale. = TRUE)
print(degree.pca)


##SCOTMI node table
node.table.scotmi <- data.frame()

names(allmetrics.scotmi) <- SPECC_rest$SPECC_ID #assumes elements of allmetrics in same order as this DF
node.table.scotmi <- do.call(rbind, lapply(1:length(allmetrics.scotmi), function (s) {
  subj <- allmetrics.scotmi[[s]]
  
  ddf <- lapply(1:length(subj), function(d) {
    density <- subj[[d]]
    
    metrics_to_grab <- names(density)[!names(density) %in% c("origgraph")]
    nvec <- lapply(1:length(metrics_to_grab), function(m) { 
      density[[metrics_to_grab[m]]]
    })
    
    names(nvec) <- metrics_to_grab
    ndf <- as.data.frame(do.call(cbind, nvec))
    ndf$density <- names(subj)[d]
    ndf$node <- atlas$vname
    ndf$anat_label <- atlas$anat_label
    ndf$scotmi <- 1
    return(ndf)
  })
  
  #names(ddf) <- names(subj)
  
  ddf <- do.call(rbind, ddf)
  ddf$subj <- names(allmetrics.scotmi)[s]
  ddf$bpd <- attr(allmetrics.scotmi, "bpd")[s]
  ddf$age <- attr(allmetrics.scotmi, "age")[s]
  
  return(ddf)
  
}))

###################
##TA DA!!!! now table with nodes, subjects, networks, and scotmi vs pearson correlations
node.table <- rbind(node.table.pearson, node.table.scotmi)
tail(node.table)

###Now: correlation matrices
library(data.table)
dcast.rs <- dcast(node.table, subj + node ~ density + scotmi, value.var = "degree")
corr.scot <- round(cor(dcast.rs[,c(-1, -2)]), 3)
melt.corr.scot <- melt(corr.scot)
melt.corr.scot$Var1 <- ordered(melt.corr.scot$Var1, levels=c(paste0("d", seq(0.01, 0.2, .01), "_0"), paste0("d", seq(0.01, 0.2, .01), "_1")))
melt.corr.scot$Var2 <- ordered(melt.corr.scot$Var2, levels=rev(c(paste0("d", seq(0.01, 0.2, .01), "_0"), paste0("d", seq(0.01, 0.2, .01), "_1"))))

#melt.corr.scot$sort <- grep('*_0', melt.corr.scot[,1])
#melt.sort0 <- melt.corr.scot[sort(grep('*_0', melt.corr.scot[,1])), ]
#melt.sort1 <- melt.corr.scot[sort(grep('*_1', melt.corr.scot[,1])), ]
#melt.sort <- rbind(melt.sort0, melt.sort1)


# pdf("scotmi_pearson_degree_allDensities.pdf", width = 15, height = 15)
# ggplot(data = melt.corr.scot, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + geom_text(aes(label=value)) 
# dev.off()

###edge comparisons on weighted matrices across subjects for scotmi and pearson

if(file.exists(paste0(basedir, "/cachedRfiles/esimmat_pearscot.RData")) == TRUE){
  esimmat  <- get(load(paste0(basedir, "/cachedRfiles/esimmat_pearscot.RData")))
} else{
nsubj <- length(allg.weighted.pearson)
esimmat <- matrix(NA, nrow=nsubj*2, ncol=nsubj*2)

for (s1 in 1:nsubj) {
  for (s2 in 1:nsubj) {
    
    ##pearson with pearson
    esimmat[s1,s2] <- cor(E(allg.weighted.pearson[[s1]])$weight, E(allg.weighted.pearson[[s2]])$weight)
    
    ##scotmi with scotmi
    esimmat[s1+nsubj,s2+nsubj] <- cor(E(allg.weighted.scotmi[[s1]])$weight, E(allg.weighted.scotmi[[s2]])$weight)
    
    esimmat[s1,s2+nsubj] <- cor(E(allg.weighted.pearson[[s1]])$weight, E(allg.weighted.scotmi[[s2]])$weight)
    esimmat[s1+nsubj,s2] <- cor(E(allg.weighted.scotmi[[s1]])$weight, E(allg.weighted.pearson[[s2]])$weight)
  }
  
}
save(esimmat, file = paste0(basedir,"/cachedRfiles/esimmat_pearscot.RData"))
}

melt.esimmat <- melt(esimmat)

melt.esimmat$Var1 <- ordered(melt.esimmat$Var1, levels = c(seq(1,158,1)))
melt.esimmat$Var2 <- ordered(melt.esimmat$Var2, levels = rev(c(seq(1,158,1))))

pdf("between_subj_edge_similarity.pdf", width = 18, height = 18)
ggplot(data = melt.esimmat, aes(x = Var1, y = Var2, fill = value)) + geom_tile()
dev.off()

#########
##check out degree distributions of scotmi and pearson subjects across 3 densities (5,10,15)
pdf("degree_dist_scot.pdf", width = 8.5, height = 11)
subs <- sample(1:79, 5)
for (i in subs){
  deg.05.p <- data.frame(density = 5, scotmi = 0, metric = allmetrics.pearson[[i]][[5]][["degree"]]) 
  deg.1.p <- data.frame(density = 10, scotmi = 0, metric = allmetrics.pearson[[i]][[10]][["degree"]]) 
  deg.15.p <- data.frame(density = 15, scotmi = 0, metric = allmetrics.pearson[[i]][[15]][["degree"]]) 
  deg.05.s <- data.frame(density = 5, scotmi = 1, metric = allmetrics.scotmi[[i]][[5]][["degree"]]) 
  deg.1.s <- data.frame(density = 10, scotmi = 1, metric = allmetrics.scotmi[[i]][[10]][["degree"]]) 
  deg.15.s <- data.frame(density= 15, scotmi = 1, metric = allmetrics.scotmi[[i]][[15]][["degree"]]) 
  
  degs.table <- rbind(deg.05.s, deg.1.s, deg.15.s, deg.05.p, deg.1.p, deg.15.p)
  g <- ggplot(degs.table, aes(x = metric)) + geom_histogram() + facet_grid(density ~ scotmi) + ggtitle((paste0("Subject ", i, " Degree Distribution"))) + xlab("degree (by SCOTMI, 1 =yes)") + ylab("Density")
  print(g)
}
dev.off()

pdf("degree_dist_scot_logTransform.pdf", width = 8.5, height = 11)
for (i in subs){
  deg.05.p <- data.frame(density = 5, scotmi = 0, metric = allmetrics.pearson[[i]][[5]][["degree"]]) 
  deg.1.p <- data.frame(density = 10, scotmi = 0, metric = allmetrics.pearson[[i]][[10]][["degree"]]) 
  deg.15.p <- data.frame(density = 15, scotmi = 0, metric = allmetrics.pearson[[i]][[15]][["degree"]]) 
  deg.05.s <- data.frame(density = 5, scotmi = 1, metric = log(allmetrics.scotmi[[i]][[5]][["degree"]])) 
  deg.1.s <- data.frame(density = 10, scotmi = 1, metric = log(allmetrics.scotmi[[i]][[10]][["degree"]])) 
  deg.15.s <- data.frame(density= 15, scotmi = 1, metric = log(allmetrics.scotmi[[i]][[15]][["degree"]])) 
  
  degs.table <- rbind(deg.05.s, deg.1.s, deg.15.s, deg.05.p, deg.1.p, deg.15.p)
  g <- ggplot(degs.table, aes(x = metric)) + geom_histogram() + facet_grid(density ~ scotmi) + ggtitle((paste0("Subject ", i, " Degree Distribution"))) + xlab("degree (by SCOTMI, 1 =yes)") + ylab("Density")
  print(g)
}
dev.off()


###Z-scored and centered metric
# metric.agg.z <- aggregate(x = metric$metric.z.val, by = list(metric$nodenum, metric$density), mean)
# colnames(metric.agg.z) <- c("nodenum", "density", "metric.z.val")
# 
# cormat.den.z <- matrix(NA, nrow = length(densities_desired), ncol = length(densities_desired))
# for(d1 in densities_desired){
#   for(d2 in densities_desired){
#     
#     cormat.den.z[d1*100,d2*100] <- cor(metric.agg.z[which(metric.agg.z$density ==d1),]$metric.z.val, metric.agg.z[which(metric.agg.z$density == d2),]$metric.z.val)
#   }
# }
# 
# cormat.den.z[7,] <- 0
# cormat.den.z[,7] <- 0
# 
# melt.cormat.den.z <- melt(cormat.den.z)
# melt.cormat.den.z$Var1 <- ordered(melt.cormat.den.z$Var1, levels = c(seq(1,20,1))) 
# melt.cormat.den.z$Var2 <- ordered(melt.cormat.den.z$Var2, levels = rev(c(seq(1,20,1))))
# colnames(melt.cormat.den.z) <- c("Density_1", "Density_2", "z.score")
# g.z <- ggplot(data = melt.cormat.den.z, aes(x = Density_1, y = Density_2, fill = z.score)) + geom_tile() + ggtitle(paste0(m, "_Z"))
# print(g.z)

 