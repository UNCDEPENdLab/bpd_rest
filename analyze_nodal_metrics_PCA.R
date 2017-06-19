#Analyses of nodal statistics using PCA to reduce across densities and metrics
#setup package dependencies and custom functions
#setwd("~/Box Sync/RS_BPD_graph")
# setwd("/Users/mnh5174/Data_Analysis/bpd_rest")
setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/")
basedir <- getwd()

require(psych)

source("functions/setup_globals.R") #this will setup details of the parcellation, conn_method, preproc_pipeline, and connection distance
source("functions/graph_util_redux.R") 
source("functions/get_subj_info.R")

#load nodal statistics for the current pipeline
nodalmetrics_dthresh_df <- load_nodal_metrics_df() #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame

#not used at the moment
#nodalmetrics_dthresh_df$density_fac <- factor(nodalmetrics_dthresh_df$density) #to overcome insanity of floating point imprecision if we use density==X approach 

metrics_to_analyze <- c("degree", "eigen.cent", "betweenness.node", "within.module.deg.zscore", "between.module.deg.zscore")

######################################################################
#PCA on latent degree/ eigenvector across densities

#stack metrics in data.frame to have subjects, nodes, densities, *and* metrics on the rows
metrics.raw <- nodalmetrics_dthresh_df %>% select_(.dots = c(metrics_to_analyze, "id", "node", "density")) %>%
    filter(density > .04) %>% gather(key="variable", value="value", degree, eigen.cent, betweenness.node, within.module.deg.zscore, between.module.deg.zscore)


#recast for PCA such that column names represent a metrics at a given density (e.g., 0.05_between.module.deg.zscore)
metrics.raw_pca <- dcast(metrics.raw, id + node ~ density + variable, value.var = "value")

metrics.pca <- prcomp(select(metrics.raw_pca, -id, -node), center = TRUE, scale. = TRUE)

pcaout1 <- pca(select(metrics.raw_pca, -id, -node), nfactors=1, rotate="varimax") #varimax (orthogonal) is default, just making that clear here
pcaout2 <- pca(select(metrics.raw_pca, -id, -node), nfactors=2, rotate="varimax")
pcaout3 <- pca(select(metrics.raw_pca, -id, -node), nfactors=3, rotate="varimax")
pcaout4 <- pca(select(metrics.raw_pca, -id, -node), nfactors=4, rotate="varimax")
pcaout5 <- pca(select(metrics.raw_pca, -id, -node), nfactors=5, rotate="varimax")
pcaout6 <- pca(select(metrics.raw_pca, -id, -node), nfactors=6, rotate="varimax")

print(pcaout4$loadings, cutoff = 0.5)
print(pcaout3$loadings, cutoff = 0.5)

pcasolution <- data.frame(pcaout3$scores) #use these, which are the rotated scores, not the raw eigenvectors from $x in prcomp
names(pcasolution) <- c("central", "within.mod", "between.node")

#this was from Nate's analysis of Pearson. For the AROMA ridge, 3 factors covers 90% of variaence
#pcasolution <- data.frame(metrics.pca$x[,1:4])
#names(pcasolution) <- c("central", "between.node", "within.mod", "between.mod")
#

toanalyze <- cbind(select(metrics.raw_pca, id, node), pcasolution)
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz")
merge.bpdage <- subj_info %>% select(SPECC_ID, BPD, AgeAtScan) %>% rename(id=SPECC_ID, Age=AgeAtScan)

toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")


##################################################################################################
##############
####run tests on nodal metrics collapsed across densities
abst <- 2.64 ##about p < .005

metrics.toanalyze <- names(select(toanalyze, -id, -node, -BPD, -Age))
sigres <- 1
results <- list()

for (m in metrics.toanalyze) {
  for (n in levels(toanalyze$node)) {
    thismetric <- toanalyze[toanalyze$node == n, c(m, "BPD", "Age")]
    thismetric$BPD <- factor(thismetric$BPD, levels = c(0,1), labels = c("control", "BPD"))
    colnames(thismetric) <- c("metric", "group", "age")
    
    node.test <- tryCatch(t.test(metric~group, thismetric, na.rm = TRUE), error = function(errorname) { print(errorname); return(NULL) })
    if (is.null(node.test)) { message("Error occurred for metric: ", m, " and node: ", n) }
    
    age.test <- tryCatch(lm(metric~age*group, thismetric, na.action = "na.exclude"), error=function(e) { print(e); return(NULL) })
    if (is.null(age.test)) { pvec <- NULL
    } else { pvec <- broom::tidy(age.test)$p.value[-1] } #p-values of age, bpd, and age x bpd. -1 to drop off intercept
    
    if ((!is.null(node.test) && !is.nan(node.test$statistic) && abs(node.test$statistic) > abst) ||
        (!is.null(age.test) && !all(is.nan(pvec)) && any(pvec < .01))) {
      
      results[[sigres]] <- list(nodename=as.character(atlas$anat_label[atlas$name == n]), ttest=node.test, 
          metric = m, agetest=age.test, nodenum=n)
      sigres <- sigres+1
    }
  }
}



results.df <- do.call(rbind, results)

####BPD main effects 
results.ttest <- lapply(results, function(node){
      df <- broom::tidy(node$ttest)
      df$nodename <- node$nodename
      df$nodenum <- node$nodenum
      df$metric <- node$metric
      return(df)
    })

results.ttest.df <- do.call(rbind, results.ttest)
a <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "central" & p.value < .005)
b <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "between.node" & p.value < .005)
c <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "within.mod" & p.value < .005)
d <- results.ttest.df %>% arrange(nodename, metric) %>% filter(metric == "between.mod" & p.value < .005)

####Age and agexbpd lm 
results.agetest <- lapply(results, function(node){
      df <- broom::tidy(node$agetest)
      df$nodename <- node$nodename
      df$nodenum <- node$nodenum
      df$metric <- node$metric
      return(df)
    })
results.agetest.df <- do.call(rbind, results.agetest)
###age main effects
e <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "central" & term == "age" & p.value < .01)
f <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "between.node" & term == "age" & p.value < .01)
g <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "within.mod" & term == "age" & p.value < .01)
h <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "between.mod" & term == "age" & p.value < .01)

###age x bpd interactions
i <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "central" & term == "age:groupBPD" & p.value < .01)
j <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "between.node" & term == "age:groupBPD" & p.value < .01)
k <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "within.mod" & term == "age:groupBPD" & p.value < .01)
l <- results.agetest.df %>% arrange(nodename, metric) %>% filter(metric == "between.mod" & term == "age:groupBPD" & p.value < .01)

library(gtools)
all.sig.nodal.pca <- bind_rows(a,b,c,d,e,f,g,h,i,j,k,l)
write.csv(all.sig.nodal.pca, file = paste0(basedir, "/output.files/all.sig.nodal.pca.",pipeline,".csv"))

#### OLD LEFTOVERS

#this is obviated by having the nodal statistics digested to data.frame inside the setup function: compute_nodal_metrics
#nodalmetrics_dthresh_df <- data.frame()
#
#names(allmetrics) <- SPECC_rest$SPECC_ID #assumes elements of allmetrics in same order as this DF
#attr(allmetrics, "bpd") <- SPECC_rest[,"BPD"]
#attr(allmetrics, "age") <- SPECC_rest[,"AgeAtScan"]
#
#nodalmetrics_dthresh_df <- do.call(rbind, lapply(1:length(allmetrics), function (s) {
#          subj <- allmetrics[[s]]
#          
#          ddf <- lapply(1:length(subj), function(d) {
#                density <- subj[[d]]
#                
#                metrics_to_grab <- names(density)[!names(density) %in% c("origgraph")]
#                nvec <- lapply(1:length(metrics_to_grab), function(m) { 
#                      density[[metrics_to_grab[m]]]
#                    })
#                
#                names(nvec) <- metrics_to_grab
#                ndf <- as.data.frame(do.call(cbind, nvec))
#                ndf$density <- names(subj)[d]
#                ndf$node <- atlas$vname
#                ndf$anat_label <- atlas$anat_label
#                ndf$scotmi <- 0
#                return(ndf)
#              })
#          
#          #names(ddf) <- names(subj)
#          
#          
#          ddf <- do.call(rbind, ddf)
#          ddf$subj <- names(allmetrics)[s]
#          ddf$bpd <- attr(allmetrics, "bpd")[s]
#          ddf$age <- attr(allmetrics, "age")[s]
#          
#          return(ddf)
#          
#        }))
