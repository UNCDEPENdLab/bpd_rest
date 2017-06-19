#Analyses of nodal statistics using PCA to reduce across densities and metrics
#setup package dependencies and custom functions
#setwd("~/Box Sync/RS_BPD_graph")
# # setwd("/Users/mnh5174/Data_Analysis/bpd_rest")
# setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/")
# basedir <- getwd()

analyze_nodal_metrics_PCA <- function(nodalmetrics_dthresh_df = NULL, allowCache = TRUE){

require(psych)

stopifnot(file.exists(file.path(basedir, "cache")))  
expectFile <- file.path(basedir, "output.files", paste0("all.sig.nodal.pca.comm.", preproc_pipeline, ".", conn_method, ".csv"))

if (file.exists(expectFile) && allowCache==TRUE) {
  message("Loading significant results for PCA analysis from file: ", expectFile)
  all.sig.nodal.pca.comm <- read.csv(expectFile, header = TRUE)
} else {  

#load nodal statistics for the current pipeline
if(is.null(nodalmetrics_dthresh_df)){
  #this returns allmetrics.nodal as nested list and allmetrics.nodal.df as flat data.frame
  nodalmetrics_dthresh_df <- load_nodal_metrics_df()
} else {

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

# pcaout1 <- pca(select(metrics.raw_pca, -id, -node), nfactors=1, rotate="varimax") #varimax (orthogonal) is default, just making that clear here
# pcaout2 <- pca(select(metrics.raw_pca, -id, -node), nfactors=2, rotate="varimax")
pcaout3 <- pca(select(metrics.raw_pca, -id, -node), nfactors=3, rotate="varimax")
pcaout4 <- pca(select(metrics.raw_pca, -id, -node), nfactors=4, rotate="varimax")
# pcaout5 <- pca(select(metrics.raw_pca, -id, -node), nfactors=5, rotate="varimax")
# pcaout6 <- pca(select(metrics.raw_pca, -id, -node), nfactors=6, rotate="varimax")

# print(pcaout4$loadings, cutoff = 0.5)
print(pcaout3$loadings, cutoff = 0.5) #if you need to look over PC loadings

pcasolution <- data.frame(pcaout3$scores) #use these, which are the rotated scores, not the raw eigenvectors from $x in prcomp
names(pcasolution) <- c("central", "within.mod", "between.node")

#this was from Nate's analysis of Pearson. For the AROMA ridge, 3 factors covers 90% of variaence
#pcasolution <- data.frame(metrics.pca$x[,1:4])
#names(pcasolution) <- c("central", "between.node", "within.mod", "between.mod")
#

toanalyze <- cbind(select(metrics.raw_pca, id, node), pcasolution)
#subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz")
merge.bpdage <- subj_info %>% select(SPECC_ID, BPD, AgeAtScan) %>% rename(id=SPECC_ID, Age=AgeAtScan)

toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")
write.csv(toanalyze, file = paste0(basedir, "/cache/toanalyze.pca.",preproc_pipeline, ".", conn_method, ".csv"))

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

community <- readRDS(paste0(getwd(), "/cache/d10_louv_MH.rds"))
membership <- data.frame(community$membership)
membership$nodenum <- row.names(membership)

all.sig.nodal.pca.comm <- dplyr::left_join(all.sig.nodal.pca, membership, by = "nodenum")
all.sig.nodal.pca.comm$community.membership <- mapvalues(all.sig.nodal.pca.comm$community.membership, from = c("1","2","3","4","5"), to = c("OCC", "SOMMOTOR", "FPN", "CO", "DMN"))


write.csv(all.sig.nodal.pca.comm, file = paste0(basedir, "/output.files/all.sig.nodal.pca.comm.",preproc_pipeline, ".", conn_method, ".csv"))

    }
  }
return(all.sig.nodal.pca.comm)
}

