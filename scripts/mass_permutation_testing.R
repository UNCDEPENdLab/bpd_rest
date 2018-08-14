#mass permutation test

#setwd("~/Data_Analysis/bpd_rest")
#setwd("/gpfs/group/mnh5174/default/Michael/bpd_rest")
library(tidyverse)
#library(doSNOW)
library(doParallel)
library(data.table)
library(igraph)
s <- makePSOCKcluster(10)
registerDoParallel(s)
library(coin)
#helper function to extract test statistic and pvalue
coinstats <- function(m) {
  tryCatch(data.frame(mdiff=as.numeric(statistic(m, type="test")), pvalue=as.numeric(pvalue(m))), error=function(e) {
    return(data.frame(mdiff=NA_real_, pvalue=NA_real_))
    })
}


#load("/Users/mnh5174/Data_Analysis/bpd_rest/cache/threshnodalmetrics_schaefer422_nosmooth_aroma_bp_nonaggr_cor.shrink_fc_binary_all.RData")
#load("cache/threshnodalmetrics_schaefer422_nosmooth_aroma_bp_nonaggr_cor.shrink_fc_binary_all.RData")
#load("cache/threshnodalmetrics_schaefer421_nosmooth_aroma_hp_cor.shrink_fc_binary_all.RData") #HP + cor.shrink + 421
#load("cache/threshnodalmetrics_schaefer421_nosmooth_aroma_bp_nonaggr_cor.shrink_fc_binary_all.RData") #BP + cor.shrink + 421

pipeline <- "braingraph_consensus_all_hp421_binary" #hp421
#pipeline <- "braingraph_consensus_all_hp421_weighted" #hp421
#pipeline <- "braingraph_consistency_hp421_binary" #hp421
#pipeline <- "braingraph_consistency_hp421_weighted" #hp421

#look at braingraph consistency basis
bg.g <- readRDS("braingraph_temp/consensus_all/consensus_all_g_binary.rds") #braingraph data for binary using the consensus all
#bg.g <- readRDS("braingraph_temp/consensus_all/consensus_all_g_weighted.rds") #braingraph data for weighted using the consensus all
#bg.g <- readRDS("braingraph_temp/consistency/consistency_g_weighted.rds") #braingraph data for weighted consistency
#bg.g <- readRDS("braingraph_temp/consistency/consistency_g_binary.rds") #braingraph data for weighted consistency

#flatten the bg stuff
allmetrics.nodal.df <- rbindlist(lapply(bg.g, function(grouplist) {
  rbindlist(lapply(grouplist, function(gthresh) {
    #gthresh has one graph per subject
    df <- lapply(gthresh, function(g) {
      data.frame(id=g$name, wthresh=g$threshold, node=V(g)$name, eigen.cent=V(g)$ev.cent, degree=V(g)$degree,
        betweenness.node=V(g)$btwn.cent, part.coeff=V(g)$PC, gateway.coeff.btw=V(g)$GC, within.module.deg.zscore=V(g)$z.score,
        local.clustering=V(g)$transitivity, leverage.cent=V(g)$lev.cent)
    })
    rbindlist(df)
  }))
}))
  


allmetrics.nodal.df$wthresh_char <- as.character(round(allmetrics.nodal.df$wthresh, 3))

#subj_info <- gdata::read.xls("/Users/mnh5174/Box Sync/DEPENd/Projects/SPECC/ID Management/SPECC_Participant_Info.xlsx")
subj_info <- read.csv("data/SPECC_Participant_Info.csv")

allmetrics.nodal.df <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID) %>% inner_join(allmetrics.nodal.df) %>%
  mutate(BPD=factor(BPD, levels=c("0", "1"), labels=c("Control", "BPD")))
str(allmetrics.nodal.df)
## lookup <- read.csv("data/schaefer422_masterlookup.csv") %>%
##   dplyr::mutate(vname=paste0("V", vname)) %>% dplyr::rename(node=vname) %>%
##   mutate(tag=paste(node, anat_label, X7_networks, sep=" - "))

lookup <- read.csv("data/schaefer421_masterlookup.csv") %>%
  dplyr::mutate(vname=paste0("roi", vname)) %>% dplyr::rename(node=vname) %>%
  mutate(tag=paste(node, anat_label, X7_networks, sep=" - "))

allmetrics.nodal.df <- allmetrics.nodal.df %>% inner_join(lookup)

## metrics <- c("local.clustering", "degree", "eigen.cent", "closeness", "betweenness.node", "page.rank",
##              "leverage.cent", "within.module.deg.zscore", "between.module.deg.zscore",
##              "within.module.deg", "between.module.deg", "part.coeff", "gateway.coeff.btw", "gateway.coeff.degree")

metrics <- c("local.clustering", "degree", "eigen.cent", "betweenness.node",
             "leverage.cent", "within.module.deg.zscore", "part.coeff", "gateway.coeff.btw")

#for speeding up testing
#allmetrics.nodal.df <- allmetrics.nodal.df %>% filter(node %in% paste0("V", 1:5))

#look at missingness for local clustering
#mdf %>% group_by(node, wthresh_char, BPD) %>% summarize(nonmiss=sum(!is.na(metric))) %>% filter(nonmiss < 2)

allmetrics_stats <- foreach(i=1:length(metrics), .packages=c("magrittr", "coin", "dplyr", "ggplot2", "tidyr", "purrr")) %dopar% {
  #for (i in 1:length(metrics)) {
  #first deal with problem of high missingness in some comparisons
  #(at the extreme, zero observations in a given group, which will crash the permutation test)
  
  #note that dplyr::rename secretly supports dynamic arguments: https://github.com/tidyverse/dplyr/issues/1600
  mdf <- allmetrics.nodal.df %>% dplyr::rename(metric := !!metrics[i]) %>% group_by(node, wthresh_char, BPD) %>% filter(!is.na(metric)) %>% filter(n() >= 20) %>% #require at least 20 per group at each node and threshold to run the permutation test
    ungroup() %>% group_by(node, wthresh_char) %>% filter(n() > 40) %>% ungroup() %>% droplevels() #now drop any wthresh + node combination with fewer that 40 obs (would just occur if one group is at zero

  #still get missingness problems if one group has 0 obs and the other group has > 20. This is fixed above by the double group_by approach
  #test <- mdf %>% group_by(node, wthresh_char) %>% summarize(nobs=n()) #select(metric, BPD, node, wthresh_char) %>% do({as.data.frame(table(.$BPD)) }) %>% ungroup()
  #range(test$Freq)
  
  stats <- mdf %>% group_by(node, wthresh_char) %>% select(metric, BPD, node, wthresh_char) %>% nest() %>%
    mutate(model=map(data, ~ independence_test(metric ~ BPD, data=., distribution=approximate(B=40000)))) %>% #use approximate null with 40k MC replications
    #mutate(model=map(data, ~ independence_test(metric ~ BPD, data=., distribution="exact"))) #the exact null requires an insane amount of RAM (like 40GB+ for one model)
    unnest(model %>% purrr::map(coinstats)) %>% select(-model, -data) %>% inner_join(lookup)
  
  pvalballpark <- stats %>% group_by(node) %>% dplyr::summarize(anyp = median(pvalue < .10)) %>%
    filter(anyp==TRUE) %>% pull(node)
  
  pdf(paste0("figures/", metrics[i], "_", pipeline, "_stats.pdf"), width=8, height=6)
  dsplit <- filter(stats, node %in% pvalballpark) %>% split(.$node)
  for(d in dsplit) {
    g <- ggplot(d, aes(x=wthresh_char,y=pvalue)) + geom_point() + ggtitle(d$tag[1])
    plot(g)
  }
  dev.off()
  
  #allstats[[metrics[i]]] <- stats
  stats
}

stopCluster(s)

names(allmetrics_stats) <- metrics
save(allmetrics_stats, file=paste0("allperms_", pipeline, ".RData"))
