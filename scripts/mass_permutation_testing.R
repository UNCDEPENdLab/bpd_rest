#mass permutation test
#setwd("~/Data_Analysis/bpd_rest")
setwd("/gpfs/group/mnh5174/default/Michael/bpd_rest")
library(tidyverse)
library(doSNOW)
s <- makeSOCKcluster(8)
registerDoSNOW(s)
library(coin)
#helper function to extract test statistic and pvalue
coinstats <- function(m) { data.frame(mdiff=as.numeric(statistic(m, type="test")), pvalue=as.numeric(pvalue(m))) }

#load("/Users/mnh5174/Data_Analysis/bpd_rest/cache/threshnodalmetrics_schaefer422_nosmooth_aroma_bp_nonaggr_cor.shrink_fc_binary_all.RData")
#load("cache/threshnodalmetrics_schaefer422_nosmooth_aroma_bp_nonaggr_cor.shrink_fc_binary_all.RData")
load("cache/threshnodalmetrics_schaefer421_nosmooth_aroma_hp_cor.shrink_fc_binary_all.RData") #HP + cor.shrink + 421
allmetrics.nodal.df$wthresh_char <- as.character(round(allmetrics.nodal.df$wthresh, 3))

#subj_info <- gdata::read.xls("/Users/mnh5174/Box Sync/DEPENd/Projects/SPECC/ID Management/SPECC_Participant_Info.xlsx")
subj_info <- read.csv("data/SPECC_Participant_Info.csv")

allmetrics.nodal.df <- subj_info %>% select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID) %>% inner_join(allmetrics.nodal.df) %>%
  mutate(BPD=factor(BPD, levels=c("0", "1"), labels=c("Control", "BPD")))

## lookup <- read.csv("data/schaefer422_masterlookup.csv") %>%
##   dplyr::mutate(vname=paste0("V", vname)) %>% dplyr::rename(node=vname) %>%
##   mutate(tag=paste(node, anat_label, X7_networks, sep=" - "))

lookup <- read.csv("data/schaefer421_masterlookup.csv") %>%
  dplyr::mutate(vname=paste0("roi", vname)) %>% dplyr::rename(node=vname) %>%
  mutate(tag=paste(node, anat_label, X7_networks, sep=" - "))

allmetrics.nodal.df <- allmetrics.nodal.df %>% inner_join(lookup)

metrics <- c("local.clustering", "degree", "eigen.cent", "closeness", "betweenness.node", "page.rank",
             "leverage.cent", "within.module.deg.zscore", "between.module.deg.zscore",
             "within.module.deg", "between.module.deg", "part.coeff", "gateway.coeff.btw", "gateway.coeff.degree")

#for speeding up testing
#allmetrics.nodal.df <- allmetrics.nodal.df %>% filter(node %in% paste0("V", 1:5))

allmetrics_stats <- foreach(i=1:length(metrics), .packages=c("magrittr", "coin", "dplyr", "ggplot2", "tidyr", "purrr")) %dopar% {
  #for (i in 1:length(metrics)) {
  mdf <- allmetrics.nodal.df
  names(mdf)[which(names(mdf) == metrics[i])] <- "metric"
  stats <- mdf %>% group_by(node, wthresh_char) %>% select(metric, BPD, node, wthresh_char) %>% nest() %>%
    mutate(model=map(data, ~ independence_test(metric ~ BPD, data=., distribution=approximate(B=40000)))) %>% #use approximate null with 40k MC replications
    #mutate(model=map(data, ~ independence_test(metric ~ BPD, data=., distribution="exact"))) #the exact null requires an insane amount of RAM (like 40GB+ for one model)
    unnest(model %>% purrr::map(coinstats)) %>% select(-model, -data) %>% inner_join(lookup)
  
  pvalballpark <- stats %>% group_by(node) %>% dplyr::summarize(anyp = median(pvalue < .1)) %>%
    filter(anyp==TRUE) %>% pull(node)
  
  pdf(paste0("figures/", metrics[i], "_hp421_stats.pdf"), width=8, height=6)
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
save(allmetrics_stats, file="allperms.RData")
