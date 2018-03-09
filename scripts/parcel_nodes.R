
# Parcel nodes script for input to SEM------------------------------------------------------------
setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
library(NbClust)
####read in package dependencies and custom functions
source("functions/setup_globals.R")
toanalyze_binary <- reduce_centrality_fa(allmetrics.nodal.df, reducemetrics = reducemetrics, allowCache = TRUE, weighted = FALSE, browse = FALSE)
toanalyze_binary <- left_join(toanalyze_binary, membership_df, by = "node")
toanalyze_binary$id <- as.character(toanalyze_binary$id)

str(toanalyze_binary)


dmn_central <- toanalyze_binary %>% dplyr::filter(membership == 7) %>% dplyr::select(node,id,central)
str(dmn_central)


# DMN ---------------------------------------------------------------------
sub_partitions <- list()
for (sub in 1:length(unique(toanalyze_binary$id))){
  subnum <- unique(toanalyze_binary$id)[sub]
  submetrics <- dplyr::filter(toanalyze_binary, id == subnum)
  
  network_partitions <- list()
  for (net in 1:length(unique(toanalyze_binary$membership))){
    net_centrals <- submetrics %>% dplyr::filter(membership == net) %>% dplyr::select(central, integration, within.mod, betweenness.node)
    net_centrals <- scale(net_centrals)
    
    set.seed(20)  
    clust <- NbClust(net_centrals, min.nc = 5, max.nc = 25, method = "kmeans")
    network_partitions[[net]] <- clust
    rm(net_centrals)
  }
  sub_partitions[[sub]] <- network_partitions
}
