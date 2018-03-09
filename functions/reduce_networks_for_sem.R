reduce_networks_forsem <- function(toanalyze_binary, network, prop.thresh){
  #prop.thresh denotes the proportion of nodes that will be retained before factor loading thresholding
  within.mod <- toanalyze_binary %>% dplyr::filter(membership == network) %>% dplyr::select(node, id, within.mod)
  
  tocor <- within.mod %>% spread(node, within.mod) 
  corr <- cor(tocor[,-1])
  diag(corr)=0
  
  rank <- data.frame(mean_corr = apply(corr, 1, function(x) {mean(x)}))
  rank$node <- rownames(rank)
  attach(rank)
  rank <- rank[order(-mean_corr),]
  row.names(rank) <- NULL
  detach(rank)
  
  within.modmean <- within.mod %>% spread(id, within.mod) 
  within.modmean <- within.modmean %>% mutate(mean_within.mod = rowMeans(dplyr::select(within.modmean, -node))) %>% dplyr::select(node, mean_within.mod)
  attach(within.modmean)
  within.modmean <- within.modmean[order(-mean_within.mod),]
  detach(within.modmean)
  
  rank_corr <- left_join(within.modmean, rank, by = "node")
  x <- cor.test(rank_corr[,2], rank_corr[,3], method = "spearman")
  
  net.names<- unique(community.names$community.name)
  message("pval for spearman correlation (mean corr and within.mod) for ",net.names[network], ": ", x$p.value)
  
  # splom(rank_corr[,-1])
  
  #single FA
  number_retain <- round(length(rank$node)*prop.thresh, 0) #keep a specified proportion of nodes in a network rather than enforcing 10 or 20 across networks
  
  retain_nodes <- as.character(head(within.modmean$node, number_retain))
  
  fa.obj <- dcast(within.mod, id~node, value.var = "within.mod")
  fa_reduced <- fa.obj[,retain_nodes]
  fa_reduced$id <- fa.obj$id
  
  faout <- fa(dplyr::select(fa_reduced, -id), nfactors = 1, fm = "ml")
  faout$loadings
  
  ###threshold nodes with factor loadings below .3
  fa.load.thresh <- names(faout$loadings[which(faout$loadings > .3),])
  fa_reduced.load.thresh <- fa.obj[,fa.load.thresh]
  fa_reduced.load.thresh$id <- fa.obj$id
  return(fa_reduced.load.thresh)
}