generate_heatmaps <- function(agg.g, graph_list){
library(dplyr)
library(viridis)
agg.adjmat <- as_adjacency_matrix(agg.g, attr = "weight", sparse = FALSE)
#str(agg.adjmat)

membership <- read.csv(file.path(getwd(), "data", "membership.yeo.csv"))
membership <- membership %>% rename(network=x) %>% mutate(node=1:nrow(membership))
reorder <- with(membership,order(network, node))
#membership[reorder,]

agg.adjmat.resort <- agg.adjmat[reorder,reorder]

####sorted heatmap on aggregate data
png(paste0(conn_method, "_heatmap.png"), width=25, height=25, units="in", res=300)
agg.melt <- melt(agg.adjmat, varnames=c("Node1", "Node2"))
agg.melt$Node1 <- ordered(agg.melt$Node1, levels=paste0("V", reorder))
agg.melt$Node2 <- ordered(agg.melt$Node2, levels=paste0("V", reorder))

ggplot(agg.melt, aes(x=Node1, y=Node2, fill=value)) + geom_tile() + scale_fill_viridis()
dev.off()
#### unsorted heatmap on aggregate data
png(paste0(conn_method, "_heatmap_unordered.png"), width=25, height=25, units="in", res=300)
agg.melt <- melt(agg.adjmat, varnames=c("Node1", "Node2"))
ggplot(agg.melt, aes(x=Node1, y=Node2, fill=value)) + geom_tile() + scale_fill_viridis()
dev.off()


##sorted heatmaps of 5 random subjects
rand.subs <- sample(1:length(graph_list), 5, replace = F)
pdf(paste0(conn_method, "_subj_rand_heatmap.pdf"), width=25, height=25)#, units="in", res=300)
for(subj in rand.subs){
  this.subj.adjmat <- as_adjacency_matrix(graph_list[[subj]], attr = "weight", sparse = FALSE)
  this.subj.melt <- melt(this.subj.adjmat, varnames = c("Node1", "Node2"))
  this.subj.melt$Node1 <- ordered(this.subj.melt$Node1, levels=paste0("V", reorder))
  this.subj.melt$Node2 <- ordered(this.subj.melt$Node2, levels=paste0("V", reorder))
  ggplot(this.subj.melt, aes(x=Node1, y=Node2, fill=value)) + geom_tile() + scale_fill_viridis() + labs(title = paste0("Subject:", dimnames(allmats)$id[subj]))
}
dev.off()

rand.subs <- sample(1:length(graph_list), 5, replace = F)
pdf(paste0(conn_method, "_subj_rand_heatmap_unordered.pdf"), width=25, height=25)#, units="in", res=300)
for(subj in rand.subs){
  this.subj.adjmat <- as_adjacency_matrix(graph_list[[subj]], attr = "weight", sparse = FALSE)
  this.subj.melt <- melt(this.subj.adjmat, varnames = c("Node1", "Node2"))
  
  ggplot(this.subj.melt, aes(x=Node1, y=Node2, fill=value)) + geom_tile() + scale_fill_viridis() + labs(title = paste0("Subject:", dimnames(allmats)$id[subj]))
}
dev.off()
}

