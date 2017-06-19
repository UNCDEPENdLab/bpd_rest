setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/BNV_nodefiles/final_louv_d10/")
c1 <- read.table("louv_d10_comm1.node.txt")
c2 <- read.table("louv_d10_comm2.node.txt")
c3 <- read.table("louv_d10_comm3.node.txt")
c4 <- read.table('louv_d10_comm4.node.txt')
c5 <- read.table("louv_d10_comm5.node.txt")

mean.graph <- apply(allmats, c(2,3), mean, na.rm = TRUE)
mean.g <- graph.adjacency(mean.graph, mode = "lower", weighted = TRUE, diag = FALSE)
mean.g <- delete.edges(mean.g, which(E(mean.g)$weight < 0))
mean.g.adj <- as_adjacency_matrix(mean.g, attr = "weight")


membership.bynode <- rbind(c1[,c(4,6)], c2[,c(4,6)], c3[,c(4,6)], c4[,c(4,6)], c5[,c(4,6)])
membership.bynode[,2] <- as.numeric(gsub("V", "", membership.bynode[,2]))
membership <- membership.bynode[order(membership.bynode[,2]), 1]
names(membership) <- c(paste0("V", seq(1,248,1)),paste0("V", seq(251,271,1)))
                                
d10_louv <- make_clusters(mean.g, membership, algorithm = "louvain")

saveRDS(d10_louv, file = "d10_louv_MH.rds")
