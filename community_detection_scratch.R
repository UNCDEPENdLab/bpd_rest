###setup community
##louvain
comm_weighted_louvain <- run_community_detection_on_agg(allmats, "louvain")
comm_d05_l <- run_community_detection_on_agg(allmats, "louvain", density=0.05)
comm_d08_l <- run_community_detection_on_agg(allmats, "louvain", density=0.08)
comm_d10_l <- run_community_detection_on_agg(allmats, "louvain", density=0.10)
comm_d15_l <- run_community_detection_on_agg(allmats, "louvain", density = 0.15)
comm_d20_l <- run_community_detection_on_agg(allmats, "louvain", density=0.20)

allg_noneg_w_louvain <- assign_communities(allg_noneg, comm_weighted_louvain, "comm_weighted_louvain")
allg_noneg_d05_louvain <- assign_communities(allg_noneg, comm_d05_l, "comm_d05_l")
allg_noneg_d10_louvain <- assign_communities(allg_noneg, comm_d10_l, "comm_d10_l")
allg_noneg_d15_louvain <- assign_communities(allg_noneg, comm_d15_l, "comm_d15_l")
allg_noneg_d20_louvain <- assign_communities(allg_noneg, comm_d20_l, "comm_d20_l")

##infomap
comm_weighted_infomap <- run_community_detection_on_agg(allmats, "infomap")
comm_d05_i <- run_community_detection_on_agg(allmats, "infomap", density=0.05)
comm_d10_i <- run_community_detection_on_agg(allmats, "infomap", density=0.10)
comm_d15_i <- run_community_detection_on_agg(allmats, "infomap", density = 0.15)
comm_d20_i <- run_community_detection_on_agg(allmats, "infomap", density=0.20)

allg_noneg_w_infomap <- assign_communities(allg_noneg, comm_weighted_infomap, "comm_weighted_infomap")
allg_noneg_d05_infomap <- assign_communities(allg_noneg, comm_d05_i, "comm_d05_i")
allg_noneg_d10_infomap <- assign_communities(allg_noneg, comm_d10_i, "comm_d10_i")
allg_noneg_d15_infomap <- assign_communities(allg_noneg, comm_d15_i, "comm_d15_i")
allg_noneg_d20_infomap <- assign_communities(allg_noneg, comm_d20_i, "comm_d20_i")


##inspect community layout:
table(get.vertex.attribute(allg_noneg_w_louvain$`001RA`)$comm_weighted_louvain)
table(get.vertex.attribute(allg_noneg_d05_louvain$`001RA`)$comm_d05_l)
table(get.vertex.attribute(allg_noneg_d10_louvain$`001RA`)$comm_d10_l)
table(get.vertex.attribute(allg_noneg_d15_louvain$`001RA`)$comm_d15_l)
table(get.vertex.attribute(allg_noneg_d20_louvain$`001RA`)$comm_d20_l)
table(get.vertex.attribute(allg_noneg_w_infomap$`001RA`)$comm_weighted_infomap)
table(get.vertex.attribute(allg_noneg_d05_infomap$`001RA`)$comm_d05_i)
table(get.vertex.attribute(allg_noneg_d10_infomap$`001RA`)$comm_d10_i)
table(get.vertex.attribute(allg_noneg_d15_infomap$`001RA`)$comm_d15_i)
table(get.vertex.attribute(allg_noneg_d20_infomap$`001RA`)$comm_d20_i)


# ##fastgreedy
# comm_weighted_fast_greedy <- run_community_detection_on_agg(allmats, "fast_greedy")
# comm_d05_fg <- run_community_detection_on_agg(allmats, "fast_greedy", density=0.05)
# comm_d10_fg <- run_community_detection_on_agg(allmats, "fast_greedy", density=0.10)
# comm_d15_fg <- run_community_detection_on_agg(allmats, "fast_greedy", density = 0.15)
# comm_d20_fg <- run_community_detection_on_agg(allmats, "fast_greedy", density=0.20)
# 
# allg_noneg_w_fast_greedy <- assign_communities(allg_noneg, comm_weighted_fast_greedy, "comm_weighted_fast_greedy")
# allg_noneg_d05_fast_greedy <- assign_communities(allg_noneg, comm_d05_fg, "comm_d05_fg")
# allg_noneg_d10_fast_greedy <- assign_communities(allg_noneg, comm_d10_fg, "comm_d10_fg")
# allg_noneg_d15_fast_greedy <- assign_communities(allg_noneg, comm_d15_fg, "comm_d15_fg")
# allg_noneg_d20_fast_greedy <- assign_communities(allg_noneg, comm_d20_fg, "comm_d20_fg")

allg_noneg <- allg_noneg_d05_louvain

setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/")

for(net in 1:length(unique(get.vertex.attribute(allg_noneg$`001RA`)$comm_d05_l))){
  allg_noneg_comm <- get.vertex.attribute(allg_noneg$`001RA`)$name[which(get.vertex.attribute(allg_noneg$`001RA`)$comm_d05_l == net)]
  
  node.file <- NodeFile(atlas = atlas, 
                      #community = rep(2, nnodes),
                      nnodes = nnodes,
                      nodestp = allg_noneg_comm,
                      filename = paste0(parcellation, "_", preproc_pipeline, "_", conn_method, "_comm_d05_l",net),
                      labels = 0,
                      outputdir = paste0(getwd(), "/BNV_nodefiles/in_consideration/finalists/pcor.ridge_d05_louv")
  )
}
