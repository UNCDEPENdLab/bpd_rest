#setwd("~/Box Sync/RS_BPD_graph")
setwd("/Users/michael/Data_Analysis/bpd_rest")
basedir <- getwd()

#this will setup details of the parcellation, conn_method, preproc_pipeline, and connection distance
#it also sources all helper functions for additional analysis
source("functions/setup_globals.R") 

#output for node files
odir=file.path(basedir, "BNV_nodefiles", "michael")

#load weighted adjacency matrices
allmats <- load_allmats()

comm_weighted_louvain <- run_community_detection_on_agg(allmats, "louvain")
sizes(comm_weighted_louvain)
nf(atlas, comm_weighted_louvain, fname="louv_weighted_n83.node", outputdir=odir, bycomm=TRUE, savebnv=TRUE)

#warning: savebnv is slow!
#louvain on densities from 3-20%
for (d in seq(.03, .20, .01)) {
  comm_l <- run_community_detection_on_agg(allmats, "louvain", density=d)
  thiso <- file.path(odir, paste0("louv_n83_d", d))
  dir.create(thiso)
  nf(atlas, comm_l, fname=paste0("louv_d", d, "_n83.node"), outputdir=thiso, bycomm=TRUE, savebnv=TRUE)
}

#warning: savebnv is slow!
#infomap on densities from 3-20%
for (d in seq(.03, .20, .01)) {
  comm_l <- run_community_detection_on_agg(allmats, "infomap", density=d)
  thiso <- file.path(odir, "infomap", paste0("infomap_n83_d", d))
  dir.create(thiso, recursive=TRUE)
  nf(atlas, comm_l, fname=paste0("infomap_d", d, "_n83.node"), outputdir=thiso, bycomm=TRUE, savebnv=TRUE)
}