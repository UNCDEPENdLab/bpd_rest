#file to import and setup braingraph files
#brainGraphFromScratch
library(doParallel)
library(data.table)
library(brainGraph)
library(dplyr)

f <- Sys.getenv('PBS_NODEFILE')
nodelist <- if (nzchar(f)) readLines(f) else rep('localhost', 3)

cat("Node list allocated to this job\n")
print(nodelist)

cl_fork <- makePSOCKcluster(nodelist, outfile='')
print(cl_fork); print(unclass(cl_fork))

registerDoParallel(cl_fork)
clusterEvalQ(cl_fork, library(igraph))
clusterEvalQ(cl_fork, library(brainGraph))

pacman::p_load(plyr, ggplot2, gridExtra)
#basedir <- "/Users/mnh5174/Data_Analysis/bpd_rest"
basedir <- "/gpfs/group/mnh5174/default/Michael/bpd_rest"
setwd(basedir)

source("functions/get_subj_info.R")
source("functions/graph_util_redux.R")

#adjmats_base <- file.path(basedir, "adjmats", "402_drop2_schaefer422_fsl_pearson_prewhitened")
adjmats_base <- file.path(basedir, "adjmats")#, "schaefer421_fsl_prewhitened_pearson")
parcellation <- "schaefer421"
# conn_method <- "cor.shrink"
conn_method <- "pearson"
#preproc_pipeline <- "nosmooth_aroma_hp"
preproc_pipeline <- "fsl_prewhitened"
#these point to the no header variants used in braingraph
file_tag <- paste(parcellation, preproc_pipeline, conn_method, sep = "_")
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub=TRUE, allowCache=FALSE) #%>% filter(SPECC_ID != "037LL")

#source("scripts/estimate_euclidean_distance.R") ##creates rmShort which will delete edges close in euclidean distance

#read covariates and use Study.ID nomenclature
# covars <- fread(file.path(getMainDir(), "bpd_rest/data/SPECC_Participant_Info.csv"))
# setnames(covars, old="NUM_ID", new="Study.ID")
# setkey(covars, BPD, Study.ID)
# head(covars)

#read custom atlas
schaefer421 <- fread(file.path(basedir, "data/Schaefer_421_final_jul2018_fsl_95_group_mask_masterlookup.csv"))
schaefer421[, hemi := as.factor(hemi)] #manipulations to make brainGraph happy with the custom atlas
schaefer421[, lobe := as.factor("LOBE")] #no conceptual attachment to this, but brainGraph wants it
schaefer421$index <- 1:nrow(schaefer421)
setkey(schaefer421, index)

#add yeo communities to objects (handled inside subject loop)
membership.file <- file.path(basedir, "data", paste0(parcellation, "_membership.yeo.csv"))
membership.yeo <- as.numeric(as.matrix(read.csv(membership.file)))

##define a master membership. Use NULL as the graph object to allow this to be a general object for passage to downstream functions
yeo7_communities <- make_clusters(NULL, membership = membership.yeo, algorithm = "Yeo_etal_2011_7Networks", modularity=FALSE)

# schaefer422 <- fread(file.path(basedir, "data/schaefer422_atlas.csv"))
# schaefer422[, hemi := as.factor(hemi)] #manipulations to make brainGraph happy with the custom atlas
# schaefer422[, lobe := as.factor("LOBE")] #no conceptual attachment to this, but brainGraph wants it
# schaefer422$index <- 1:nrow(schaefer422)
# setkey(schaefer422, index)

clusterExport(cl_fork, c("schaefer421")) #export doesn't seem to be listening...

#datadir <- file.path(basedir, "adjmats/schaefer421_nosmooth_aroma_hp_cor.shrink")
datadir <- file.path("gpfs/group/mnh5174/default/MMClock/adjmats_mnh_qsub_all")

#I don't think I really need this
#init.vars <- brainGraph_init(atlas="custom", densities=seq(.08, .25, .01), custom.atlas="schaefer421", datadir="/Users/mnh5174/Data_Analysis/bpd_rest/adjmats/schaefer421_nosmooth_aroma_hp_cor.shrink", covars=covars)

#define covars based on subset of subj_info
covars <- data.table(subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan, Female) %>% 
                       mutate(
                         Group=factor(BPD, levels=c(0,1), labels=c("Control", "BPD")),
                         Female=factor(Female, levels=c(0,1), labels=c("Male", "Female")),
                         Group_all=factor("All") #for consensus all approach
                       ) %>% dplyr::select(-BPD)
)

#saveRDS(covars, file = file.path(basedir, "cache", "covars_prewhite402_drop2_84.rds")) #for the global checks, dont need to do this every time
setnames(covars, old="SPECC_ID", new="Study.ID")
setkey(covars, Group, Study.ID)
head(covars)

#thresh_methods <- c("consensus_all", "consensus_all_25", "consensus_all_30", "consensus_all_25_30")
thresh_methods <- c("consensus_all_25_dense", "consensus_all_25_40")
#thresh_methods <- c("fc") #others already complete
#thresh_methods <- c("consensus_all_25", "consensus_all", "consensus_all_25_30", "consensus_all_30")
#thresh_methods <- c("consistency", "consensus_all") #, "mean") #others have already completed
#thresh_methods <- c("consensus_all") #, "mean") #others have already completed

#corresponding thresholds for each method
if(preproc_pipeline == "fsl_prewhitened"){
  #prewhitened data has much lower FC values
  
  mat_thresh <- list(
    #density=seq(.1, .25, .01), #density
    #consensus=pracma::logspace(log10(.1), log10(.44), 20), #consensus (by group) works with FC values
    #fc=pracma::logspace(log10(.1), log10(.44), 20), #we can hack a raw FC threshold if we use sub.thresh = 0
    #consistency=seq(.1, .25, .01), #consistency works with densities
    #consensus_all=pracma::logspace(log10(.1), log10(.44), 20), #consensus_all (one group) works with FC values
    #consensus_all_25=pracma::logspace(log10(.1), log10(.44), 20), #variant of consensus all with 25% presence threshold
    #consensus_all_30=pracma::logspace(log10(.1), log10(.44), 30), #consensus_all (one group) works with FC values
    #consensus_all_25_30=pracma::logspace(log10(.1), log10(.44), 30) #variant of consensus all with 25% presence threshold
    consensus_all_25_dense = pracma::logspace(log10(.1), log10(.25), 15),
    consensus_all_25_40 = pracma::logspace(log10(.1), log10(.44),40)
  )
} else {
  
  mat_thresh <- list(
    density=seq(.1, .25, .01), #density
    consensus=pracma::logspace(log10(.44), log10(.58), 20), #consensus (by group) works with FC values
    fc=pracma::logspace(log10(.44), log10(.58), 20), #we can hack a raw FC threshold if we use sub.thresh = 0
    consistency=seq(.1, .25, .01), #consistency works with densities
    consensus_all=pracma::logspace(log10(.44), log10(.58), 20), #consensus_all (one group) works with FC values
    consensus_all_25=pracma::logspace(log10(.44), log10(.58), 20), #variant of consensus all with 25% presence threshold
    consensus_all_30=pracma::logspace(log10(.12), log10(.58), 30), #consensus_all (one group) works with FC values
    consensus_all_25_30=pracma::logspace(log10(.12), log10(.58), 30) #variant of consensus all with 25% presence threshold
  )  
}


sub_thresh <- 0.5 #keep connections present in at least 50% of subjects

#need to strip header from cor mats for braingraph import: find . -iname "*.txt" | xargs -I {} -n 1 sed -i '' '1d' {}
inds <- list(BPD=which(subj_info$BPD==1), Control=which(subj_info$BPD==0))
inds_all <- list(all=1:nrow(subj_info)) #for the consensus_all approach in which the edge must be present in > 0.5 of all subjects
today <- format(Sys.Date(), '%Y-%m-%d')
# savedir <- file.path(basedir, "braingraph_objects", "20_thresh_prewhite402")
savedir <- file.path(basedir, "braingraph_objects", "prewhite402_drop")
atlas <- "schaefer421"
groups <- c("BPD", "Control")
#roiMat <- read.table(file.path(basedir,  "data", "schaefer422_roiMat.txt"), col.names=c("x.mni", "y.mni", "z.mni", "roinum"))
roiMat <- read.csv(file.path(basedir, "data", "Schaefer_421_final_jul2018_fsl_95_group_mask_roiMat.csv"), header = TRUE); colnames(roiMat)[4] <- "roinum"

#function to compute a 1/0 matrix that handles ROIs who CoMs fall within a given distance of each other (default 20mm)
#expect input as data.frame or matrix where columns 1, 2, 3 contain coordinates in x, y, z, respectively
remove_short_dist <- function(roiMat, dist_thresh = 20) {
  nnodes <- nrow(roiMat)
  roiDist <- matrix(NA, nrow=nnodes, ncol=nnodes) 
  
  for (i in 1:nnodes) {
    for (j in 1:nnodes) {
      #populate lower triangle only for speed, then transpose onto upper triangle
      if (i > j) { roiDist[i,j] <- sqrt((roiMat[i,1] - roiMat[j,1])^2 + (roiMat[i,2] - roiMat[j,2])^2 + (roiMat[i,3] - roiMat[j,3])^2) }
    }
  }
  
  #copy lower triangle to upper triangle
  roiDist[upper.tri(roiDist)] <- t(roiDist)[upper.tri(roiDist)]
  diag(roiDist) <- 0
  
  ####quick QA
  # hist(roiDist)
  # vecDist <- roiDist[lower.tri(roiDist)]
  # sum(vecDist < 20)/length(vecDist)
  
  #roi.dist can be changed at the front end of the script
  rmShort <- roiDist > dist_thresh
  #creates binary matrix, in which 0 denotes a short distanced connection that is to be removed. 
  rmShort <- apply(rmShort, c(1,2), as.numeric)

  return(rmShort)
}

#15mm removes 631 edges out of a possible 88410 (about 0.7%)
#this seems like a good call since many of the censored nodes are the subnuclei of amygdala, thalamus, striatum
rm15 <- remove_short_dist(roiMat, dist_thresh = 15)
sum(rm15[lower.tri(rm15)] < 1)
#which(rm15==0, arr.ind=TRUE)

fix_neg_edges <- function(A, rmShort=NULL) {
  require(abind)
  #A is a list of adjacency matrices by threshold. Each element is a 3-d array ROIs x ROIs x subjects
  #This function removes negative edges (since these are not tolerated in many graph metrics)
  #It also records how many positive, negative, and zero edges are present. Zero edges likely reflect the thresholding method (create_mats).
  #rmShort is an options ROIs x ROIs 1/0 matrix that multiplies each subject's adjancency matrix. Thus, if it is passed in, this censors short connections.
  
  lapply(A, function(thresh_array) {

    mats <- apply(thresh_array, 3, function(mat) {
      mat[mat < 0] <- 0
      return(mat)
    })
    dim(mats) <- dim(thresh_array) #apply screws up the dimensions of the return -- fix them

    if (is.null(rmShort)) {
      #don't remove connections based on distance
      rmShort <- array(1, dim=dim(thresh_array))
    } else {
      #duplicate ROI x ROI censoring along the subjects dimension
      rmShort <- abind(replicate(dim(thresh_array)[3], rmShort, simplify=FALSE), along=3) 
    }

    mats <- mats * rmShort #apply short-range censoring, if specified
    
    #use machine precision to ensure that we are binning floating point numbers properly near zero
    attr(mats,"n_neg_edges") <- apply(thresh_array, 3, function(mat) { sum(mat[lower.tri(mat)] < -sqrt(.Machine$double.eps)) })
    attr(mats,"n_pos_edges") <- apply(thresh_array, 3, function(mat) { sum(mat[lower.tri(mat)] > sqrt(.Machine$double.eps)) })
    attr(mats,"n_zero_edges") <- apply(thresh_array, 3, function(mat) { sum(abs(mat[lower.tri(mat)]) < sqrt(.Machine$double.eps)) })
    
    return(mats)
  })
}

#import raw adjacency matrices
for (tm in 1:length(thresh_methods)) {
  dir.create(file.path(savedir, thresh_methods[tm], "weighted"), recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(savedir, thresh_methods[tm], "binary"), recursive=TRUE, showWarnings=FALSE)
  thresholds <- mat_thresh[[ thresh_methods[tm] ]]

#  if (thresh_methods[tm] == "consensus_all" || thresh_methods[tm] == "consensus_all_25") {
  if(grepl("consensus_all", thresh_methods[tm], fixed = TRUE)){   #allows for detection of the 30 threshold cases
    groups <- "All" #just one
    inds <- inds_all #use the no group structure below
    setkey(covars, Group_all, Study.ID) #to get the subject id lookup below to work
    #if (thresh_methods[tm] == "consensus_all_25") { ss <- 0.25 #use tmp var to avoid overwriting sub_thresh if operating in a loop over methods
    if (grepl("consensus_all_25", thresh_methods[tm], fixed = TRUE)) { ss <- 0.25
    } else { ss <- sub_thresh }
    allmats <- create_mats(subj_info$file, modality="fmri", threshold.by="consensus", mat.thresh=thresholds, sub.thresh=ss, inds=inds_all)
  } else if (thresh_methods[tm] == "fc") {
    allmats <- create_mats(subj_info$file, modality="fmri", threshold.by="consensus", mat.thresh=thresholds, sub.thresh=0, inds=inds) #use 0 similarity threshold to get raw FC
  } else {
    allmats <- create_mats(subj_info$file, modality="fmri", threshold.by=thresh_methods[tm], mat.thresh=thresholds, sub.thresh=sub_thresh, inds=inds)
  }    

  g.group <- g <- fnames <- vector('list', length=length(groups))
  
  allmats$A.norm.sub.noneg <- fix_neg_edges(allmats$A.norm.sub) #remove negative edges (usually only a few) before calculating graphs
  allmats$A.norm.sub.noneg.rmshort <- fix_neg_edges(allmats$A.norm.sub, rmShort=rm15) #remove negative edges (usually only a few) and short connections before calculating graphs

  #currently applying the rmShort logic to all graphs
  A.norm.sub.noneg.rmshort <- allmats$A.norm.sub.noneg.rmshort

  #this is the right idea (using the noneg from subjects, but it would mean parsing out the groups myself. and I don't really need the mean mats anyhow)
  #all.mean <- lapply(A.norm.sub.noneg.rmshort, function(x) { rowMeans(x, dims=2) })
  A.norm.mean <- lapply(allmats$A.norm.mean, function(glist) { lapply(glist, function(g) { g[g < 0] <- 0; return(g) }) })

  #single subject fitting
  for (i in seq_along(groups)) {
    for (j in seq_along(thresholds)) {
      print(paste0('Method ', thresh_methods[tm], '; threshold ', j, '/', length(thresholds), '; group ', i, '; ', format(Sys.time(), '%H:%M:%S')))

      #compute metrics on binary and weighted variants
      foreach (k=seq_along(inds[[i]]), .noexport=c("allmats")) %dopar% { #.export=c("schaefer421", "covars")
        thisadj <- A.norm.sub.noneg.rmshort[[j]][, , inds[[i]][k]] #the norm.sub.noneg.rmshort field has the positive thresholded correlation matrix without short-range connections
        
        g.tmp <- graph_from_adjacency_matrix(thisadj, mode='undirected', diag=FALSE, weighted=TRUE)
        g.tmp <- set_vertex_attr(g.tmp, "community", value=yeo7_communities$membership) #add community attribute to vertices to get module-based statistics
        g.bin <- remove.edge.attribute(g.tmp, "weight")

        g.tmp <- set_brainGraph_attr(g.tmp, atlas, modality='fmri',
          weighting='pearson', threshold=thresholds[j],
          subject=covars[groups[i], Study.ID[k]], group=groups[i],
          use.parallel=FALSE, A=thisadj, clust.method="apriori")

        #wi bw strength
        wibw <- wibw_module_degree(g.tmp, community_attr="community")

        V(g.tmp)$wi_strength_z <- wibw$z_within
        V(g.tmp)$bw_strength_z <- wibw$z_between
        V(g.tmp)$wi_strength <- wibw$Ki
        V(g.tmp)$bw_strength <- wibw$Bi
                
        g.tmp.bin <- set_brainGraph_attr(g.bin, atlas, modality='fmri',
          weighting='pearson', threshold=thresholds[j],
          subject=covars[groups[i], Study.ID[k]], group=groups[i],
          use.parallel=FALSE, clust.method="apriori")

        #wi bw degree
        wibw <- wibw_module_degree(g.tmp.bin, community_attr="community")

        V(g.tmp.bin)$wi_degree_z <- wibw$z_within
        V(g.tmp.bin)$bw_degree_z <- wibw$z_between
        V(g.tmp.bin)$wi_degree <- wibw$Ki
        V(g.tmp.bin)$bw_degree <- wibw$Bi
        
        saveRDS(g.tmp.bin, file=file.path(savedir, thresh_methods[tm], "binary", paste0(sprintf('g%i_thr%02i_subj%03i%s', i, j, k, '.rds'))))
        saveRDS(g.tmp, file=file.path(savedir, thresh_methods[tm], "weighted", paste0(sprintf('g%i_thr%02i_subj%03i%s', i, j, k, '.rds'))))
      }
    }
    
    # Group mean weighted graphs
    print(paste0('Group', i, '; ', format(Sys.time(), '%H:%M:%S')))
    g.group[[i]] <- lapply(seq_along(thresholds), function(x)
      graph_from_adjacency_matrix(A.norm.mean[[x]][[i]], mode='undirected', diag=FALSE, weighted=TRUE))
    
    g.group[[i]] <- llply(seq_along(thresholds), function(x) {
      set_brainGraph_attr(g.group[[i]][[x]], atlas, modality='fmri', weighting='pearson',
        threshold=thresholds[x], group=groups[i],
        A=A.norm.mean[[x]][[i]], use.parallel=FALSE) }, .parallel=FALSE)

  }

  #save these
  saveRDS(allmats, file=file.path(savedir, thresh_methods[tm], paste0("allmats_", thresh_methods[tm], ".rds")))
  rm(allmats)

  #aggregate graphs into sample-level objects
  for (gtype in c("binary", "weighted")) {
    
    for (i in seq_along(groups)) {
      g[[i]] <- fnames[[i]] <- vector('list', length=length(thresholds)) 
      for (j in seq_along(thresholds)) {
        fnames[[i]][[j]] <- list.files(file.path(savedir, thresh_methods[tm], gtype), sprintf('g%i_thr%02i.*', i, j), full.names=T)
        g[[i]][[j]] <- lapply(fnames[[i]][[j]], readRDS)
        #add missing ID variable... capitalization problem! don't want to wait to recompute
        #g[[i]][[j]] <- lapply(fnames[[i]][[j]], function(f) {
        #  snum <- as.numeric(sub(".*_subj(\\d+).rds", "\\1", f, perl=TRUE))
        #  #Study.ID has to be char for some subsidiary functions to run
        #  x <- readRDS(f) %>% set_graph_attr("name", as.character(covars[groups[i], Study.ID[snum]])) %>%
        #    set_graph_attr("Group", groups[i])
        #})
      }
      
      x <- all.equal(sapply(g[[i]][[1]], graph_attr, 'name'),
        covars[groups[i], Study.ID])
      cat("match between IDs in graphs and covars is: ", x, "\n")
      #browser()
      # if (isTRUE(x)) lapply(fnames[[i]], file.remove) 
    }
  
    saveRDS(g, file=file.path(savedir, thresh_methods[tm], paste0(thresh_methods[tm], '_g_', gtype, '.rds')))
    saveRDS(g.group, file=file.path(savedir, thresh_methods[tm], paste0(thresh_methods[tm], '_g_', gtype, '.group.rds')))

  }

}
stopCluster(cl_fork)
