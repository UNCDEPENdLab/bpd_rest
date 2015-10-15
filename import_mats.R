setwd(file.path(getMainDir(), "bpd_rest"))

corrfiles=system("find /Volumes/Serena/Raj/Preprocess_Rest -iname corr_rois_pearson_new_r_v2.txt -type f", intern=TRUE)
ids <- sub("/Volumes/Serena/Raj/Preprocess_Rest/(?:SPECC/)*([^/]+)/corr_rois_pearson_new_r_v2.txt$", "\\1", corrfiles, perl=TRUE)
bpd <- grep("\\d{3}[A-z]{2}_.*", ids, perl=TRUE)
bpd <- bpd[!ids[bpd] %in% c("001RA_07DEC2013", "005AI_06NOV2013", "023ds_07May2014", "050ai_06Nov2014", "0531lw_16Dec2014")]
controlmatch <- read.table("neil/demographics/limited_controls.txt")$V1
controls <- which(ids %in% controlmatch)

allmats <- array(NA, c(length(corrfiles), 264, 264), dimnames=list(id=ids, NULL, NULL))

for (f in 1:length(corrfiles)) {
  m <- as.matrix(read.table(corrfiles[f]))
  allmats[f,,] <- m
}

#compute averaged adjacency matrix
avgmat <- apply(allmats, c(2), function(x) { colMeans(x, na.rm=TRUE) })
hist(avgmat[lower.tri(avgmat)])
summary(avgmat[lower.tri(avgmat)])

#create graphs from adjacency matrices
library(igraph)

#original graphs with correlation weights 
allg <- apply(allmats, 1, function(sub) {
      graph.adjacency(sub, mode="undirected", weighted=TRUE, diag=FALSE)
    })

#remove negative correlations since these are extremely uncommon (without global signal regression) and not handled well by many weighted metrics 
allg_noneg <- lapply(allg, function(g) {
      delete.edges(g, which(E(g)$weight < 0))
    })

#binarize graphs at densities ranging from 1-20%
densities_desired <- seq(.01, .2, .01) #1-20%

allg_density <- lapply(allg_noneg, function(g) {
      nnodes <- length(V(g))
      maxedges <- (nnodes*(nnodes-1))/2
      
      #compute densities for all correlation thresholds up to the max
      maxcorr <- round(max(E(g)$weight)*100) + 1
      densities_observed <- sapply(0:maxcorr, function(x) { sum(E(g)$weight*100 >= x)/maxedges })      
      
      #interpolate the correlation threshold to achieve the precise density required
      dgraphs <- lapply(densities_desired, function(d) {
            ddiscrep <- densities_observed - d
            above_density <- min(densities_observed[ddiscrep > 0])
            below_density <- max(densities_observed[ddiscrep < 0])
            above_corr <- (which(densities_observed==above_density) - 1)/100 #subtract 1 since the 0:maxcorr above is zero-based
            below_corr <- (which(densities_observed==below_density) - 1)/100 
            thresh_corr = below_corr + ((d - below_density)/(above_density - below_density))*(above_corr - below_corr)
            gthresh <- delete.edges(g, which(E(g)$weight <= thresh_corr)) #return binary weighted density-thresholded graph
            gthresh <- remove.edge.attribute(gthresh, "weight")
          })
      #lapply(dgraphs, edge_density) #just for verification
      names(dgraphs) <- paste0("d", densities_desired)
      return(dgraphs)
    })

#each element of allg_density is a list of 20 binary graphs for that subject at 1-20% density

#for basic checks, compute degree, eigenvector centrality, local clustering, betweenness, and pagerank at each density threshold
all_centrality <- lapply(allg_density, function(subj) {
      lapply(subj, function(dgraph) {
            #each of these returns a nnodes-length vector of values 
            deg <- degree(dgraph)
            evc <- evcent(dgraph)$vector
            bet <- betweenness(dgraph, normalize=TRUE)*100
            locclust <- transitivity(dgraph, type="local")
            pr <- page.rank(dgraph, algo="prpack")$vector
            cbind(degree=deg, evcent=evc, betweenness=bet, locclust=locclust, pagerank=pr) #return a 264 x nmetrics matrix for analysis/reduction            
          })
    })

#nodal analysis of centrality measures at each density
controlstats <- all_centrality[controls]
bpdstats <- all_centrality[bpd]
for (d in 1:length(densities_desired)) {
   #conduct two-sample t-tests for each node and each centrality measure and each density (ugly hack!)
   #build a matrix of subjects and nodes
   bpdmat <- do.call(rbind, lapply(1:length(bpdstats), function(subnum) { data.frame(id=names(bpdstats)[subnum], bpd=1, node=1:264, bpdstats[[subnum]][[d]]) }))
   controlmat <- do.call(rbind, lapply(1:length(controlstats), function(subnum) { data.frame(id=names(controlstats)[subnum], bpd=0, node=1:264, controlstats[[subnum]][[d]]) }))
   
   metrics <- names(bpdmat)[!names(bpdmat) %in% c("id", "bpd", "node")]
   abst <- 3.0 #about p = .005
   for (m in metrics) {
     for (n in unique(bpdmat$node)) {
       bpdvec <- bpdmat[bpdmat$node==n, m]
       controlvec <- controlmat[controlmat$node==n, m]
       if (any(is.na(bpdvec)) || any(is.na(controlvec))) {
         next #should print a warning here with details
       }
       res <- t.test(bpdvec, controlvec)
       if (!is.nan(res$statistic) && abs(res$statistic) > abst) {
         cat("For density ", densities_desired[d], ",", m, ", node:", n, "group t =", round(res$statistic, 3), "p =", round(res$p.value, 3), "means are:", round(res$estimate, 3), "\n")
       }
     }
   }
}

#also a quick brute force edge analysis on fisher-transformed correlations
#allmats <- atanh(allmats) #fisher z transform

edgepvals <- matrix(NA, nrow=264, ncol=264)
abst <- 4.3 #about p = .0001 (two-tailed)
for (i in 1:264) {
  for (j in 1:264) {
    if (i > j) {
      #lower triangle
      controlcorrs <- allmats[controls,i,j]
      bpdcorrs <- allmats[bpd, i, j]
      res <- t.test(atanh(controlcorrs), atanh(bpdcorrs))
      edgepvals[i,j] <- res$p.value
      edgepvals[j,i] <- res$statistic #t stat on upper triangle
      if (!is.nan(res$statistic) && abs(res$statistic) > abst) {
        cat("Significant edge difference: ", i, j, ", t=", round(res$statistic, 3), "p =", round(res$p.value, 3), "means are:", round(res$estimate, 5), "\n")
      }
    }
  }
}

#try fdr correction
pvec <- edgepvals[lower.tri(edgepvals)]
#tvec <- edgepvals[upper.tri(edgepvals)]
lookup <- which(lower.tri(edgepvals) == TRUE, arr.ind=TRUE)
padj <- p.adjust(pvec, "fdr")
lookup[which(pvec < .0001),]
padj[which(padj < .05)]
lookup[which(padj < .05),]