membership_df <- get(load("/Users/nth7/Box Sync/bpd_rest/cache_protected/membership_df.RData"))
community <- membership_df$membership

nodes.file = file.path(resultsdir,  "nodal_stats_nbs_low_df.csv")
nodevals <- read.csv(nodes.file, row.names = NULL)
nodevals <- as.numeric(nodevals$strength)
nodestp <- read.csv(nodes.file, row.names = NULL)
nodestp <- as.character(nodestp$node)

nodevals 

output <- nodevals[,c(2,4,5,6)]


atlas_exp <- atlas[,1:4]
names(atlas_exp) <- c("node", "x","y","z")

more_exp <- atlas_exp %>% inner_join(output, by = "node")

more_exp <- left_join(more_exp, x, by = "node")


exp_file <- more_exp[,2:4]
exp_file$net <- more_exp$net
exp_file$strength <- more_exp$strength
exp_file$lab <- more_exp$node


write.table(exp_file, file = "nbs_results_with_t-weighted_strength.node", col.names = FALSE, row.names = FALSE)

x$node <- paste0("V",seq(1,422,1))
factor(membership_df$membership)


NodeFile <- function(atlas, community = NULL, nodestp = NULL, nodevals = NULL, nnodes, labels, filename,outputdir){
  #nodestp provides vector of node numbers that need to be plotted
  #community can be set to community membership, whether higher in BPD/ controls, or anything else  
  #labels set equal to either 1 (full anatomical label) or 2 (node number [i.e. V_1...])
  # browser()
  nf <- atlas
  
  #if want to read in a csv rather than a predfined atlas in R
  #nf <- read.csv(atlas, header = TRUE)        #####csv file should have at least mni attributes and anatomical labels
  vnames <- data.frame(atlas[,"name"])
  names(vnames) = "vname"
  stopifnot(all(c("x.mni", "y.mni", "z.mni", "anat_label") %in% names(nf)))
  # row.names(nf) <- c(seq(1, 248, 1), seq(251, 271, 1))
  row.names(nf) <- sort(as.numeric(gsub("V", "", vnames[,1])))
  
  if(!is.null(community)){
    commvals <- community
    names(commvals) <- vnames$vname
  } else {
    commvals = rep(1,length(vnames[,1]))
    names(commvals) <- vnames$vname
  }
  
  
  if(!is.null(nodestp)) {
    nf <- nf[which(nf$name %in% nodestp),]
    commvals <- commvals[which(names(commvals) %in% nodestp)]
    # nf$anat_label <- paste0("V_", nf$vname)
    # nf <- nf %>% select(-vname)
    #names(nodevals) <- vnames$vname
    
    if(is.null(nodevals)){
      nodevals <- rep(1, length(nodestp))
    } else {
      nodevals <- nodevals[which(names(nodevals) %in% nodestp)]
    }
  } else  {
    if(is.null(nodevals)){
      nodevals <- rep(1, nnodes)
    } else {
      nodevals <- nodevals
    }
  }
  
  mni <- subset(nf, select = c(x.mni, y.mni, z.mni))
  if(labels == 1){
    node_label <- subset(nf, select = anat_label) 
  } else if(labels == 0) {node_label <- paste0("V", row.names(nf))}
  
  nf.exp <- cbind(mni, commvals, nodevals, node_label)
  
  write.table(nf.exp, file = file.path(outputdir, paste0(filename, ".node.txt")), row.names = FALSE, col.names = FALSE)
  return(nf.exp)
}
