#run and parse deltacon output for two binary graphs
#run_parse_deltacon <- function(g1, g2, workingdir=tempdir(), deltacon_binary="/Users/michael/Box_Sync/RS_BPD_graph/DeltaConv1.0/deltacon_mod", propnodes=0.8,
run_parse_deltacon <- function(g1, g2, workingdir=tempdir(), deltacon_binary="/Users/nth7/Box' 'Sync/RS_BPD_graph/DeltaConv1.0/deltacon_mod", propnodes=0.8, 
                               normalize_node_attrib=0) {
#run_parse_deltacon <- function(g1, g2, workingdir=tempdir(), deltacon_binary="/Users/michael/Box_Sync/RS_BPD_graph/DeltaConv1.0/deltacon") {
  require(igraph)
  #pass in two binary graphs as igraph objects
  
  g1edges <- get.edgelist(g1) #assumes that graphs are undirected
  g1edges[,c(1,2)] <- sub("V", "", g1edges[,c(1,2)], fixed=TRUE)
  g1edges <- cbind(g1edges, rep(1, nrow(g1edges)))
  #g1edges <- apply(g1edges, c(1,2), as.numeric) #this is slow
  mode(g1edges) <- "numeric" #faster
  g1edges <- rbind(g1edges, g1edges[,c(2,1,3)]) #transpose edges to force directed symmetric format
  
  g2edges <- get.edgelist(g2)
  g2edges[,c(1,2)] <- sub("V", "", g2edges[,c(1,2)], fixed=TRUE)
  g2edges <- cbind(g2edges, rep(1, nrow(g2edges)))
  #g2edges <- apply(g2edges, c(1,2), as.numeric)
  mode(g2edges) <- "numeric" #faster
  g2edges <- rbind(g2edges, g2edges[,c(2,1,3)]) #transpose edges to force directed symmetric format
  
  write.table(g1edges, file=file.path(workingdir, "g1edges.txt"), row.names=FALSE, col.names=FALSE)
  write.table(g2edges, file=file.path(workingdir, "g2edges.txt"), row.names=FALSE, col.names=FALSE)
  cat("g1edges.txt\ng2edges.txt\n", file=file.path(workingdir, "filelist"))
  
  #deltacon seems to require that the file list be in the current directory (doesn't support absolute paths for the filelist unless you're already in that dir). 
  #this is the segfault 11 error. To workaround, we cd into the working dir, then use the absolute path to deltacon to run
  
  doutput <- system(cmd <- paste0("cd ", workingdir, " && ", deltacon_binary, " filelist 1 naive 0.1 ", propnodes, " ", normalize_node_attrib), intern=TRUE)
  cat("Running command: ", cmd, "\n")
  #doutput <- system(paste0("cd ", workingdir, " && ", deltacon_binary, " filelist 1 naive 0.1"), intern=TRUE)

  #nodal attribution statistics come first. The occurrence of a single blank line indicates transition to edges
  blanklines <- which(doutput=="")
  stopifnot(length(blanklines) >= 2L)
  
  deltacon <- as.numeric(sub("^Similarity:\\s+([\\d\\.-]+)", "\\1", doutput[blanklines[2]+1], perl=TRUE)) #row after second blank
  deltacon_runtime <- as.numeric(sub("^DeltaCon runtime:\\s+([\\d\\.-]+)", "\\1", doutput[blanklines[2]+2], perl=TRUE)) #second row after second blank
  node_runtime <- as.numeric(sub("^Node Attribution runtime:\\s+([\\d\\.-]+)", "\\1", doutput[blanklines[2]+3], perl=TRUE)) #third row after second blank
  edge_runtime <- as.numeric(sub("^Edge Attribution runtime:\\s+([\\d\\.-]+)", "\\1", doutput[blanklines[2]+4], perl=TRUE)) #fourth row after second blank
  
  
  if (all(blanklines==c(1,2))) {
    message("Unable to compute deltacon node and edge stats")
    #nodestats <- data.frame(NA)
    nodestats <- data.frame(node=NA, attr_stat=NA)
    edgestats <- data.frame(node1=NA, node2=NA, sign=NA, attr_stat=NA)
    #edgestats <- data.frame(NA, nrow = 1, ncol = 8, colnames("node1", "node2", "sign", "attr_stat", "subj1", "subj2", "subj1.bpd", "subj2.bpd"))
    #nodestats <- NULL
    #edgestats <- NULL
  } else {
  nodestats <- doutput[1:(blanklines[1]-1)]
  nodestats <- strsplit(nodestats, "\\s+", perl=TRUE)
  nodestats <- data.frame(node=as.numeric(sapply(nodestats, "[[", 1)), attr_stat=as.numeric(sapply(nodestats, "[[", 2)))
  
  #elegant but slow
  # nodestats <- as.data.frame(do.call(rbind, gsubfn::strapply(nodestats, "^(\\d+)\\s+([\\d\\.-]+)$", function(node, attr_stat) { 
  #   return(c(node=as.numeric(node), attr_stat=as.numeric(attr_stat)))
  # })))
  
  edgestats <- doutput[(blanklines[1]+1):(blanklines[2]-1)]
  edgestats <- strsplit(edgestats, "\\s+", perl=TRUE)
  stopifnot(all(sapply(edgestats, length) == 4L)) #each row has 4 columns: node1, node2, sign, attr_stat
  
  #this is fast, but ugly (keeping it because of speed)
  edgestats <- data.frame(node1=as.numeric(sapply(edgestats, "[[", 1)), node2=as.numeric(sapply(edgestats, "[[", 2)),
                           sign=sapply(edgestats, "[[", 3), attr_stat=as.numeric(sapply(edgestats, "[[", 4)))
  
  #horribly slow
  # edgestats <- do.call(rbind, lapply(edgestats, function(el) {
  #   data.frame(node1=as.numeric(el[[1]]), node2=as.numeric(el[[2]]),
  #              sign=el[[3]], attr_stat=as.numeric(el[[4]]))
  # }))

  #clunky and requires more wranling...
  #edgestats <- do.call(rbind, edgestats)
  #colnames(edgestats) <- c("node1", "node2", "sign", "attr_stat")
  #edgestats <- as.data.frame(edgestats)
  
  #though elegant, this is way slower than the above
  #edgestats <- gsubfn::strapply(edgestats, "^(\\d+)\\s+(\\d+)\\s+([-+])\\s+([\\d\\.-]+)$", function(node1, node2, sign, attr_stat) { 
  #  return(list(node1=as.numeric(node1), node2=as.numeric(node2), sign=sign, attr_stat=as.numeric(attr_stat)))
  #})
  
  }
  return(list(deltacon=deltacon, nodestats=nodestats, edgestats=edgestats, deltacon_runtime=deltacon_runtime, node_runtime=node_runtime, edge_runtime=edge_runtime))
}





# #setwd("/Users/nth7/Box Sync/RS_BPD_graph/output.node")
# setwd("/Users/michael/Box_Sync/RS_BPD_graph/output.node")
# #save(bpd.graph10.1, control.graph10.1, file="example_binary_graphs_fordeltacon.RData")
# load(file="example_binary_graphs_fordeltacon.RData")
# orig <- run_parse_deltacon(bpd.graph10.1, control.graph10.1)
# mod <- run_parse_deltacon(bpd.graph10.1, control.graph10.1)
# identical(mod$nodestats, orig$nodestats) #checks out for 0.8 0
# mod <- run_parse_deltacon(bpd.graph10.1, control.graph10.1, normalize_node_attrib = 1, propnodes = 1)
# 
# ###
# #optimize function by profiling code
# install.packages("profvis")
# library(profvis)
# profvis(testout <- run_parse_deltacon(bpd.graph10.1, control.graph10.1))
# 
# Rprof("profout2")
# ## some code to be profiled
# xx <- replicate(100, run_parse_deltacon(bpd.graph10.1, control.graph10.1))
# Rprof(NULL)
# 
# summaryRprof("profout2")
# summaryRprof("profout")
