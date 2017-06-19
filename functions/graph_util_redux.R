###################################################################################
##Graph Utility Functions: NodeFile(), tagGraph(), plotMetricQuant(), getMotionInfo()
##MotionInfoSPEC: function pulls motion info from specified dir and exports a list containing: 
#directories not found, fd.txt files not found, raw FD data, and a table with relevant FD information based on specifications

filter_movement <- function(subj_info, dir, thresh=0.5, spikes.exclude=0.2, max.exclude=10) {
  ##thresh indicates what specifies a significant head movement (in mm)
  ##spikes.exclude specifies the cutoff for the proportion of movements above the threshold at which to exclude a given subject
  ##max.exclude is similar to spikes.exclude, but sets a cutoff for any large head movements above a set point
  require(plyr); require(dplyr); require(tidyr)
  options(dplyr.width=200)
  return_list <- list()
  
  ## MH has now converted all SPECC MR directory names to all lower case to allow for match on case-sensitive filesystem and to make the naming consistent
  
  ##pull subject folders
  subj_info <- subj_info %>% select(-Notes) %>% rowwise() %>% mutate(mr_dir=ifelse(LunaMRI==1,
          paste0(dir, "/MMClock/MR_Proc/", Luna_ID, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d")), #convert to Date, then reformat YYYYMMDD
          paste0(dir, "/SPECC/MR_Proc/", tolower(SPECC_ID), "_", tolower(format((as.Date(ScanDate, format="%Y-%m-%d")), "%d%b%Y")))))
  
  ##pull fd info into dataframe
  fd.info.raw <- data.frame()
  fd.info <- data.frame()
  fd.notfound <- data.frame()
  for(r in 1:nrow(subj_info)) {
    thissub <- subj_info[r,]
    if(file.exists(paste0(thissub$mr_dir, "/mni_5mm_", preproc_pipeline,"/rest1/motion_info/fd.txt"))){
      this.sub.fd <- data.frame(Subj=as.character(thissub$SPECC_ID), FD=read.table(paste0(thissub$mr_dir, "/mni_5mm_", preproc_pipeline, "/rest1/motion_info/fd.txt"))$V1)
      
      meanFD = mean(this.sub.fd$FD)
      maxFD = max(this.sub.fd$FD)
      max.exclude.this = if (maxFD > max.exclude) 1 else 0
      prop.thresh = sum(this.sub.fd$FD > thresh)/length(this.sub.fd$FD)
      thresh.exclude.this = if (prop.thresh > spikes.exclude) 1 else 0
      num_volumes = length(this.sub.fd$FD)
      
      thissub.desc <- data.frame(Subj = as.character(thissub$SPECC_ID), meanFD, maxFD, max.exclude.this, prop.thresh, thresh.exclude.this, num_volumes, stringsAsFactors=FALSE)
      
      if(!(num_volumes==300)) { cat("Subject ", as.character(thissub$SPECC_ID), " does not have the standard number of volumes", "\n") }
      ##subjects 003 and 008 have been taken care of 6/14/17 - NH
      
      fd.info <- rbind(fd.info, thissub.desc)
      fd.info.raw <- rbind(fd.info.raw, this.sub.fd)
    } else{
      fd.notfound <- rbind(fd.notfound, data.frame(thissub))
      cat("FD.txt file not found for subject: ", as.character(thissub$SPECC_ID), "\n")
    }
  }
  
  #verify that mr_dir is present as expected, if not, store in dir.notfound
  subj_info$dirfound <- file.exists(subj_info$mr_dir)
  if (!all(subj_info$dirfound)) {
    dir.notfound <- subset(subj_info, dirfound==FALSE)
    cat("MR_dir not found for the following subjects:\n\n")
    print(dir.notfound)
  }
  
  fd.exclude <- fd.info[,c("Subj", "max.exclude.this", "thresh.exclude.this")]
  colnames(fd.exclude) <- c("SPECC_ID", "fd.max.exclude", "fd.thresh.exclude")

  subj_info_exclude <- inner_join(subj_info, fd.exclude, by = "SPECC_ID")
  subj_info_exclude <- filter(subj_info_exclude, fd.max.exclude == 0 & fd.thresh.exclude == 0)
  return(subj_info_exclude)
}
#a <- MotionInfoSPECC(lab_ics, 0.5, .20, 10)


##NodeFile: function outputs node file readable by BrainNet Viewer in MATLAB
NodeFile <- function(atlas, community = NULL, nodestp = NULL, nodevals = NULL, nnodes, labels, filename,outputdir){
  #nodestp provides vector of node numbers that need to be plotted
  #community can be set to community membership, whether higher in BPD/ controls, or anything else  
  #labels set equal to either 1 (full anatomical label) or 2 (node number [i.e. V_1...])
  
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
  
  # browser()
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
  } else {
    if(is.null(nodevals)){
      nodevals <- rep(1, nnodes)
    } else {
      nodevals <- nodevals[which(names(nodevals) %in% nodestp)]
    }
  }
  
  mni <- subset(nf, select = c(x.mni, y.mni, z.mni))
  if(labels == 1){
    node_label <- subset(nf, select = anat_label) 
  } else if(labels == 0) {node_label <- paste0("V", row.names(nf))}
  
  # browser()
  nf.exp <- cbind(mni, commvals, nodevals, node_label)
  
  write.table(nf.exp, file = file.path(outputdir, paste0(filename, ".node.txt")), row.names = FALSE, col.names = FALSE)
  return(nf.exp)
  
}

#/Users/nth7/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/BNV_nodefiles/in_consideration

#######################################################################################################################
##tagGraph: tags igraph object with mni coordinates, brodmanns areas
tagGraph <- function(gobj, atlas, nodefile) {
  require(dplyr)
  
  if (!is.igraph(gobj)) {
    stop(sprintf("%s is not a graph object", deparse(substitute(gobj))))
  }
  
  stopifnot(inherits(atlas, "data.frame")) #atlas should be a data.frame
  
  #these are required fields
  stopifnot(all(c("x.mni", "y.mni", "z.mni", "anat_label", "name") %in% names(atlas)))
  #match merge the observed graph against the master lookup
  #this is because the number of nodes in the graph may be a subset of the master (e.g., 240 of 264 nodes)
  gobj_vertices <- data.frame(name = V(gobj)$name, number=1:length(V(gobj)), stringsAsFactors = FALSE) # $name is convention in igraph for node naming
  vmerge <- merge(gobj_vertices, atlas, by="name", all.x=TRUE) #merge atlas attributes onto available nodes
  
  #must ensure that vmerge is sorted in the same order as gobj_vertices (merge will alpha sort the name field)
  #otherwise, assignment of attributes below could be completely out of order!
  vmerge <- vmerge %>% arrange(number) %>% select(-number) #drop number to skip adding as vertex attribute
  
  #read in the MNI coords (264 x 3) where 3 is x y z in MNI space
  #currently exported hemisphere into the atlas csv file
  #hemi <- sapply(vmerge$x.mni, function(x) { if (x < 0) { "L" } else { "R" } })
  #V(gobj)$hemi <- hemi
  
  #loop over master attributes in node file and add them as graph attributes
  toadd <- names(vmerge)
  
  #factors in the data.frame can screw up the attribute assignment (become numbers)
  #explicitly type cast to character
  atlas <- atlas %>% mutate_if(is.factor, as.character)
  
  for (i in 1:length(toadd)) {
    gobj <- set_vertex_attr(gobj, toadd[i], value=vmerge[[ toadd[i] ]])
  }
  
  return(gobj)
  
  #in progress for braingraph compatibility
  # V(gobj)
  # 
  #   name <- x.mni <- y.mni <- z.mni <- NULL
  #   coords <- eval(parse(text = g$atlas))[, list(name, x.mni, y.mni, z.mni)]
  #   setkey(coords, name)
  #   es <- get.edgelist(g)
  #   browser()
  #   dists <- sqrt(rowSums((coords[es[, 2], list(x.mni, y.mni, 
  #                                               z.mni)] - coords[es[, 1], list(x.mni, y.mni, z.mni)])^2))
  # }
}

#######################################################################################################################
##plotMetricQuant:  gives high nodal metric values to specifications for export

plotMetricQuant <- function(obj, q, metric, atlas) {
  
  v <- obj[[metric]]
  
  quant.val <- quantile(v, probs = q)
  high.val <- v[v >= quant.val]
  names(high.val) <- as.numeric(substring(names(high.val), 2, nchar(names(high.val))))
  atlas.quant <- atlas[match(names(high.val), atlas$vname),]
  atlas.quant$metric <- unname(high.val)
  #if (any(is.na(atlas.quant$vname))) { browser() }
  
  # high.val.names <- V(obj$cutgraph)$anat_label[match(names(high.val), V(obj$cutgraph)$name)]
  # highvallabeled <- high.val
  # names(highvallabeled) <- high.val.names
  return(atlas.quant)
}   

#simple function to apply density thresholding to a weighted graph
density_threshold <- function(g, d) {
  stopifnot(is_igraph(g))
  stopifnot(is.numeric(d) && d <= 1.0)
  #Obtains desired density given graph diameter
  weights <- sort(E(g)$weight, decreasing=TRUE)
  threshold <- weights[length(V(g))*(length(V(g))-1)/2 * d]
  gthresh <- delete.edges(g, which(E(g)$weight < threshold))
  gthresh <- remove.edge.attribute(gthresh, "weight")
  gthresh$density <- d #copy density into object for tracking
  return(gthresh)  
}


#small helper function to just load the nodal metrics data.frame (see compute_nodal_metrics for code that creates this structure)
load_nodal_metrics_df <- function() {
  expectFile <- file.path(basedir, "cache", paste0("dthreshnodemetrics_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
  if (file.exists(expectFile)) {
    message("Loading density-thresholded nodal statistics from file: ", expectFile)
    load(expectFile)
    return(allmetrics.nodal.df)
  } else {
    warning("Cannot find file: ", expectFile, ". You should run rs_initialize_graphs.R for this pipeline.")
    return(NULL)
  }
}

#small helper function to just load the nodal metrics data.frame (see compute_nodal_metrics for code that creates this structure)
load_allmats <- function() {
  expectFile <- file.path(basedir, "cache", paste0("adjmats_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
  if (file.exists(expectFile)) {
    message("Loading raw adjacency matrices from file: ", expectFile)
    load(expectFile)
    return(allmats)
  } else {
    warning("Cannot find file: ", expectFile, ". You should run rs_initialize_graphs.R for this pipeline.")
    return(NULL)
  }
}

nf <- function(atlas, membership, fname, nodestat=NULL, toplot=NULL, bycomm=FALSE, outputdir=getwd(), savebnv=FALSE, 
    bnvconfig=file.path(basedir, "BNV_nodefiles", "bnv_config.mat"),
    matlab="/Applications/MATLAB_R2016a.app/bin/matlab")
{
  #require(R.matlab)
  if (is.null(toplot)) { toplot <- atlas$name } #all nodes
  if (inherits(x=membership, what="communities")) {
    membership <- membership$membership #igraph communities object membership is a vector inside the object
  }
  df <- atlas
  df$membership <- membership
  df <- subset(df, name %in% toplot)
  if (is.null(nodestat)) { df$nodestat <- 1 } #replicate a dummy column
  
  df <- df[,c("x.mni", "y.mni", "z.mni", "membership", "nodestat", "name")]
  outfiles <- c()
  if (bycomm) {
    df <- split(df, df$membership)
    outfiles <- sapply(df, function(thiscomm) {
          fn <- file.path(outputdir, sub("\\.node(\\.txt)*", paste0("_comm", thiscomm$membership[1], ".node\\1"), fname, perl=TRUE))  
          write.table(file=fn, thiscomm, row.names=FALSE, col.names=FALSE)
          return(fn)
        })
  } else {
    fn <- file.path(outputdir, fname)
    write.table(file=, df, row.names=FALSE, col.names=FALSE)
    outfiles <- fn
  }
  
  if (savebnv) {
    #options(matlab="/Applications/MATLAB_R2016a.app/bin/matlab")
    #Matlab$startServer()
    #matlab <- Matlab()
    #isOpen <- open(matlab)
    #if (!isOpen) { throw("MATLAB server is not running: waited 30 seconds.") }
    
    allbnv <- sapply(outfiles, function(of) {
          bnvstring <- paste0("BrainNet_MapCfg('BrainMesh_ICBM152.nv',", #surface
              "'", bnvconfig, "',", #configuration for bnv display
              "'", of, "',", #node file
              "'", sub("\\.node(\\.txt)*", ".jpg", of, perl=TRUE), "');") #jpg file
          return(bnvstring)
        })
    allbnv <- paste(paste(allbnv, collapse=" "), "close all;", "exit;") #close all figures and exit after plotting
    
    #evaluate(matlab, bnvstring) #using R.matlab just stalls indefinitely here
    
    #cat("Calling BNV with string: ", allbnv, "\n")
    system(paste0(matlab, " -nosplash -r \"", allbnv, "\""))
    
    #close(matlab)    
  }
  return(invisible(NULL))
}
