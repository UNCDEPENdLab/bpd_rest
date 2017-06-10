###################################################################################
##Graph Utility Functions: NodeFile(), tagGraph(), plotMetricQuant(), getMotionInfo()
##MotionInfoSPEC: function pulls motion info from specified dir and exports a list containing: 
#directories not found, fd.txt files not found, raw FD data, and a table with relevant FD information based on specifications 

MotionInfoSPECC <- function(dir, thresh, spikes.exclude, max.exclude){
  ##thresh indicates what specifies a significant head movement
  ##spikes.exclude specifies the cutoff for the proportion of movements above the threshold at which to exclude a given subject
  ##max.exclude is similar to spikes.exclude, but sets a cutoff for large head movements
  require(plyr); require(dplyr); require(tidyr)
  options(dplyr.width=200)
  return_list <- list()
  idfile <- paste0(dir, "/SPECC/SPECC_Participant_Info.csv")
  idinfo <- read.csv(idfile)
  ## MH has now converted all SPECC MR directory names to all lower case to allow for match on case-sensitive filesystem
  ## and to make the naming consistent
  
  ##pull subject folders
  idinfo <- idinfo %>% select(-Notes) %>% rowwise() %>% mutate(mr_dir=ifelse(LunaMRI==1,
                                                                             paste0(dir, "/MMClock/MR_Proc/", Luna_ID, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d")), #convert to Date, then reformat YYYYMMDD
                                                                             paste0(dir, "/SPECC/MR_Proc/", tolower(SPECC_ID), "_", tolower(format((as.Date(ScanDate, format="%Y-%m-%d")), "%d%b%Y")))))
  ##pull fd info into dataframe
  fd.info.raw <- data.frame()
  fd.info <- data.frame()
  fd.notfound <- data.frame()
  for(id in idinfo$NUM_ID){
    thissub <- idinfo[which(idinfo$NUM_ID == id),]
    if(file.exists(paste0(thissub$mr_dir, "/mni_5mm_wavelet/rest1/motion_info/fd.txt"))){
      #this.sub.fd <- t(data.frame(read.table(paste0(thissub$mr_dir, "/mni_5mm_wavelet/rest1/motion_info/fd.txt"))))
      #rownames(this.sub.fd) <- as.character(thissub$SPECC_ID)
      this.sub.fd <- data.frame(as.character(thissub$SPECC_ID), read.table(paste0(thissub$mr_dir, "/mni_5mm_wavelet/rest1/motion_info/fd.txt")))
      colnames(this.sub.fd) <- c("Subj", "FD")
      # browser()
      
      meanFD = mean(this.sub.fd$FD)
      maxFD = max(this.sub.fd$FD)
      maxExclude = if(maxFD > max.exclude) 1 else 0
      prop.thresh = sum(this.sub.fd$FD > thresh)/length(this.sub.fd$FD)
      thresh.exclude = if(prop.thresh > spikes.exclude) 1 else 0
      num_volumes = length(this.sub.fd$FD)
      
      thissub.desc <- data.frame(Subj = as.character(thissub$SPECC_ID),meanFD, maxFD, maxExclude, prop.thresh, thresh.exclude, num_volumes)
                          
      if(!(num_volumes==300)){cat("Subject ", as.character(thissub$SPECC_ID), " does not have the standard number of volumes", "\n")}
      fd.info <- rbind(fd.info, thissub.desc)
      fd.info.raw <- rbind(fd.info.raw, this.sub.fd)
    } else{
      fd.notfound <- rbind(fd.notfound, data.frame(thissub))
      cat("FD.txt file not found for subject: ", as.character(thissub$SPECC_ID), "\n")
    }
  }
  #verify that mr_dir is present as expected, if not, store in dir.notfound
  idinfo$dirfound <- file.exists(idinfo$mr_dir)
  dir.notfound <- subset(idinfo, dirfound==FALSE)
  cat("MR_dir not found:", as.character(dir.notfound$SPECC_ID))
  
  return_list[["dir.notfound"]] <- data.frame(dir.notfound)
  return_list[["fd.txt.notfound"]] <- fd.notfound
  return_list[["fd.info.raw"]] <- fd.info.raw
  return_list[["fd.info"]] <- fd.info
  

  return(return_list)
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
  vnames <- data.frame(atlas[,"vname"])
  names(vnames) = "vname"
  stopifnot(all(c("x.mni", "y.mni", "z.mni", "anat_label") %in% names(nf)))
  row.names(nf) <- as.numeric(vnames$vname)
  
  if(!is.null(community)){
    commvals <- community
    names(commvals) <- vnames$vname
  } else {
    commvals = rep(1,length(vnames[,1]))
    names(commvals) <- vnames$vname
  }
  
  if(!is.null(nodestp)) {
    nf <- nf[which(nf$vname %in% nodestp),]
    commvals <- commvals[which(as.numeric(names(commvals)) %in% nodestp)]
    # nf$anat_label <- paste0("V_", nf$vname)
    # nf <- nf %>% select(-vname)
    names(nodevals) <- vnames$vname
    
    if(is.null(nodevals)){
      nodevals <- rep(1, length(nodestp))
    } else {
      nodevals <- nodevals[which(as.numeric(names(nodevals)) %in% nodestp)]
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



#######################################################################################################################
##tagGraph: tags igraph object with mni coordinates, brodmanns areas
      
tagGraph <- function(gobj, mergeby="vname", nodefile, atlasname) {
  if (!is.igraph(gobj)) {
    stop(sprintf("%s is not a graph object", deparse(substitute(gobj))))
  }
  
  f <- read.csv(nodefile) #this should be a csv file containing minimally MNI x, y, z, and Label. Could also contain a given community structure (e.g., Power's)
  
  #these are required
  stopifnot(all(c("x.mni", "y.mni", "z.mni", "anat_label", mergeby) %in% names(f)))
  
  #match merge the observed graph against the master lookup
  #this is because the number of nodes in the graph may be a subset of the master (e.g., 240 of 264 nodes)
  vnames <- data.frame(v = V(gobj)$name)
  names(vnames) <- mergeby
  f <- merge(f, mergeby, all.x=TRUE)
  
  #read in the MNI coords (264 x 3) where 3 is x y z in MNI space
  hemi <- sapply(f$x.mni, function(x) { if (x < 0) { "L" } else { "R" } })
  V(gobj)$hemi <- hemi
  
  #loop over master attributes in node file and add them as graph attributes
  toadd <- names(f)
  
  for (i in 1:length(toadd)) {
    gobj <- set_vertex_attr(gobj, toadd[i], value=as.character(f[[toadd[i] ]]))
    #V(gobj)[[toadd[i]]] <- f[[toadd[i] ]]
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


