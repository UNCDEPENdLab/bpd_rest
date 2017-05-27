###################################################################################
##Graph Utility Functions: NodeFile(), tagGraph(), plotMetricQuant()

##NodeFile: function outputs node file readable by BrainNet Viewer in MATLAB
NodeFile <- function(graph, nodestp = NULL, nodestat = NULL, community = NULL, atlas, outputdir=getwd(), edgestat = NULL, name = NULL){
  
  nf <- atlas
   #if want to read in a csv rather than a predfined atlas in R
   #nf <- read.csv(atlas, header = TRUE)        #####csv file should have at least mni attributes and anatomical labels
   
        vnames <- data.frame(atlas[,"vname"])
        names(vnames) = "vname"
        
       
        #match merge the observed graph against the master lookup
        #this is because the number of nodes in the graph may be a subset of the master (e.g., 240 of 264 nodes)
        # nf <- merge(nf, vnames, by = "vname")
        # nf$vname <- as.numeric(sub("V", "", nf$vname, fixed=TRUE))
        # nf <- nf[order(nf$vname),] #keep the ordered version
        # row.names(nf) <- NULL #remove confusing names
      browser()
        #these are required in atlas variable
        stopifnot(all(c("x.mni", "y.mni", "z.mni", "anat_label") %in% names(nf)))
        if(!is.null(nodestp)) {
          nodevals <- nodestat
          #nodevals <- nodestp[,"metric"]
          names(nodevals) = NULL
        
        }else {
            nodevals = rep(1,length(names(vnames)))  
            }
        
        if(!is.null(community)){commvals <- community}else {
          commvals = rep(1,length(names(vnames)))  
        }
        
        #pull mni coordinates
        mni <- subset(nf, select = c(x.mni, y.mni, z.mni))
        anat_label <- subset(nf, select = anat_label)
        nf.exp <- cbind(mni, commvals, nodevals, anat_label)
        nf.exp[,"anat_label"] <- gsub(" ", "_", nf.exp[,"anat_label"], perl=TRUE)
        nodestp[,"anat_label"] <- gsub(" ", "_", nodestp[,"anat_label"], perl=TRUE)
        
        #nf.exp <- nf.exp[c(10, 84, 85, 107, 109, 196, 230),]
        
        #cut_adj.exp <- as.matrix(as_adj(graph, attr = edgestat))
        
        ##if we have a subset of nodes we would like to plot
        
        if(!is.null(nodestp)) {
          nf.exp <- nf.exp[match(as.character(nodestp[,"anat_label"]), nf.exp[, "anat_label"]),]
          row.names(nf.exp) <- row.names(nodestp)
          }

        
        # write.table(file = file.path(outputdir), x = nf.exp, row.names = FALSE, col.names = FALSE)
        write.table(file=file.path(outputdir, paste0(name, ".node.txt")), x=nf.exp, row.names=FALSE, col.names=FALSE)
        
        
      
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

