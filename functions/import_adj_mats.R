import_adj_mats <- function(subj_info, allowCache=TRUE, rmShort=NULL, densities_desired = NULL) {
  #NB. several variables are (implicitly) pulled from the global environment here:
  #   nnodes, parcellation, conn_method, basedir
  
  if(is.null(densities_desired)){densities_desired <- seq(.05,.25,.01)} #create densities desired list if nothing is passed in

  # browser()
  stopifnot(file.exists(file.path(basedir, "cache")))
  
  expectFile <- file.path(basedir, "cache", paste0("adjmats_", file_tag, ".RData"))
  
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading raw adjacency matrices from file: ", expectFile)
    allmats <- get(load(expectFile))
  } else {
    sep <- ""
    
    if(conn_method == "dens.clime_partial"){
      # adapted script for dens.clime to export a list with subjs as elements containing densityxnodexnode 3d arrays
      allmats <- get(load("adjmats/schaefer422_aroma_dens.clime_partial/dens.clime.all.RData"))
      names(allmats) <- toupper(sub("_.*", "", names(allmats)))
      allmats <- allmats[names(allmats) %in% subj_info$Luna_ID | names(allmats) %in% subj_info$SPECC_ID]
    
      allmats.list <- list()
      
      for (subj in 1:nrow(subj_info)){
        thrdim_allmats <- array(NA, c(length(densities_desired), nnodes, nnodes),
                                dimnames=list(den = densities_desired, roi1=atlas$name, roi2=atlas$name))
        for (den in 1:length(densities_desired)){
          
          thrdim_allmats[den,,] <- allmats[[subj]][[den]]
          }
        allmats.list[[subj]] <- thrdim_allmats
      }
      
      names(allmats.list) <- subj_info$SPECC_ID
    } else {
    ##preallocate empty array to read adjacency matrices into
    allmats <- array(NA, c(nrow(subj_info), nnodes, nnodes), 
                     dimnames=list(id = subj_info$SPECC_ID, roi1=atlas$name, roi2=atlas$name))

    for (f in 1:nrow(subj_info)) {
      if(conn_method == "dens.clime_partial_plateau"){
        x <- load(file = subj_info[f, "file"])        
      } else {
        cat("Reading file:", as.character(subj_info[f,"file"]), "\n")
        m <- as.matrix(read.table(as.character(subj_info[f,"file"]), sep=sep, header=TRUE))
      }
      
      allmats[f,,] <- m
      }
    }
    
    
    #remove short distance connections using 0/1 matrix passed in
    if (!is.null(rmShort)) {
      
      message("Zeroing ", sum(rmShort[lower.tri(rmShort)] == 0), " connections less than ", roi.dist, "mm from raw adjacency matrices. ")

      #apply the short-range removal by replicating the 2D rmShort to a 3D array of the same size as
      #allmats, then elementwise multiplication (since this is vectorized and fast)
      if(!exists("allmats.list")){
      rmtest <- aperm(replicate(n=dim(allmats)[1], rmShort), c(3,1,2))
      allmats_rm <- allmats * rmtest
      
      # pdf("which_edges_go.pdf", width = 11, height = 8)
      # histogram(allmats, breaks = 50)
      # hist(allmats_rm)
      # 
      } else {
        ##in the case of dens.clime, allmats list has nsubj list elements with a density x nodes x nodes 3d array. 
        #so, what needs to be done here is create 3D array of the same size and lapply over the list to perform elementwise multiplication
        rmtest <- aperm(replicate(n=length(densities_desired), rmShort), c(3,1,2))
        
        allmats.list <- lapply(allmats.list, function(sub){
          allmats.short.remove <- sub * rmtest
          return(allmats.short.remove)
        })
      }
      
     
      
      
      #slightly different tack using abind to create array (not faster)
      #system.time(rmtest2 <- do.call(abind, list(along=0, replicate(n=dim(allmats)[1], rmShort, simplify=FALSE))))
      
      #interesting approach using afill from the abind package, but not really faster than above
      #rmtest <- array(NA, dim=dim(allmats), dimnames=dimnames(allmats))
      #dimnames(rmShort) <- dimnames(allmats)[c(2,3)]
      #system.time(afill(rmtest,TRUE,,) <- rmShort)
      
      #plyr aaply keeps the dimensions intact, but is slow
      #allmats <- plyr::aaply(allmats, 1, function(m) { m <- m * rmShort; return(m) }))
      
      #basic apply over rows, apply(allmats, 1), throws away array dimensions and dimnames
    }
  }
  
  ##remove diagonal regardless of you remove short edges
  rmdiag <- array(1, c(nrow(subj_info), nnodes, nnodes), 
                   dimnames=list(id = subj_info$SPECC_ID, roi1=atlas$name, roi2=atlas$name))
  
  for(i in 1:nrow(subj_info))(
     diag(rmdiag[i,,]) <- 0
  )
  
  allmats <- allmats *rmdiag
  
  
  #return 3d array to caller
  if(exists("allmats.list")){
    save(allmats.list, file = expectFile)
    } else {save(allmats, file=expectFile) }#cache results
  return(allmats)
  }
