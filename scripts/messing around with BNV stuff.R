
membership.file <- file.path(getwd(), "data", "membership.yeo.csv")
stopifnot(file.exists(membership.file))
community <- as.numeric(as.matrix(read.csv(membership.file)))

bnv_dir <- "/Users/nth7/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/BNV_nodefiles"
Node.File <- NodeFile(atlas = atlas_fancy, 
                      community = community, 
                      nnodes = nnodes, 
                      labels = 1, 
                      filename = "schaefer422_all", 
                      outputdir = bnv_dir)

atlas$anat_label <- gsub(" ", ".", atlas$anat_label)

#nice looking anats


atlas_fancy <- atlas
atlas_fancy$anat_label <-  gsub(" ", ".",  paste0(atlas$name, ":", atlas$anat_label))
