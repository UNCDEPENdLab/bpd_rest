# Specify Pipeline Inputs send specification message and caluclate euclidean distances to fix to zero -------------------------------------------------
specify_inputs <- function(nnodes = 422, 
                           parcellation = "schaefer422", 
                           roiFile = file.path(basedir,  "data", "schaefer422_roiMat.txt"),
                           roiMat = read.table(roiFile, header=FALSE, col.names=c("x", "y", "z", "roi")), 
                           atlasFile = paste0(parcellation, "_atlas.csv"),
                           atlas = read.csv(file.path(basedir, "data", atlasFile), header = TRUE, stringsAsFactors = FALSE), 
                           use.infomap = 0, 
                           use.yeo = 1,
                           data.reduce = "fa",
                           weighted.reduce = NULL,
                           reduce.method = "a_priori",
                           node.file.dir = file.path(basedir, "Viz_files"),
                           include.weighted = 1,
                           thresh = "fc",
                           thresh_weighted = "weighted",
                           preproc_pipeline = "aroma", 
                           fc_out_rm = 1, 
                           conn_method = "ridge.net_partial",
                           adjmat_extension = ".txt.gz",
                           roi.dist = 20,
                           metricstorun.nodal = c("eigen.cent","degree", "closeness", "betweenness.node", "page.rank",  "part.coeff", "within.module.deg.zscore", "local.clustering", "gateway.coeff.btw", "gateway.coeff.degree", "between.module.deg.zscore", "within.module.deg", "between.module.deg", "leverage.cent"),
                           metricstorun.global = c("characteristic_path_length", "clustering_coefficient", "small_worldness", "modularity"),
                           reducemetrics = c("eigen.cent","degree", "betweenness.node", "page.rank",  "part.coeff", "gateway.coeff.btw", "gateway.coeff.degree", "within.module.deg"),
                           adjmats_base = file.path(basedir, "adjmats"),
                           results_dir = file.path(basedir, "results", paste0(parcellation, "_", conn_method)),
                           figure.dir = file.path(basedir, "figures", paste0(parcellation, "_", conn_method)),
                           densities_desired = seq(.05, .25, .01),
                           rs_desired_log = logspace(log10(.01), log10(.02), 20),
                           neg_edges = FALSE
){


# Specify inputs ----------------------------------------------------------
##empty list of inputs (to be extracted later)
inputs <- list()
    
inputs[["nnodes"]] <- nnodes  #varying this will require you to change the .txt coordinate file you read in for roiMat (should be 4 X nnodes txt file labelled:"x", "y", "z", "roi") and the masterlookup atlas
inputs[["parcellation"]] <- parcellation
inputs[["roiFile"]] <- roiFile
inputs[["roiMat"]] <- roiMat
inputs[["atlasFile"]] <- atlasFile
inputs[["atlas"]] <- atlas #eventually would be nice to have a lockstep name here and re-use parcellationName
inputs[["use.infomap"]] <- use.infomap #from earlier community detection efforts, has since been decommissioned
inputs[["use.yeo"]] <- use.yeo
inputs[["data.reduce"]] <- data.reduce #can be fa or pca
inputs[["weighted.reduce"]] <-weighted.reduce #run fa/pca on weighted graphs?
inputs[["reduce.method"]] <- reduce.method #can be set to "all", "a_priori", or "metrics". This specifies how data reduction step is performed across metrics and densities (for binary analyses)
inputs[["node.file.dir"]] <- node.file.dir
inputs[["include.weighted"]] <- include.weighted #0 = no 1 = yes. Shall you compute weighted (non-thresholded) graphs?
inputs[["thresh"]] <- thresh #none, fc or prop    
inputs[["thresh_weighted"]] <- thresh_weighted #"weighted" or "binary" Shall you compute thresholded graphs and retain edge weights?
inputs[["preproc_pipeline"]] <- preproc_pipeline #method for data preprocessing. This corresponds to mni_5mm_aroma data. 
inputs[["fc_out_rm"]] <- fc_out_rm #Remove 2 mean FC outliers? 1 = yes, 0 = no
inputs[["neg_edges"]] <- neg_edges

##UPDATE 8/8/17. no 5mm spatial smoothing re: Alakorkko et al 2017

inputs[["conn_method"]] <- conn_method

###Below remains for user to look at possible conn_method options from earlier iterations to test

#conn_method <- "pearson" #Jun2017: aroma preprocessing, pearson correlations
#conn_method <- "cor.shrink" #Jun2017: aroma preprocessing, shrinkage estimator of correlation
#conn_method <- "pcor.shrink_partial" #Jun2017: aroma preprocessing, shrinkage estimator of *partial* correlation
#conn_method <- "ridge.net_partial" #Jun2017: aroma preprocessing, shrinkage estimator of *partial* correlation
#conn_method <- "dens.clime_partial" #Sept2017 aroma preprocessing, density-based approach for partial correlation estimation (uses Constrained L1-Minimization (CLIME)) from Wang et al (2016)
#conn_method <- "dens.clime_partial_plateau" #Sept2017 aroma preprocessing, density-based approach for partial correlation estimation (uses Constrained L1-Minimization (CLIME)) from Wang et al (2016). based on the plateau of the lambda parameter
#conn_method <- "quic" #coming soon
# conn_method <- "pearson_fisherz" #uses older files (from wavelet 5mm on power 269)
#conn_method  <- "scotmi"  #decommissioned due to not enough time points in our data

inputs[["adjmat_extension"]] <- adjmat_extension #used to ID adjmats to import into pipeline

inputs[["roi.dist"]] <- roi.dist #distance in mm to zero connections when setting up graphs (cf. Power 2011)
inputs[["metricstorun.nodal"]] <- metricstorun.nodal
inputs[["metricstorun.global"]] <- metricstorun.global
inputs[["reducemetrics"]] <- reducemetrics

inputs[["adjmats_base"]] <- adjmats_base
inputs[["results_dir"]] <- results_dir
inputs[["figure.dir"]] <- figure.dir


inputs[["densities_desired"]] <- densities_desired #used globally for binary graphs

#for FC thresholding
inputs[["rs_desired_log"]] <- rs_desired_log

# Send Specification Message ----------------------------------------------


# message("Initializing RS graph analysis with the following settings: ")
# message("Parcellation: ", parcellation)
# message("roiMat: ", roiFile)
# message("atlas: ", atlasFile)
# message("preprocessing pipeline: ", preproc_pipeline)
# message("adj mat calculation conn_method: ", conn_method)
# message("adjmats directory: ", adjmats_base)
# message("results directory: ", results_dir)
# message("figures directory: ", figure.dir)
# message("Data reduction method: ", data.reduce, " + ", reduce.method)
# message("Nodal metrics calculated: ", paste(metricstorun.nodal, collapse = ", "))
# message("Global metrics calculated: ", paste(metricstorun.global, collapse = ", "))
# if (thresh_weighted ==0){message("Thresholding applied: ", paste0(thresh, " + binary"))} else{
#   message("Thresholding applied: ", paste0(thresh, " + weighted"))
# }

#super ugly script to give file tag
if (include.weighted == 1){
  if(fc_out_rm == 1){
    inputs[["file_tag"]] <- paste(parcellation, preproc_pipeline, conn_method, thresh, thresh_weighted, "out_rm", sep = "_")
    inputs[["file_tag_nothresh"]] <- paste(parcellation, preproc_pipeline, conn_method, "weighted_out_rm", sep = "_")
  } else {
    inputs[["file_tag"]] <- paste(parcellation, preproc_pipeline, conn_method, thresh, thresh_weighted, "all", sep = "_")
    inputs[["file_tag_nothresh"]] <- paste(parcellation, preproc_pipeline, conn_method, "weighted_all", sep = "_")}
} else {
  if(fc_out_rm == 1){
    inputs[["file_tag"]] <- paste(parcellation, preproc_pipeline, conn_method, thresh, thresh_weighted, "out_rm", sep = "_")
  } else {inputs[["file_tag"]] <- paste(parcellation, preproc_pipeline, conn_method, thresh, thresh_weighted, "all", sep = "_")}
}

if(neg_edges) {inputs[["file_tag"]] <- paste0(inputs[["file_tag"]], "_negative")}

message("File tag for current thresholded specifications: ", inputs[["file_tag"]])
if(include.weighted == 1) {message("File tag for current weighted no threshold specifications: ", inputs[["file_tag_nothresh"]])}



return(inputs)
# (nnodes, parcellation, roiFile, roiMat, atlasFile, atlas, use.infomap, use.yeo, data.reduce,
#        weighted.reduce, reduce.method, node.file.dir, include.weighted, thresh, thresh_weighted,
#        preproc_pipeline, fc_out_rm, conn_method, adjmat_extension, roi.dist, metricstorun.nodal,
#        metricstorun.global, reducemetrics, adjmats_base, results_dir, figure.dir, densities_desired,
#        rs_desired_log)
}
