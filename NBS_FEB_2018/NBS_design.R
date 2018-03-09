# Generate NBS design matrices --------------------------------------------

##IN PROGRESS: IMPLEMENTATION OF NBS AND DELTACON-ATTR

##this attempts to use brainGraph code, instead, let's rework for input into MATLAB
# allmats.nbs <- aperm(allmats, c(2,3,1))
# covars.NBS <- data.table(Group = factor(subj_info$BPD))
# con.vec <- c(1,-1)
# nbs.out <- NBS(allmats.df, covars.NBS, con.vec, p.init = .001, N = 1000, symmetric = TRUE)

##another failed attempt:
# design_matrix_nbs <- subj_info$BPD
# design_matrix_nbs <- cbind(design_matrix_nbs, ifelse(design_matrix_nbs == 0, 1, 0))
# write.table(design_matrix_nbs, file = “/Users/nth7/Documents/MATLAB/NBS1.2/schaefer422_aroma_ridge.net_partial/design_matrix_group.txt”, row.names = FALSE, col.names = FALSE)

setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
#initialize the graph analysis pipeline (includes sourcing pipeline inputs, creating graphs, thresholding, binarization, community assignment, and calculation and reduction of global and nodal graph metrics)

source("scripts/setup_globals.R")

#specify pipeline inputs.leave blank for defaults, consult with function for deviations from default#RIDGE
# inputs <- specify_inputs(thresh_weighted = "binary", 
#                          fc_out_rm = FALSE, 
#                          preproc_pipeline = "nosmooth_aroma_bp_nonaggr",
#                          reducemetrics =  c( "degree", "page.rank", "part.coeff", "eigen.cent", "gateway.coeff.btw", "gateway.coeff.degree", "within.module.deg"),
#                          rs_desired_log = logspace(log10(.01), log10(.02), 20)) 
#PEARSON
inputs <- specify_inputs(thresh_weighted = "binary", 
                         conn_method = "pearson",
                         fc_out_rm = FALSE, 
                         preproc_pipeline = "nosmooth_aroma_bp_nonaggr",
                         reducemetrics =  c( "degree", "page.rank", "part.coeff", "eigen.cent", "gateway.coeff.btw", "gateway.coeff.degree", "within.module.deg"),
                         rs_desired_log = logspace(log10(.01), log10(.02), 20)) 
for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

##if you want to add on an additional tag
# file_tag <- paste0(file_tag, "_normshort")

##load subj info
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)

#reexport files for NBS
subj_info_nbs <- arrange(subj_info, BPD, SPECC_ID)
subj_info$NBS_ID <- 1:nrow(subj_info)
subj_info <- as.data.frame(subj_info) #get rid of tibble weirdness

for (i in 1:nrow(subj_info)) {
  x <- read.table(subj_info[i,"file"])
  write.table(file=paste0("/Users/nth7/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/NBS_schaefer422_nosmooth_aroma_bp_nonaggr_pearson/matrices_numeric/subj_", sprintf("%03d", i), ".txt"), x=x, row.names=FALSE, col.names=FALSE)
}
# 
# nbs.mat.ixn <- model.matrix(lm(NBS_ID ~ BPD*AgeAtScan, subj_info))
# subj_info$BPD
# write.table(file="/Users/nth7/Documents/MATLAB/NBS1.2/schaefer422_aroma_ridge.net_partial/design_matrix_bpdxage.txt", x=nbs.mat.ixn, row.names=FALSE, col.names=FALSE)
# 
# 
# model.matrix(lm(NBS_ID~BPD, subj_info))
# 
# design_matrix_nbs <- subj_info$NBS_ID
# design_matrix_nbs <- cbind(design_matrix_nbs, ifelse(design_matrix_nbs == 0, 1, 0))
# write.table(design_matrix_nbs, file ="/Users/nth7/Documents/MATLAB/NBS1.2/schaefer422_aroma_ridge.net_partial/design_matrix_group.txt", row.names = FALSE, col.names = FALSE)
# 
# 
