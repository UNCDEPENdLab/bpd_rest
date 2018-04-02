# initialize brainGraph pipeline ------------------------------------------

# inputs ------------------------------------------------------------------

setwd("~/Box Sync/bpd_rest"); basedir <- getwd()
source("scripts/setup_globals.R")

inputs <- specify_inputs(
  nnodes=421,
  parcellation = "schaefer421",
  roiFile = file.path(basedir,  "data", "schaefer421_roiMat.txt"),
  thresh_weighted = "binary",
  #conn_method = "ridge.net_partial",
  conn_method="cor.shrink",
  fc_out_rm = FALSE,
  #thresh="fc",
  #thresh="prop",
  #preproc_pipeline = "nosmooth_aroma_bp_nonaggr",
  preproc_pipeline = "nosmooth_aroma_hp",
  reducemetrics =  c("degree", "page.rank", "part.coeff", "eigen.cent", "gateway.coeff.btw", "gateway.coeff.degree", "within.module.deg"),
  #rs_desired_log = logspace(log10(.01), log10(.024), 20))
  #rs_desired_log = logspace(log10(.51), log10(.68), 20) #for cor.shrink on bp data
  rs_desired_log = logspace(log10(.44), log10(.58), 20) #for cor.shrink on hp data
)

for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

add_tag <- "update2018"
file_tag <- paste0(file_tag, "_", add_tag)
file_tag_nothresh <- paste0(file_tag_nothresh, "_", add_tag)

source("scripts/estimate_euclidean_distance.R") ##creates rmShort which will delete edges close in euclidean distance


# load data and setup graphs ---------------------------------------------------------

#covars <- read.csv(paste0(basedir, "/data/SPECC_Participant_Info.csv"))
subj_info <- get(load(paste0(basedir, "/cache/subjinfo_schaefer421_nosmooth_aroma_hp_cor.shrink_fc_binary_all_update2018.RData"))) 

covars <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan, Female, file)
colnames(covars) <- c("Study.ID", "Group", "Age", "Sex", "adjmat")
covars <- covars %>% mutate(Group=factor(Group, levels=0:1, labels=c("Control", "BPD")), Sex = factor(Sex, levels = 0:1, labels = c("Male", "Female")))

sch421 <- read.csv(paste0(basedir, "/data/schaefer421_masterlookup.csv"))
sch421$index <- seq(1,length(sch421$vname),1); sch421$lobe <- sch421$X7_networks
# str(sch421)
colnames(sch421) <- c("name", "x.mni", "y.mni", "z.mni", "name.full", "hemi", "network", "index", "lobe") #to play well with any future brainGraph weirdness. WARNING right now lobe =network 
sch421 <- data.table(sch421)

#3d array of adjmay
allmats <- import_adj_mats(subj_info, rmShort = rmShort, allowCache=TRUE)

gobjs <- setup_graphs(allmats, file_tag = file_tag, file_tag_nothresh = file_tag_nothresh, fc_out_rm = fc_out_rm, allowCache=TRUE)
#pull vars into environment
allg <- gobjs$allg; allg_noneg <- gobjs$allg_noneg; allg_density <- gobjs$allg_density; agg.g <- gobjs$agg.g; allg_density_fc <- gobjs$allg_density_fc
agg.g.controls <- gobjs$agg.g.controls; agg.g.bpd <- gobjs$agg.g.bpd

#str(allg_density_fc)
yeo7 <- yeo7_community(agg.g)
allg_density_fc <- assign_communities(allg_density_fc, yeo7, "comm")


# convert to something brainGraph likes -----------------------------------

allg_density_bg <- format_to_bg(allg_density_fc, covars, rs_desired_log, allowCache = TRUE)


# MTPC nodal group differences --------------------------------------------

con.mat <- matrix(rbind(c(0, 1)), nrow=1, dimnames=list(conname=c("Group"), NULL))

mtpcVars <- data.table(level = rep("vertex", 5), #rep 5 bc five graph metrics
                       outcome = c("degree", "PC", "ev.cent", "btwn.cent", "within.mod"),
                       N = rep(500, 5))
covars_test <- data.table(covars)[,c(-3,-4,-5)]
#design_matrix_group <- model.matrix(~Group , covars_test)


mtpc_nodal_group <- list()
for(metric in mtpcVars$outcome){
  mtpc_nodal_group[[metric]] <- mtpc(g.list = allg_density_bg, 
                           thresholds = rs_desired_log, covars = covars_test, 
                           measure = metric,
                           con.mat = con.mat,
                           con.type = 't', level = 'vertex',
                           X = design_matrix_group)

}

sig.nodes <- rbindlist(lapply(mtpc_nodal_group, function(x){
  x$DT[A.mtpc > A.crit, .SD[1], by = region]
}))

mtpc_nodal_group$degree$DT
