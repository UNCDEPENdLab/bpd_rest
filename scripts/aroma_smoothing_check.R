
# Comparision of smoothed and unsmoothed aroma ----------------------------


setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
#initialize the graph analysis pipeline (includes sourcing pipeline inputs, creating graphs, thresholding, binarization, community assignment, and calculation and reduction of global and nodal graph metrics)

source("scripts/setup_globals.R")
#RIDGE
inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE) #leave blank for defaults, consult with function for deviations from default
for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

source("scripts/estimate_euclidean_distance.R")
# get subject info for old and new ----------------------------------------

SPECC_rest <- read.csv(file.path(basedir, "data", "SPECC_Participant_Info.csv"), header = TRUE, stringsAsFactors = FALSE) 
SPECC_rest$file <- NA_character_
SPECC_rest$Luna_ID <- as.character(SPECC_rest$Luna_ID) #force to character type to align with SPECC_ID
SPECC_rest$SPECC_ID <- as.character(SPECC_rest$SPECC_ID) #force to character type to align with Luna_ID
SPECC_rest$ScanDate <- as.Date(SPECC_rest$ScanDate, format="%m/%d/%y")
SPECC_rest$DOB <- as.Date(SPECC_rest$DOB, format="%m/%d/%y")
SPECC_rest <- subset(SPECC_rest, HasRest==1 & FMRI_Exclude==0)  
SPECC_rest <- SPECC_rest[!SPECC_rest$SPECC_ID == "023DS",] #remove excluded participant with NA field for ScanDate

SPECC_rest_bp_withsmooth <- SPECC_rest
SPECC_rest_bp_nosmooth <- SPECC_rest

file_extension_withsmooth <- ".txt"
file_extension_nosmooth <- ".txt.gz"
#expects a base directory (e.g., "adjmats") with subdirectories of <parcellation>_<conn_method> (e.g., power269_pearson)
adjexpect_bp_withsmooth <- file.path(adjmats_base, paste(parcellation, preproc_pipeline, conn_method, "bp", sep="_"))
adjexpect_bp_nosmooth <- file.path(adjmats_base, paste(parcellation, preproc_pipeline, conn_method, "unsmoothed_aroma_bp", sep="_"))
if (!file.exists(adjexpect_bp_withsmooth) | !file.exists(adjexpect_bp_nosmooth)) { stop("Cannot find expected directory: ", adjexpect_bp_withsmooth, " ", adjexpect_bp_nosmooth)}  #check existence of expected directory

#populate file field of subj info and verify file existence
#figure out all proper IDs and date formats up front (Luna or SPECC)
expectid <- with(SPECC_rest, ifelse(LunaMRI == 1, Luna_ID, SPECC_ID))
expectdate <- with(SPECC_rest, ifelse(LunaMRI == 1, format(ScanDate, "%Y%m%d"), format(ScanDate, "%d%b%Y")))


# highpass adjmats ---------

for (i in 1:nrow(SPECC_rest)) {
  fname_bp_nosmooth <- file.path(adjexpect_bp_nosmooth, paste0(tolower(paste0(expectid[i], "_", expectdate[i], "_", parcellation, "_", conn_method)), file_extension_nosmooth))
  if (!file.exists(fname_bp_nosmooth)) { 
    message("Cannot find expected adjacency matrix: ", fname_bp_nosmooth) 
  } else {
    SPECC_rest_bp_nosmooth$file[i] <- fname_bp_nosmooth
  }    
}

# bandpass adjmats ---------

for (i in 1:nrow(SPECC_rest)) {
  fname_bp_withsmooth <- file.path(adjexpect_bp_withsmooth, paste0(tolower(paste0(expectid[i], "_", expectdate[i], "_", parcellation, "_", conn_method)), file_extension_withsmooth))
  if (!file.exists(fname_bp_nosmooth)) { 
    message("Cannot find expected adjacency matrix: ", fname_bp_withsmooth) 
  } else {
    SPECC_rest_bp_withsmooth$file[i] <- fname_bp_withsmooth
  }    
}



#scrub subjects with lots of movement
if(!dir.exists("/mnt/ics/SPECC")) { message("ICS folder not mounted, unable to read FD.txt files") } else {
  #NOTE: dir command leads to ICS directory, which needs to be mounted on your computer    
  SPECC_rest_bp_withsmooth <- filter_movement(SPECC_rest_bp_withsmooth, "/mnt/ics", 0.5, .20, 10) #0.5mm FD threshold, 20% of volumes at that threshold, or any 10mm+ movement
  SPECC_rest_bp_nosmooth <- filter_movement(SPECC_rest_bp_nosmooth, "/mnt/ics", 0.5, .20, 10) #0.5mm FD threshold, 20% of volumes at that threshold, or any 10mm+ movement
}


# create new + old adjmats ------------------------------------------------------

adjmats_bp_withsmooth <- import_adj_mats(SPECC_rest_bp_withsmooth, allowCache = FALSE, rmShort = rmShort)

adjmats_bp_nosmooth <- import_adj_mats(SPECC_rest_bp_nosmooth, allowCache = FALSE, rmShort = rmShort)

# New + old correlation ---------------------------------------------------

x <- array(NA, 83)
for (i in 1:83){
  subj_old <- adjmats_bp_withsmooth[i,,]
  subj_new <- adjmats_bp_nosmooth[i,,]
  
  x[i] <- cor(subj_new[lower.tri(subj_new)], subj_old[lower.tri(subj_old)])
}
x
mean(x)