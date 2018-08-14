
# Script to examine voxel missingness across subjects and ROIs, then decide what to do about it

library("brainGraph")
library(tidyverse)
setwd("/gpfs/group/mnh5174/default/Michael/bpd_rest"); basedir <- getwd()
groupdir <- "/gpfs/group/mnh5174/default"
mr_subdir <- "mni_aroma_minimal_fsl/rest1"
expect_mr_file <- "rnawuktm_rest1.nii.gz"


#worker function to get statistics for one file
fd_stats <- function(fd_file, fd_spike=0.5, max_prop_spikes=.1, max_spike=5) {
  stopifnot(file.exists(fd_file))
  fd <- read.table(fd_file)$V1 #vector of FDs
  n_spikes=sum(fd >= fd_spike)
  p_spikes <- n_spikes/length(fd)
  bad_fd <- p_spikes > max_prop_spikes || any(fd > max_spike)
  
  ret <- list(mean_fd=mean(fd), max_fd=max(fd), nvol=length(fd),
    prop_spikes=p_spikes, bad_fd=bad_fd)

  return(ret)
}

#read in and process data
specc_info <- read.csv(file.path(basedir, "data", "SPECC_Participant_Info.csv"), stringsAsFactors=FALSE) %>%
  filter(HasRest==1 & FMRI_Exclude==0) %>%
  mutate(ScanDate = as.Date(ScanDate, format="%m/%d/%y"), DOB = as.Date(DOB, format="%m/%d/%y"), Luna_ID=as.character(Luna_ID)) %>%
  rowwise() %>% mutate(mr_dir=ifelse(LunaMRI==1,
    paste0(groupdir, "/MMClock/MR_Proc/", Luna_ID, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d")), #convert to Date, then reformat YYYYMMDD
    paste0(groupdir, "/SPECC/MR_Proc/", tolower(SPECC_ID), "_", tolower(format((as.Date(ScanDate, format="%Y-%m-%d")), "%d%b%Y")))),
    mr_file=file.path(mr_dir, mr_subdir, expect_mr_file), fd_file=file.path(mr_dir, mr_subdir, "motion_info", "fd.txt"), mr_exists=file.exists(mr_file), file_id=ifelse(LunaMRI==1, Luna_ID, SPECC_ID)) %>%
  ungroup() %>%
  bind_cols(map_dfr(.$fd_file, fd_stats, fd_spike=0.5, max_prop_spikes=.20, max_spike=10)) %>%   # %>% unnest(fd_info) #nesting isn't really effective for list returns...
  filter(!bad_fd)

#to get fd data.frame alone
#vv <- map_dfr(specc_info$fd_file, fd_stats, fd_spike=0.5, max_prop_spikes=.20, max_spike=10)


# creating binary 95% mask ------------------------------------------------

#find masks
specc_info$mask <- sub(expect_mr_file, "subject_mask.nii.gz", specc_info$mr_file, fixed=TRUE)
stopifnot(all(file.exists(specc_info$mask)))

fsl_workdir <- tempdir()
fslcmd <- paste0("fslmerge -t ", file.path(fsl_workdir, "mask_merge"), " ", paste(specc_info$mask, collapse=" "))
system(fslcmd)
system(paste0("fslmaths ", file.path(fsl_workdir, "mask_merge"), " -Tmean ", file.path(fsl_workdir, "mask_prop_present")))

library(oro.nifti)
full <- readNIfTI(file.path(fsl_workdir, "mask_merge.nii.gz"), reorient=FALSE)
mni_mask <- readNIfTI("/gpfs/group/mnh5174/default/lab_resources/standard/fsl_mni152/MNI152_T1_2.3mm_brain_mask.nii", reorient=FALSE)
mask_log <- as.logical(mni_mask) #for selecting voxels in mask

p_present <- apply(full, c(1,2,3), mean) #mean over subjects will give proportion present (since all voxels are 0/1)
p_present <- p_present * mni_mask #apply MNI mask
#p_present[which(p_present == 0)] <- NA

df_present <- reshape2::melt(p_present, varnames=c("x", "y", "z"), value.name="p_present")
miss_info <- df_present %>% filter(p_present < 1 & p_present > 0)

mask_display <- full
mask_display@.Data <- p_present #depict things < 1 and > 0
mask_display@dim_[1] <- 3
mask_display@dim_[5] <- 1
mask_display@datatype <- 16 #switch to floating point, not byte
mask_display@bitpix <- 32 #switch to floating point, not byte
mask_display[mask_display==1] <- 0
summary(p_present[p_present > 0 & p_present < 1])

#write an image where the proportion of voxels is more than 0, but less than 1 (i.e., where there is between-subject heterogeneity)
#this is just used for visual diagnostics
writeNIfTI(mask_display, file.path(basedir, "Parcellation", "SPECC_n84_prop_present_less_than_1"))

#proportions identical using -Tmean from fslmaths or custom calculation above
#from_fsl <- readNIfTI(file.path(fsl_workdir, "mask_prop_present.nii.gz"), reorient=FALSE)
#all.equal(from_fsl@.Data, p_present)

#look at missingness by subject
specc_info$miss_by_subj <- apply(full, 4, function(subj) {
  sum(mni_mask[!(mni_mask & subj)]) #number of voxels from MNI mask that are not present in subject
})
hist(specc_info$miss_by_subj)

filter(specc_info, miss_by_subj > 2000) %>% arrange(desc(miss_by_subj)) %>% select(NUM_ID, SPECC_ID, Luna_ID, Wrong_fmap_dims, mean_fd, miss_by_subj)

#QA check: correlate mean functionals in MNI space before voxelwise intensity normalization
specc_info$mean_img <- sub(expect_mr_file, "awuktm_mean_float.nii.gz", specc_info$mr_file, fixed=TRUE)

fslcmd <- paste0("fslmerge -t ", file.path(fsl_workdir, "mean_merge"), " ", paste(specc_info$mean_img, collapse=" "))
system(fslcmd)

mean_imgs <- readNIfTI(file.path(fsl_workdir, "mean_merge.nii.gz"), reorient=FALSE)

similarity <- apply(mean_imgs, 4, function(mat) {
  vv <- as.vector(mat[mask_log]) #vector of voxels in MNI mask
  vv <- psych::winsor(vv, trim=0.05) #cut down on odd/extreme intensity values
  return(vv)
})

subj_sims <- cor(similarity)
diag(subj_sims) <- NA
subj_agg <- rowMeans(subj_sims, na.rm=TRUE)
hist(subj_agg)
specc_info$mean_func_similarity <- subj_agg
filter(specc_info, mean_func_similarity < .65) #folks we may want to check, though on first glance, nothing seems troubling

##ROI Diagnostics, round 1: look at missingness in the full Schaefer 422 parcellation (slightly masked by FSL MNI)
roi_dir <- file.path(basedir, "roi_diagnostics_mnh")
alldf <- bind_rows(lapply(Sys.glob(paste0(roi_dir, "/", tolower(specc_info$file_id), "*roidiagnostics_*fsl.csv")), read_csv, col_types = cols()))

alldf$subj <- sub(".*/MR_Proc/([^\\/]+)/.*$", "\\1", alldf$dataset, perl=TRUE)
#note that masking at this stage does very little since we did not apply a group threshold to the 422. Thus, any masking
#reflects voxels in Schaefer that aren't in the tight FSL MNI mask

summaries <- alldf %>% group_by(maskval) %>% dplyr::summarize(m_masked=mean(prop_masked), m_missing=mean(prop_missing), max_missing=max(prop_missing),
  min_nvox_good=min(nvox_good), min_nvox_masked=min(nvox_observed), sd_missing=sd(prop_missing))

missdf <- dplyr::filter(summaries, m_missing > .02 | m_masked > .05) %>% arrange(desc(m_missing)) %>% print(n=Inf)
filter(summaries, max_missing > .5) #ROIs with > 50% missing in one subject

bysubj <- alldf %>% group_by(subj) %>% summarise(missing=sum(nvox_total - nvox_good)) %>% arrange(desc(missing))
hist(bysubj$missing)

#Conclusions: 119, 120, 324, 325, 326 are the only ROIs that are significantly damaged
#These are all inferior temporal regions

#ROI Diagnostics, round 2: re-run ROI_TempCorr diagnostics using a Schaefer 422 95% group mask
mask95 <- full
mask95_matrix <- array(0, dim(mask95)[1:3])
mask95_matrix[p_present >= .95] <- 1 #include only voxels where at least 95% of the group is observed/present
mask95@.Data <- mask95_matrix
mask95@dim_[1] <- 3
mask95@dim_[5] <- 1
writeNIfTI(mask95, file.path(basedir, "Parcellation", "SPECC_n84_95_groupmask"))

##heavy on duplicate code here...
#data.table version
#alldf_postmask <- rbindlist(lapply(Sys.glob(paste0(roi_dir, "/", tolower(specc_info$file_id), "*roidiagnostics_*fsl_95_groupmask.csv")), fread))

#dplyr version
alldf_postmask <- bind_rows(lapply(Sys.glob(paste0(roi_dir, "/", tolower(specc_info$file_id), "*roidiagnostics_*fsl_95_groupmask.csv")), read_csv, col_types = cols()))

alldf_postmask$subj <- sub(".*/MR_Proc/([^\\/]+)/.*$", "\\1", alldf_postmask$dataset, perl=TRUE)

summaries <- alldf_postmask %>% group_by(maskval) %>% dplyr::summarize(m_masked=mean(prop_masked), m_missing=mean(prop_missing), max_missing=max(prop_missing),
  min_nvox_good=min(nvox_good), min_nvox_masked=min(nvox_observed), roi_size=min(nvox_total), avg_size=mean(nvox_good), sd_missing=sd(prop_missing))

missdf <- dplyr::filter(summaries, m_missing > .02 | m_masked > .08) %>% arrange(desc(m_missing)) %>% print(n=Inf)
filter(summaries, max_missing > .5) #ROIs with > 50% missing in one subject after the 95% group mask. Just 119
filter(summaries, m_masked > .5) #ROIs that are > 50% masked using the group 95 approach: 119, 325, 326

bysubj <- alldf_postmask %>% group_by(subj) %>% summarise(missing=sum(nvox_total - nvox_good)) %>% arrange(desc(missing))
hist(bysubj$missing)

#conclusions: use 95% mask, drop 119 altogether
schaefer <- readNIfTI(file.path(basedir, "Parcellation", "Schaefer_422_final_jul2018_fsl.nii.gz"), reorient=FALSE)
schaefer[schaefer==119] <- 0
schaefer <- schaefer * mask95
writeNIfTI(schaefer, file.path(basedir, "Parcellation", "Schaefer_421_final_jul2018_fsl_95_group_mask"))
