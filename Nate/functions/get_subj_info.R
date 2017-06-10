get_subj_info <- function(adjmats_base, parcellation, pipeline, file_extension=".txt.gz") {
  ################################################################################
  #################read in demographic info and combine with proper files in correct directory
  SPECC_rest <- read.csv("SPECC_Participant_Info.csv", header = TRUE, stringsAsFactors = FALSE)
  SPECC_rest$file <- NA_character_
  SPECC_rest$Luna_ID <- as.character(SPECC_rest$Luna_ID) #force to character type to align with SPECC_ID
  SPECC_rest$SPECC_ID <- as.character(SPECC_rest$SPECC_ID) #force to character type to align with Luna_ID
  SPECC_rest$ScanDate <- as.Date(SPECC_rest$ScanDate, format="%m/%d/%y")
  SPECC_rest$DOB <- as.Date(SPECC_rest$DOB, format="%m/%d/%y")
  SPECC_rest <- subset(SPECC_rest, HasRest==1 & FMRI_Exclude==0)
  
  #expects a base directory (e.g., "adjmats") with subdirectories of <parcellation>_<pipeline> (e.g., power269_pearson)
  adjexpect <- file.path(adjmats_base, paste(parcellation, pipeline, sep="_"))
  stopifnot(file.exists(adjexpect)) #check existence of expected directory
  
  #populate file field of subj info and verify file existence
  #figure out all proper IDs and date formats up front (Luna or SPECC)
  expectid <- with(SPECC_rest, ifelse(LunaMRI == 1, Luna_ID, SPECC_ID))
  expectdate <- with(SPECC_rest, ifelse(LunaMRI == 1, format(ScanDate, "%Y%m%d"), format(ScanDate, "%d%b%Y")))
  
  #files should be named <ID>_<DATE>_<PARCELLATION>_<PIPELINE><FILE_EXTENSION>
  for (i in 1:nrow(SPECC_rest)) {
    fname <- file.path(adjexpect, tolower(paste0(expectid[i], "_", expectdate[i], "_", parcellation, "_", pipeline, file_extension)))
    if (!file.exists(fname)) { 
      stop("Cannot find expected adjacency matrix: ", fname) 
    } else {
      SPECC_rest$file[i] <- fname
    }
    
  }
  # ################################################################################
  # ####Framewise displacement
  # #####filter subjects with over .20 brain volumes displaced .5mm or more
  # 
  # ##In progress (make sure ics is mounted): get motion info (notes on how to implement this in RS notes folder in OneNote)
  # ####this should include mean FD and max FD at the very least, standard script removes subjects with proportion of FD >.5mm of 20% or more 
  # #####currently no safeguard against very large head movements
  # 
  # ##standard FD script
  # SPECC_rest <- filter(SPECC_rest, pr_over5mm <= .15)
  # # SPECC_rest <- filter(SPECC_rest, pr_over5mm <= .2)
  # table(SPECC_rest[,c(3,5)])
  # 
  # describe(SPECC_rest[,c(1:6, 8)])
  
  return(SPECC_rest)
}