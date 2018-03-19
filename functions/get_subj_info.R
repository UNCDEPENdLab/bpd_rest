get_subj_info <- function(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache=TRUE) {
  ################################################################################
  #################read in demographic info and combine with proper files in correct directory
  
  stopifnot(file.exists(file.path(basedir, "cache")))
  expectFile <- file.path(basedir, "cache", paste0("subjinfo_", parcellation, "_", preproc_pipeline, "_", conn_method, ".RData"))
  if (file.exists(expectFile) && allowCache==TRUE) {
    message("Loading subject info from file: ", expectFile)
    load(expectFile)
  } else {
    SPECC_rest <- read.csv(file.path(basedir, "data", "SPECC_Participant_Info.csv"), header = TRUE, stringsAsFactors = FALSE) 
    SPECC_rest$file <- NA_character_
    SPECC_rest$Luna_ID <- as.character(SPECC_rest$Luna_ID) #force to character type to align with SPECC_ID
    SPECC_rest$SPECC_ID <- as.character(SPECC_rest$SPECC_ID) #force to character type to align with Luna_ID
    SPECC_rest$ScanDate <- as.Date(SPECC_rest$ScanDate, format="%m/%d/%y")
    SPECC_rest$DOB <- as.Date(SPECC_rest$DOB, format="%m/%d/%y")
    SPECC_rest <- subset(SPECC_rest, HasRest==1 & FMRI_Exclude==0)  
    SPECC_rest <- SPECC_rest[!SPECC_rest$SPECC_ID == "023DS",] #remove excluded participant with NA field for ScanDate
    
    if(conn_method != "dens.clime_partial"){
      #expects a base directory (e.g., "adjmats") with subdirectories of <parcellation>_<conn_method> (e.g., power269_pearson)
      adjexpect <- file.path(adjmats_base, paste(parcellation, preproc_pipeline, conn_method, sep="_"))
      if (!file.exists(adjexpect)) { stop("Cannot find expected directory: ", adjexpect) } #check existence of expected directory
      
      #populate file field of subj info and verify file existence
      #figure out all proper IDs and date formats up front (Luna or SPECC)
      expectid <- with(SPECC_rest, ifelse(LunaMRI == 1, Luna_ID, SPECC_ID))
      expectdate <- with(SPECC_rest, ifelse(LunaMRI == 1, format(ScanDate, "%Y%m%d"), format(ScanDate, "%d%b%Y")))
      
      #files should be named <ID>_<DATE>_<PARCELLATION>_<CONN_METHOD><FILE_EXTENSION>
      for (i in 1:nrow(SPECC_rest)) {
        fname <- file.path(adjexpect, paste0(tolower(paste0(expectid[i], "_", expectdate[i], "_", parcellation, "_", conn_method)), file_extension))
        if (!file.exists(fname)) { 
          message("Cannot find expected adjacency matrix: ", fname) 
        } else {
          SPECC_rest$file[i] <- fname
        }    
      }
    }
    browser()
    
    #scrub subjects with lots of movement
    if (fd.scrub == TRUE) {
      if(!dir.exists("/mnt/ics/SPECC")) { message("ICS folder not mounted, unable to read FD.txt files") }
      #NOTE: dir command leads to ICS directory, which needs to be mounted on your computer    
      SPECC_rest <- filter_movement(SPECC_rest, "/mnt/ics", 0.5, .20, 10) #0.5mm FD threshold, 20% of volumes at that threshold, or any 10mm+ movement
    }
    
    save(file=expectFile, SPECC_rest)
  }
  
  
  return(SPECC_rest)
}
