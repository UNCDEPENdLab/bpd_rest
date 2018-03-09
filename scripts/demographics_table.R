#####################generate subjects demographics #####################
setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
source("scripts/setup_globals.R")
inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = 0) #leave blank for defaults, consult with function for deviations from default
for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

#get_subj info, includes motion scrubbing procedure
subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)

#pull demogs vars of interest for merging
demogs_table <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan, Female)

demogs_tomerge <- data.frame(read.csv("data/Demogs_all_RS.csv", header = TRUE, na.strings = c(".", "NA"), stringsAsFactors = FALSE))
# str(demogs_tomerge)


# df_test <- demogs_tomerge %>% dplyr::select(SPECC.ID, ETHNICITY, RACE1, RACE2, RELATIONSHIP, SEXUAL.ATTRACTION, SEXUALITY, RELATIONSHIP.LENGTH, PAST.30.PHYSICAL.HEALTH, PAST.30.MENTAL.HEALTH, OVERALL.HEALTH.)
df_test <- demogs_tomerge %>% dplyr::select(-DATE, -SEX, -GROUP, -BIRTHDATE, -AGE)

# head(df_test)
# df_test[df_test == "."] <- NA

colnames(df_test)[1] <- "SPECC_ID"
df_test$SPECC_ID <- as.character(df_test$SPECC_ID)
for(i in 1:nrow(df_test)){
  df_test$SPECC_ID[i] <- gsub("_", "",df_test$SPECC_ID[i])
}

demogs_table <- left_join(demogs_table, df_test, by = "SPECC_ID")

demogs_table$RACE1 <- as.numeric(demogs_table$RACE1)
demogs_table$RACE2 <- as.numeric(demogs_table$RACE2)


RACE_ALL <- rep(NA, length(demogs_table$RACE1))
for(i in 1:length(RACE_ALL)){
  if(is.na(demogs_table$RACE1[i]) & is.na(demogs_table$RACE2[i])){
    RACE_ALL[i] <- 99
  } else if(!is.na(demogs_table$RACE1[i]) & is.na(demogs_table$RACE2[i])){
    RACE_ALL[i] <- demogs_table$RACE1[i]
  } else if(is.na(demogs_table$RACE1[i]) & !is.na(demogs_table$RACE2[i])){
    RACE_ALL[i] <- demogs_table$RACE2[i]
  } else if(!is.na(demogs_table$RACE1[i]) & !is.na(demogs_table$RACE2[i])){
    RACE_ALL[i] <- 6
  }
}

demogs_table <- data.frame(demogs_table, RACE_ALL = RACE_ALL) %>% dplyr::select(-RACE1, -RACE2)

##hacky-split
BPD_demo <- demogs_table %>% dplyr::filter(BPD == 1) 
control_demo <- demogs_table %>% dplyr::filter(BPD == 0)

#ethnicity info
sum(is.na(demogs_table$ETHNICITY))
sum(is.na(BPD_demo$ETHNICITY))
table(demogs_table$BPD, demogs_table$ETHNICITY)

#race info
sum(is.na(demogs_table$RACE_ALL))
table(demogs_table$BPD, demogs_table$RACE_ALL) 
##99 = not provided/missing
##6 = bi/multiracial

#Income info
sum(is.na(demogs_table$TOTAL.INCOME))
sum(is.na(BPD_demo$TOTAL.INCOME))
table(demogs_table$BPD, demogs_table$TOTAL.INCOME)

#Sexuality info
sum(is.na(demogs_table$SEXUALITY))
sum(is.na(BPD_demo$SEXUALITY))
table(demogs_table$BPD, demogs_table$SEXUALITY)



