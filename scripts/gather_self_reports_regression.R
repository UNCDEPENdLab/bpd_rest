# Gather self reports and score for regression analyses ---------------------------------------------
setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
source("scripts/setup_globals.R")
source("functions/meanreplace.R")
subj_info <- as.data.frame(get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)) %>% dplyr::select(SPECC_ID, Luna_ID, BPD, AgeAtScan)

colnames(subj_info) <- c("SPECC_ID", "LUNA_ID", "BPD", "Age")
basedir_surveys <- "~/Box Sync/DEPENd/Projects/SPECC/SelfReports/data/surveys/";setwd(basedir_surveys)


# UPPS scoring and imputation ---------------------------------------------

#code for UPPS scoring 
upps <- read.csv("UPPS.csv", header=TRUE)
upps$SPECC_ID <- sub("_", "", as.character(upps$SPECC_ID)) 
upps$LUNA_ID <- as.character(upps$LUNA_ID)

upps <- subj_info %>% left_join(upps, by = c( "SPECC_ID" = "SPECC_ID", "LUNA_ID" = "LUNA_ID"))

# Updated 12/3/17 NH ---

#reverse score items
reverseItems <- paste("UPPS", c(2, 7, 12, 17, 22, 29, 34, 39, 44, 50, 53, 58, 9, 47, 3, 8, 13, 18, 23, 26, 31, 36, 41, 46, 51, 56, 5, 10, 15, 20, 25, 30, 35, 40, 45, 49,52, 54, 57, 59), sep=".")
upps[,reverseItems] <- lapply(upps[,reverseItems], function(x) { 5 - x })

#mean imputation
uppsItems <- paste("UPPS", 1:59, sep=".")

for (i in 1:nrow(upps)) {
  miss <- which(is.na(upps[i, uppsItems]))
  if (length(miss) > 0) {
    cat("For subject id: ", paste(sapply(upps[i, c("LUNA_ID", "SPECC_ID")], as.character), collapse="/"), ", these items are missing: ", paste(miss, collapse=", "), "\n", sep="")
  }
}

#SPECC ID 066DW: missing item 49
upps <- meanreplace(upps, "066DW",varprefix = "UPPS.", subscaleitems = c(5,10,15,20,25,30,35,40,45,49,52,54,57,59))

#SPECC ID 074DS: missing items 4,6,7,8,11,27,49
upps <- meanreplace(upps, "074DS",varprefix = "UPPS.", subscaleitems = c(5,10,15,20,25,30,35,40,45,49,52,54,57,59))
upps <- meanreplace(upps, "074DS",varprefix = "UPPS.", subscaleitems = c(4,9,14,19,24,27,32,37,42,47))
upps <- meanreplace(upps, "074DS",varprefix = "UPPS.", subscaleitems = c(1,6,11,16,21,28,33,38,43,48,55))
upps <- meanreplace(upps, "074DS",varprefix = "UPPS.", subscaleitems = c(2,7,12,17,22,29,34,39,44,50,53,58))
upps <- meanreplace(upps, "074DS",varprefix = "UPPS.", subscaleitems = c(3,8,13,18,23,26,31,36,41,46,51,56))
#rest of missing data is completely missing from subjs

#remove one extraneous 040_BS
upps <- upps[-22,]

#define variables and score items
upps$NegUrg <- with(upps, UPPS.2 + UPPS.7 + UPPS.12 + UPPS.17 + UPPS.22 + UPPS.29 + UPPS.34 + UPPS.39 + UPPS.44 + UPPS.50 + UPPS.53 + UPPS.58)
upps$LackPrem <- with(upps, UPPS.1 + UPPS.6 + UPPS.11 + UPPS.16 + UPPS.21 + UPPS.28 + UPPS.33 + UPPS.38 + UPPS.43 + UPPS.48 + UPPS.55)
upps$LackPers <- with(upps, UPPS.4 + UPPS.9 + UPPS.14 + UPPS.19 + UPPS.24 + UPPS.27 + UPPS.32 + UPPS.37 + UPPS.42 + UPPS.47)
upps$SenSeek <- with(upps, UPPS.3 + UPPS.8 + UPPS.13 + UPPS.18 + UPPS.23 + UPPS.26 + UPPS.31 + UPPS.36 + UPPS.41 + UPPS.46 + UPPS.51 + UPPS.56)
upps$PosUrg <- with(upps, UPPS.5 + UPPS.10 + UPPS.15 + UPPS.20 + UPPS.25 + UPPS.30 + UPPS.35 + UPPS.40 + UPPS.45 + UPPS.49 + UPPS.52 + UPPS.54 + UPPS.57 + UPPS.59)
upps$upps_total <- with(upps, NegUrg + LackPrem + LackPers, SenSeek, PosUrg)

upps.complete <- complete.cases(dplyr::select(upps,-SPECC_ID, -LUNA_ID, -BPD, -Age, -ACID, -CompletionDate))
upps.complete <- upps[upps.complete,]
upps.complete <- dplyr::select(upps.complete,-SPECC_ID, -LUNA_ID, -BPD, -Age, -ACID, -CompletionDate)




# DERS --------------------------------------------------------------------

ders <- read.csv("DERS.csv", header=TRUE)


ders$SPECC_ID <- sub("_", "", as.character(ders$SPECC_ID)) 
ders$LUNA_ID <- as.character(ders$LUNA_ID)

ders <- subj_info %>% left_join(ders, by = c( "SPECC_ID" = "SPECC_ID", "LUNA_ID" = "LUNA_ID"))
ders <- ders[-22,]

dersItems <- paste("DERS", 1:36, sep=".")

#check for missing items
for (i in 1:nrow(ders)) {
  miss <- which(is.na(ders[i, dersItems]))
  if (length(miss) > 0) {
    cat("For subject id: ", paste(sapply(ders[i, c("LUNA_ID", "SPECC_ID")], as.character), collapse="/"), ", these items are missing: ", paste(miss, collapse=", "), "\n", sep="")
  }
}
#reverse score items
reverseItems <-paste("DERS", c(1, 2, 6, 7, 8, 10, 17, 20, 22, 24, 34), sep=".")
ders[,reverseItems] <- lapply(ders[,reverseItems], function(x) { 6 - x })

#SPECC ID 074DS: missing item 22
ders <- meanreplace(ders, "074DS",varprefix = "DERS.", subscaleitems = c(16,15,31,35,28,22,36,30))

#SPECC ID 089VF: missing item 36
ders <- meanreplace(ders, "089VF", varprefix = "DERS.", subscaleitems = c(16,15,31,35,28,22,36,30))

#SPECC ID 140JP: missing item 26
ders <- meanreplace(ders, "140JP", varprefix = "DERS.", subscaleitems = c(26,18,13,33,20))



#define variables and score items
ders$NonAccept <- with(ders, DERS.25 + DERS.21 + DERS.12 + DERS.11 + DERS.29 + DERS.23)
ders$Goals <- with(ders, DERS.26 + DERS.18 + DERS.13 + DERS.33 + DERS.20)
ders$Impulse <- with(ders, DERS.32 + DERS.27 + DERS.14 + DERS.19 + DERS.3 + DERS.24)
ders$Aware <- with(ders, DERS.6 + DERS.2 + DERS.10 + DERS.17 + DERS.8 + DERS.34)
ders$Strategy <- with(ders, DERS.16 + DERS.15 + DERS.31 + DERS.35 + DERS.28 + DERS.22 + DERS.36 + DERS.30)
ders$Clarity <- with(ders, DERS.5 + DERS.4 + DERS.9 + DERS.7 + DERS.1)
ders$ders_total <- with(ders, NonAccept + Goals, Impulse, Aware, Strategy, Clarity)


ders.complete <- complete.cases(dplyr::select(ders,-SPECC_ID, -LUNA_ID, -BPD, -Age, -ACID, -CompletionDate))
ders.complete <- ders[ders.complete,]
ders.complete <- dplyr::select(ders.complete,-SPECC_ID, -LUNA_ID, -BPD, -Age, -ACID, -CompletionDate)


# BPQ ---------------------------------------------------------------------

BPQ <- read.csv("BPQScored.csv", header = TRUE)


BPQ$SPECC_ID <- sub("_", "", as.character(BPQ$SPECC_ID)) 
BPQ$LUNA_ID <- as.character(BPQ$LUNA_ID)

BPQ <- subj_info %>% left_join(BPQ, by = c( "SPECC_ID" = "SPECC_ID", "LUNA_ID" = "LUNA_ID"))
bpq.complete <- complete.cases(dplyr::select(BPQ,-SPECC_ID, -LUNA_ID, -BPD, -Age, -ACID, -CompletionDate))
bpq.complete <- BPQ[bpq.complete,]
bpq.complete <- dplyr::select(bpq.complete,-SPECC_ID, -LUNA_ID, -BPD, -Age, -ACID, -CompletionDate)



# PANAS -------------------------------------------------------------------

panas <- read.csv("PANAS.csv", header=TRUE)


panas$SPECC_ID <- sub("_", "", as.character(panas$SPECC_ID)) 
panas$LUNA_ID <- as.character(panas$LUNA_ID)

panas <- subj_info %>% left_join(panas, by = c( "SPECC_ID" = "SPECC_ID", "LUNA_ID" = "LUNA_ID"))


panasItems <- paste("PANAS", 1:20, sep=".")

#check for missing items
for (i in 1:nrow(panas)) {
  miss <- which(is.na(panas[i, panasItems]))
  if (length(miss) > 0) {
    cat("For subject id: ", paste(sapply(panas[i, c("LUNA_ID", "SPECC_ID")], as.character), collapse="/"), ", these items are missing: ", paste(miss, collapse=", "), "\n", sep="")
  }
}


#SPECC ID 051KM: missing item 7
panas <- meanreplace(panas, "051KM",varprefix = "PANAS.", subscaleitems = c(2,4,6,7,8,11,13,15,18,20))

#define variables and score items. doesn't get used bc constructing single factor with no subscales
panas$negaff <- with(panas, PANAS.2+ PANAS.4 + PANAS.6 + PANAS.7 +PANAS.8 + PANAS.11 + PANAS.13 + PANAS.15 + PANAS.18 + PANAS.20) 


panas.negaff.complete <- complete.cases(dplyr::select(panas,-SPECC_ID, -LUNA_ID, -BPD, -Age, -ACID, -CompletionDate))
panas.negaff.complete <- panas[panas.negaff.complete,]
panas.negaff.complete <- dplyr::select(panas.negaff.complete,-SPECC_ID, -LUNA_ID, -BPD, -Age, -ACID, -CompletionDate)

# IIP ---------------------------------------------------------------------


IIP <- read.csv("IIP.csv", header = TRUE)

IIP$SPECC_ID <- sub("_", "", as.character(IIP$SPECC_ID)) 
IIP$LUNA_ID <- as.character(IIP$LUNA_ID)

IIP <- subj_info %>% left_join(IIP, by = c( "SPECC_ID" = "SPECC_ID", "LUNA_ID" = "LUNA_ID"))


IIPItems <- c(paste("IIP", 1:10, sep="."), paste("IIP", 12:90, sep = "."))

#check for missing items
for (i in 1:nrow(IIP)) {
  miss <- which(is.na(IIP[i, IIPItems]))
  if (length(miss) > 0) {
    cat("For subject id: ", paste(sapply(IIP[i, c("LUNA_ID", "SPECC_ID")], as.character), collapse="/"), ", these items are missing: ", paste(miss, collapse=", "), "\n", sep="")
  }
}

IIP <- meanreplace(IIP, "015CW",varprefix = "IIP.", subscaleitems = c(seq(1,10,1), seq(12,90,1)))
IIP <- meanreplace(IIP, "066DW",varprefix = "IIP.", subscaleitems = c(seq(1,10,1), seq(12,90,1)))
IIP <- meanreplace(IIP, "074DS",varprefix = "IIP.", subscaleitems = c(seq(1,10,1), seq(12,90,1)))
IIP <- meanreplace(IIP, "118ZB",varprefix = "IIP.", subscaleitems = c(seq(1,10,1), seq(12,90,1)))
IIP <- meanreplace(IIP, "138AH",varprefix = "IIP.", subscaleitems = c(seq(1,10,1), seq(12,90,1)))


IIP$IIP_PD1=with(IIP, (IIP.1 +  IIP.35 + IIP.36 + IIP.42 + IIP.51 + IIP.55 +IIP.60 +IIP.78 +IIP.79 +IIP.81 + IIP.86)/11)
IIP$IIP_PD2=with(IIP, (IIP.13 + IIP.14 +IIP.26 + IIP.28 + IIP.32 +IIP.34 +IIP.38 +IIP.40 +IIP.41+IIP.84)/10)
IIP$IIP_PD3=with(IIP, (IIP.50 + IIP.53 + IIP.58 +IIP.63 +IIP.77 +IIP.80 +IIP.88)/7)
IIP$IIP_bpd=with(IIP, (IIP.51 +IIP.53 +IIP.55 +IIP.66 +IIP.77 +IIP.80+IIP.89+IIP.90)/8)

IIP$PA <- with(IIP, (IIP.21 + IIP.40 + IIP.57 + IIP.58 + IIP.65 + IIP.68 + IIP.76 + IIP.80))
IIP$BC=with(IIP, (IIP.1 + IIP.26 + IIP.28 + IIP.38 + IIP.41 + IIP.50 + IIP.73 + IIP.88))
#81341 has missing value at 11 (DE) and 74 which does not need to be imputed (not needed for octant scales)
IIP$DE=with(IIP, (IIP.18 +IIP.20 + IIP.24 + IIP.27 + IIP.31 + IIP.46 + IIP.82))  #+IIP.11  
IIP$FG=with(IIP, (IIP.3 + IIP.7 + IIP.17 + IIP.22 + IIP.43 + IIP.45 + IIP.71 + IIP.85))
#81341 has missing value at 5 which needs to be imputed
IIP$HI=with(IIP, (IIP.5 + IIP.6 + IIP.8 + IIP.9 + IIP.12 + IIP.15 + IIP.23 + IIP.49))
IIP$JK=with(IIP, (IIP.2 +IIP.10 + IIP.29 + IIP.44 + IIP.48 + IIP.54 +IIP.69 + IIP.83)) 
#80430, impute q 72, q64
IIP$LM=with(IIP, (IIP.25 +  IIP.37+ IIP.47 +  IIP.59 + IIP.64 +  IIP.67 + IIP.70 + IIP.87))
IIP$NO=with(IIP, (IIP.4 + IIP.30 + IIP.39 + IIP.52 + IIP.56 + IIP.61 + IIP.62 + IIP.78))
#81011, impute 42 
IIP$IIP_PD1=with(IIP, (IIP.1 +  IIP.35 + IIP.36 + IIP.42 + IIP.51 + IIP.55 +IIP.60 +IIP.78 +IIP.79 +IIP.81 + IIP.86)/11)
IIP$IIP_PD2=with(IIP, (IIP.13 + IIP.14 +IIP.26 + IIP.28 + IIP.32 +IIP.34 +IIP.38 +IIP.40 +IIP.41+IIP.84)/10)
IIP$IIP_PD3=with(IIP, (IIP.50 + IIP.53 + IIP.58 +IIP.63 +IIP.77 +IIP.80 +IIP.88)/7)
IIP$IIP_bpd=with(IIP, (IIP.51 +IIP.53 +IIP.55 +IIP.66 +IIP.77 +IIP.80+IIP.89+IIP.90)/8)
#8090 and 8134 need to impute
IIP$IIP_agency <- with(IIP, .25*(PA - HI + .707*(BC + NO - FG - JK)))
#8043 and 8143 need to impute
IIP$IIP_communion <- with(IIP, .25*(LM - DE + .707*(NO + JK - BC - FG)))
#8043, 8090, and 8143 need to impute
IIP$IIP_elevation <- with(IIP, (PA + BC + DE + FG + HI + JK + LM + NO)/8)



IIP.complete <- complete.cases(dplyr::select(IIP,-SPECC_ID, -LUNA_ID, -BPD, -Age, -ACID, -CompletionDate))
IIP.complete <- IIP[IIP.complete,]
IIP.complete <- dplyr::select(IIP.complete,-SPECC_ID, -LUNA_ID, -BPD, -Age, -ACID, -CompletionDate)


# compile -----------------------------------------------------------------

selfreports <- BPQ %>% left_join(ders, by = c( "SPECC_ID" = "SPECC_ID", "LUNA_ID" = "LUNA_ID", "BPD" = "BPD")) %>% 
  left_join(panas, by = c( "SPECC_ID" = "SPECC_ID", "LUNA_ID" = "LUNA_ID", "BPD" = "BPD")) %>% 
  left_join(upps, by = c( "SPECC_ID" = "SPECC_ID", "LUNA_ID" = "LUNA_ID", "BPD" = "BPD")) %>%
  left_join(IIP, by = c( "SPECC_ID" = "SPECC_ID", "LUNA_ID" = "LUNA_ID", "BPD" = "BPD"))

selfreports <- selfreports %>% select(SPECC_ID, LUNA_ID, BPD, Age, impulsivity, instability, abandonment, relationships, selfimage, suicide, emptiness, 
                                      anger, psychotic, totalscore, NonAccept, Goals, Impulse, Aware, Strategy, Clarity,ders_total,
                                      NegUrg, LackPrem, LackPers, SenSeek, PosUrg, upps_total, negaff, IIP_PD1, IIP_PD2, IIP_PD3, IIP_bpd, PA, BC, 
                                      DE, FG, HI, JK, LM, NO, IIP_agency, IIP_communion, IIP_elevation)
setwd(basedir)
