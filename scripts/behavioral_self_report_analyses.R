
# behavioral and self-report regressions -------------------------------------
library(foreign)
suppressMessages(require(ggplot2)); suppressMessages(require(multcomp)); suppressMessages(require(wesanderson)); suppressMessages(require(ggpubr))
setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
#initialize the graph analysis pipeline (includes sourcing pipeline inputs, creating graphs, thresholding, binarization, community assignment, and calculation and reduction of global and nodal graph metrics)

source("scripts/setup_globals.R")
#RIDGE
inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, preproc_pipeline = "nosmooth_aroma_bp_nonaggr") #leave blank for defaults, consult with function for deviations from default
#PEARSON
#inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, conn_method = "pearson", rs_desired_log = logspace(log10(.2), log10(.4), 20)) #leave blank for defaults, consult with function for deviations from default
for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)

###load brain data
toanalyze_thresh <- get(load("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/cache/toanalyze.fa.thresh_schaefer422_nosmooth_aroma_bp_nonaggr_ridge.net_partial_fc_binary_all.RData"))
for(i in 1:length(toanalyze_thresh)) assign(names(toanalyze_thresh)[i], toanalyze_thresh[[i]])
path <- paste0(basedir, "/cache/membership_df.RData")
membership_df <- get(load("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/cache/membership_df.RData"))

membership_df <- membership_df %>%
  mutate(membership=factor(membership, levels=1:7, labels=c("VIS", "SOMMOT", "DORSATTN", "SALVENTATTN", "LIMBIC", "FPN", "DMN")))

toanalyze <- left_join(toanalyze, membership_df, by="node") %>%
  mutate(BPD=factor(BPD, levels=0:1, labels=c("Control", "BPD")))
str(toanalyze)

fa.metrics <- colnames(dplyr::select(toanalyze, -id, -node, -BPD, -Age, -membership))
# Crossed-effects MLM allowing for Empirical Bayes Estimates to be estimated for each subject at each network for each connectivity metric--------

network_sums <- data.frame(id = subj_info$SPECC_ID)

for(i in fa.metrics){
  f <- as.formula(paste0(i, " ~ 1 + membership + (1 + membership|id) + (1|node)"))
  metric <- lmer(f, toanalyze)
  emp_bayes <- ranef(metric)$id
  
  metric_scores <- emp_bayes
  for(j in 2:7) {
    metric_scores[,j] <- metric_scores[,1] + metric_scores[,j] #get rid of offsets
  }
  colnames(metric_scores) <- paste0(i, "_", seq(1,7,1)); metric_scores <- add_rownames(metric_scores, "id")
  
  
  network_sums <- left_join(network_sums, metric_scores, by = "id")
  
}

####if column names bomb out:
# central_scores <- paste0("central_", seq(1,7,1))
# integration_scores <- paste0("integration_", seq(1,7,1))
# within.mod_scores <- paste0("within.mod_", seq(1,7,1))
# 
# allnames <- c("id", central_scores, integration_scores, within.mod_scores)
# colnames(network_sums) <- allnames
# head(network_sums)

# self-report data --------------------------------------------------------
selfreports <- read.csv("data/allselfreports.csv")

#dmn_toanalyze <- toanalyze %>% dplyr::filter(membership == 7) %>% left_join(selfreports, by = c("id", "BPD", "Age"))

#combine self-report and network sums
network_sums <- network_sums %>% left_join(selfreports, by = "id")

# load task data ----------------------------------------------------------
emo_gender_data <- read.csv("~/Box Sync/DEPENd/Projects/Emotion and Gender Conflict BPD/readytofit_noinacc_no_outliers_no_na_centered_interaction_invage.csv")
bandito <- read.spss("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/data/hallquist_bandit_data_n95.sav", to.data.frame = TRUE)



# test mediation models DMN and IIP---------------------------------------------------
library(lavaan)

dmn_iip <- "
IIP_elevation ~ b*central_7 + c*BPD
central_7 ~ a*BPD
ab:=a*b
"
med_dmn_iip <- sem(model = dmn_iip, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_iip)

dmn_integration_iip <- "
IIP_elevation ~ b*integration_7 + c*BPD
integration_7 ~ a*BPD
ab:=a*b
"
med_dmn_integration_iip <- sem(model = dmn_integration_iip, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_integration_iip)

dmn_central_iip_bpd <- "
IIP_bpd ~ b*central_7 + c*BPD
central_7 ~ a*BPD
ab:=a*b
"
med_dmn_central_iip_bpd <- sem(model = dmn_central_iip_bpd, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_central_iip_bpd)

dmn_within_iip_bpd <- "
IIP_bpd ~ b*within.mod_7 + c*BPD
within.mod_7 ~ a*BPD
ab:=a*b
"
med_dmn_within_iip_bpd <- sem(model = dmn_within_iip_bpd, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_within_iip_bpd)

dmn_integration_iip_bpd <- "
IIP_bpd ~ b*integration_7 + c*BPD
integration_7 ~ a*BPD
ab:=a*b
"
med_dmn_integration_iip_bpd <- sem(model = dmn_integration_iip_bpd, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_integration_iip_bpd)

dmn_integration_iip_bpd <- "
IIP_PD1 ~ b*within.mod_7 + c*BPD
within.mod_7 ~ a*BPD
ab:=a*b
"
med_dmn_integration_iip_bpd <- sem(model = dmn_integration_iip_bpd, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_integration_iip_bpd)


# test mediation models DMN and BPQ---------------------------------------------------
dmn_integration_bpq_total <- "
totalscore ~ b*integration_7 + c*BPD
integration_7 ~ a*BPD
ab:=a*b
"
med_dmn_integration_bpq_total <- sem(model = dmn_integration_bpq_total, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_integration_bpq_total)

dmn_central_bpq_total <- "
totalscore ~ b*central_7 + c*BPD
central_7 ~ a*BPD
ab:=a*b
"
med_dmn_central_bpq_total <- sem(model = dmn_central_bpq_total, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_central_bpq_total)

dmn_within.mod_bpq_total <- "
totalscore ~ b*within.mod_7 + c*BPD
within.mod_7 ~ a*BPD
ab:=a*b
"
med_dmn_within.mod_bpq_total <- sem(model = dmn_within.mod_bpq_total, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_within.mod_bpq_total)

dmn_central_bpq_selfimage <- "
selfimage ~ b*central_7 + c*BPD
central_7 ~ a*BPD
ab:=a*b
"
med_dmn_central_bpq_selfimage <- sem(model = dmn_central_bpq_selfimage, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_central_bpq_selfimage)

dmn_integration_bpq_selfimage <- "
selfimage ~ b*integration_7 + c*BPD
integration_7 ~ a*BPD
ab:=a*b
"
med_dmn_integration_bpq_selfimage <- sem(model = dmn_integration_bpq_selfimage, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_integration_bpq_selfimage)

dmn_within.mod_bpq_selfimage <- "
selfimage ~ b*within.mod_7 + c*BPD
within.mod_7 ~ a*BPD
ab:=a*b
"
med_dmn_within.mod_bpq_selfimage <- sem(model = dmn_within.mod_bpq_selfimage, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_within.mod_bpq_selfimage)

dmn_central_bpq_emptiness <- "
emptiness ~ b*central_7 + c*BPD
central_7 ~ a*BPD
ab:=a*b
"
med_dmn_central_bpq_emptiness <- sem(model = dmn_central_bpq_emptiness, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_central_bpq_emptiness)

dmn_within.mod_bpq_emptiness <- "
emptiness ~ b*central_7 + c*BPD
central_7 ~ a*BPD
ab:=a*b
"
med_dmn_central_bpq_emptiness <- sem(model = dmn_central_bpq_emptiness, data = network_sums, mimic= "Mplus", se = "bootstrap", bootstrap = 5000, estimator = "ML", missing = "ML", meanstructure = TRUE)
summary(med_dmn_central_bpq_emptiness)


colnames(network_sums)

# simple linear associations for testing ----------------------------------------------
summary(lm(IIP_elevation ~ within.mod_7, data = network_sums))
summary(lm(IIP_bpd ~ within.mod_7, data = network_sums))
summary(lm(selfimage ~ within.mod_7, data = network_sums))
summary(lm(IIP_agency ~ within.mod_7, data = network_sums))
summary(lm(IIP_communion ~ within.mod_7, data = network_sums))
summary(lm(negaff ~ within.mod_7, data = network_sums))
summary(lm(Clarity ~ within.mod_7, data = network_sums))
summary(lm(instability ~ within.mod_7, data = network_sums))
summary(lm(suicide ~ within.mod_7, data = network_sums))
summary(lm(anger ~ within.mod_7, data = network_sums))
summary(lm(psychotic ~ within.mod_7, data = network_sums))
summary(lm(ders_total ~ within.mod_7, data = network_sums))
summary(lm(Goals ~ within.mod_7, data = network_sums))
summary(lm(Impulse ~ within.mod_7, data = network_sums))
summary(lm(Aware ~ within.mod_7, data = network_sums))
summary(lm(Strategy ~ within.mod_7, data = network_sums))
summary(lm(IIP_PD1 ~ within.mod_7, data = network_sums))
summary(lm(IIP_PD2 ~ within.mod_7, data = network_sums))
summary(lm(IIP_PD3 ~ within.mod_7, data = network_sums))

# summary(lm(central_5 ~ 1 +BPD*Age, data = network_sums))
# summary(lm(IIP_bpd ~ BPD, data = network_sums))
# 
# summary()
