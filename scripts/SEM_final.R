# read in graph theory data and clean up for entry into single FAs --------
library("lavaan"); library("semPlot")
setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()

# Read in self report and combine with graph data -----------------

setwd(basedir)
source("scripts/gather_self_reports_sem.R")
setwd(basedir)
source("scripts/compile_network_factor_scores.R")

#combine selfreports with connectivity estimates
selfreports$DMN_conn <- DMN_scores
selfreports$FPN_conn <- FPN_scores
selfreports$LIMBIC_conn <- LIMBIC_scores
selfreports$DORSAL_conn <- DORSAL_scores

#drop 001_RA who is missing all self report data
selfreports <- selfreports[-1,]


toanalyze_sem <- selfreports %>% dplyr::select(SPECC_ID, 
                                               impulsivity, instability, abandonment, relationships, selfimage, suicide,emptiness,anger, #BPQ
                                               NegUrg, LackPrem, LackPers, SenSeek, #UPPS
                                               NonAccept, Goals, Impulse, Aware, Strategy, Clarity, #DERS
                                               PANAS.2, PANAS.4, PANAS.6, PANAS.7, PANAS.8, PANAS.11, PANAS.13, PANAS.15, PANAS.18, PANAS.20, #PANAS
                                               DMN_conn, FPN_conn, LIMBIC_conn, DORSAL_conn #graph connectivity
                                               )


# Regress network factor scores on self reports of interest ---------------
colnames(selfreports)
dmn_bpq <- lm(DMN_conn ~ totalscore, data = selfreports)
dmn_na <- lm(DMN_conn ~ negaff, data = selfreports)
summary(dmn_bpq)



# examine single factor models --------------------------------------------
#BPQ
BPQ_mod <- '
BPD =~  abandonment + relationships + selfimage + emptiness + suicide 
selfimage ~~ emptiness
'
bpq.data <- toanalyze_sem %>% dplyr::select(abandonment, relationships, selfimage, suicide, emptiness)


bpq.sem <- sem(BPQ_mod, data = bpq.data, estimator = "MLR")#,  missing = "FIML") 
semPaths(bpq.sem)
summary(bpq.sem, fit.measures = TRUE)
# modificationindices(bpq.sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)
rowMeans(cor(bpq.data, use = "complete.obs"))

##DERS
DERS_mod <- '
DERS =~ NonAccept + Goals  + Aware + Strategy 
'

# ders.data <- toanalyze_sem %>% dplyr::select(NonAccept, Goals, Aware, Strategy, Clarity)


ders.sem <- sem(DERS_mod, data = toanalyze_sem, estimator = "MLR")#,  missing = "FIML") 
semPaths(ders.sem)
summary(ders.sem, fit.measures = TRUE)


##UPPS
upps_mod <- '
upps =~ NegUrg + LackPrem + LackPers 
'

upps.sem <- sem(upps_mod, data = toanalyze_sem, estimator = "MLR")#,  missing = "FIML") 
semPaths(upps.sem)
summary(upps.sem, fit.measures = TRUE)
# modificationindices(bpq.sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)
#rowMeans(cor(upps.data, use = "complete.obs"))
#resid(upps.sem, type = "cor")


##PANAS_NA
panas_mod <- '
PANAS_NA =~ PANAS.2 + PANAS.4  + PANAS.8 + PANAS.15 + PANAS.18 + PANAS.20
'
panas.sem <- sem(panas_mod, data = toanalyze_sem, estimator = "MLR")#,  missing = "FIML") 
semPaths(panas.sem)
summary(panas.sem, fit.measures = TRUE)

#diagnostics
# modificationindices(panas.sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)
# resid(panas.sem, type = "cor")

##IIP
IIP_bpd_mod <- '
IIP_bpd =~ IIP.51  +IIP.55 +IIP.66 +IIP.80
'
IIP_bpd.sem <- sem(IIP_bpd_mod, data = selfreports, estimator = "MLR")#,  missing = "FIML") 
semPaths(IIP_bpd.sem)
summary(IIP_bpd.sem, fit.measures = TRUE)

#diagnostics
modificationindices(IIP_bpd.sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)
resid(IIP_bpd.sem, type = "cor")

##ends up in final model
IIP_bpd_mod <- '
IIP_bpd =~ IIP.51  +IIP.55  +IIP.80
'
IIP_bpd.sem <- sem(IIP_bpd_mod, data = selfreports, estimator = "MLR")#,  missing = "FIML") 
semPaths(IIP_bpd.sem)
summary(IIP_bpd.sem, fit.measures = TRUE)

##IIP_PD1 (interpersonal hypersensitivity)
IIP_pd1_mod <- '
IIP_pd1 =~ IIP.1  + IIP.51 + IIP.55 +IIP.60 +IIP.78 +IIP.79  
IIP.51 ~~IIP.55
'
IIP_pd1.sem <- sem(IIP_pd1_mod, data = selfreports, estimator = "MLR")#,  missing = "FIML") 
# semPaths(IIP_bpd.sem)
summary(IIP_pd1.sem, fit.measures = TRUE)

#diagnostics
modificationindices(IIP_pd1.sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)
resid(IIP_pd1.sem, type = "cor")


# test full model ---------------------------------------------------------

model1_full <- '
BPD =~  abandonment + relationships + selfimage + emptiness + suicide 
selfimage ~~ emptiness
DERS =~ NonAccept + Goals  + Aware + Strategy 
UPPS =~ NegUrg + LackPrem + LackPers 
PANAS_NA =~ PANAS.2 + PANAS.4 + PANAS.6 + PANAS.8 + PANAS.15 + PANAS.18 + PANAS.20
UPPS ~ FPN_conn
DERS ~ FPN_conn
PANAS_NA ~ DMN_conn
PANAS_NA ~ LIMBIC_conn
DERS ~ LIMBIC_conn
DERS ~ PANAS_NA
BPD ~ DERS
BPD ~ UPPS
UPPS ~ DERS
LackPrem ~~ LackPers
PANAS.15 ~~ PANAS.20
Goals ~~ LackPers

'

mod1.sem <- sem(model1_full, data = toanalyze_sem, estimator = "MLR",  missing = "FIML") 
semPaths(mod1.sem)
summary(mod1.sem, fit.measures = TRUE)
modificationindices(mod1.sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)


# add IIP ---------------------------------------------------------

model2_full <- '
BPD =~  abandonment + relationships + selfimage + emptiness + suicide 
selfimage ~~ emptiness
DERS =~ NonAccept + Goals + Strategy 
UPPS =~ NegUrg + LackPrem + LackPers 
PANAS_NA =~ PANAS.2 + PANAS.4  + PANAS.8 + PANAS.15 + PANAS.18 + PANAS.20
IIP_bpd =~ IIP.51  +IIP.55  +IIP.80 
UPPS ~ FPN_conn
DERS ~ FPN_conn
PANAS_NA ~ DMN_conn
PANAS_NA ~ LIMBIC_conn
DERS ~ LIMBIC_conn
IIP_bpd ~ DMN_conn
IIP_bpd ~ FPN_conn
IIP_bpd ~ PANAS_NA
IIP_bpd ~ UPPS
BPD ~ IIP_bpd
DERS ~ PANAS_NA
BPD ~ DERS
BPD ~ UPPS
UPPS ~ DERS
LackPrem ~~ LackPers
PANAS.15 ~~ PANAS.20
Goals ~~ LackPers
'

mod2.sem <- sem(model2_full, data = selfreports, estimator = "MLR",  missing = "FIML") 
#semPaths(mod2.sem)
summary(mod2.sem, fit.measures = TRUE)
1247.99modificationindices(mod2.sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)


# change to IIP interpersonal sensitivity ---------------------------------

model3_full <- '
BPD =~  abandonment + relationships + selfimage + emptiness + suicide 
selfimage ~~ emptiness
DERS =~ NonAccept + Goals + Strategy 
UPPS =~ NegUrg + LackPrem + LackPers 
PANAS_NA =~ PANAS.2 + PANAS.4  + PANAS.8 + PANAS.15 + PANAS.18 + PANAS.20
IIP_pd1 =~ IIP.1  + IIP.55 +IIP.60 +IIP.78 +IIP.79  
UPPS ~ FPN_conn
DERS ~ FPN_conn
PANAS_NA ~ DMN_conn
PANAS_NA ~ LIMBIC_conn
DERS ~ LIMBIC_conn
IIP_pd1 ~ DMN_conn
IIP_pd1 ~ FPN_conn
IIP_pd1 ~ PANAS_NA
BPD ~ IIP_pd1
DERS ~ IIP_pd1
DERS ~ PANAS_NA
BPD ~ DERS
BPD ~ UPPS
UPPS ~ DERS
LackPrem ~~ LackPers
PANAS.15 ~~ PANAS.20
Goals ~~ LackPers
'

mod3.sem <- sem(model3_full, data = selfreports, estimator = "MLR",  missing = "FIML") 
#semPaths(mod3.sem)
summary(mod3.sem, fit.measures = TRUE)
modificationindices(mod3.sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)

anova(mod1.sem,mod2.sem,mod3.sem)
