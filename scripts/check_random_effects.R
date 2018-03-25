# Initialize --------------------------------------------------------------

setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
#initialize the graph analysis pipeline (includes sourcing pipeline inputs, creating graphs, thresholding, binarization, community assignment, and calculation and reduction of global and nodal graph metrics)

source("scripts/setup_globals.R")

#specify pipeline inputs.leave blank for defaults, consult with function for deviations from default#RIDGE
inputs <- specify_inputs(thresh_weighted = "binary", 
                        conn_method = "cor.shrink",                         
                         fc_out_rm = FALSE, 
                         preproc_pipeline = "nosmooth_aroma_bp_nonaggr",
                         reducemetrics =  c( "degree", "page.rank", "part.coeff", "eigen.cent", "gateway.coeff.btw", "gateway.coeff.degree", "within.module.deg")) 

for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

# Prep toanalyze ------------------------------------
# toanalyze_thresh <- get(load("/Users/nth7/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/cache/toanalyze.fa.thresh_schaefer422_nosmooth_aroma_bp_nonaggr_ridge.net_partial_fc_binary_all.RData"))
toanalyze_thresh <- get(load("toanalyze_cor.shrink_final_questionmark.RData"))
for(i in 1:length(toanalyze_thresh)) assign(names(toanalyze_thresh)[i], toanalyze_thresh[[i]])

toanalyze$Age_c <- toanalyze$Age - mean(toanalyze$Age)


library(MASS)
library(car)
library(lme4)
library(multcomp)
library(lattice)

membership_df <- data.frame(node = rownames(do.call(cbind,yeo7)), membership = as.numeric(do.call(cbind,yeo7)[,1]))
membership_df <- membership_df %>%
  mutate(membership=factor(membership, levels=1:7, labels=c("VIS", "SOMMOT", "DORSATTN", "SALVENTATTN", "LIMBIC", "FPN", "DMN")))

toanalyze <- left_join(toanalyze, membership_df, by="node") %>%
  mutate(BPD=factor(BPD, levels=0:1, labels=c("Control", "BPD")))


# int.t<- Winsorize(toanalyze$integration, minval = quantile(toanalyze$integration, .01), maxval = quantile(toanalyze$integration, .9999))
# int.t <- int.t +4
# toanalyze$integration <- BoxCox_extract(int.t, seq(-3,25,.1))

# check random effects structure ------------------------------------------


attach(toanalyze)
m1 <- lmer(within.mod ~1 + (1|id))
m2 <- lmer(within.mod ~1 +(1|node) + (1|id))
#m3 <- lmer(central ~1 +(1|node) + (1|id) + (1|membership))
m4 <- lmer(within.mod ~1 +(1 | node) + (1 + membership|id))
m5 <- lmer(within.mod ~1 + (1 + membership|id))
# m5 <- lmer(within.mod ~1 +(1|node) + (1 + node|id))
# m6 <- lmer(within.mod ~1 +(1|node) + (1 + membership*node|id))
#m5 <- lmer(central ~ 1 + (1|node) + (1+ membership*BPD|id))
anova(m2, m4)

ranef(m4)$id


# calculate ICC -----------------------------------------------------------

summary(m2)

ICC <- .5484/(.5484+.1168+.3525)
ICC

library(sjPlot); library(sjmisc)

# plot random slopes ------------------------------------------------------
sjp.lmer(m4, type = "rs.ri")#, sample.n = 15)

library(lme4)
m.int <- lmer(integration ~1 +(1|node) + (1 + membership|id), toanalyze)
sjp.lmer(m.int, type = "rs.ri")

m.within <- lmer(within.mod ~1 +(1|node) + (1 + membership|id), toanalyze)
sjp.lmer(m.within, type = "rs.ri")
