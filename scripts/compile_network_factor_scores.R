setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()
#initialize the graph analysis pipeline (includes sourcing pipeline inputs, creating graphs, thresholding, binarization, community assignment, and calculation and reduction of global and nodal graph metrics)
if(!exists("toanalyze_binary")){source("scripts/rs_initialize_pipeline.R")}



# DMN ---------------------------------------------------------------------
dmn_reduced.thresh <- reduce_networks_forsem(toanalyze_binary, 7, .15)

mod_network <- '
DMN =~ V183 + V164 + V184 + V389 +V385 + V165 + V386 + V178 + V197
V385 ~~ V178
V183 ~~ V389
'

dmn_sem <- sem(mod_network, dmn_reduced.thresh, estimator = "MLR")
semPaths(dmn_sem)
summary(dmn_sem, fit.measures = TRUE)
resid(dmn_sem, type = "cor")
modificationindices(dmn_sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)
inspect(dmn_sem)

DMN_scores <- predict(dmn_sem)


# FPN ---------------------------------------------------------------------
fpn_reduced.thresh <- reduce_networks_forsem(toanalyze_binary, 6, .15)

mod_network <- '
FPN =~ V349 + V336 + V345 + V333 + V137 + V337 + V130
'

fpn_sem <- sem(mod_network, fpn_reduced.thresh, estimator = "MLR")
semPaths(fpn_sem)
summary(fpn_sem, fit.measures = TRUE)
# resid(dmn_sem, type = "cor")
# modificationindices(dmn_sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)
# inspect(dmn_sem)

FPN_scores <- predict(fpn_sem)


# LIMBIC ------------------------------------------------------------------
limbic_reduced.thresh <- reduce_networks_forsem(toanalyze_binary, 5, prop.thresh = .2)

mod_network <- '
LIMBIC =~ V326 + V119 + V404 + V322 + V403 + V327
V404 ~~ V403
'

limbic_sem <- sem(mod_network, limbic_reduced.thresh, estimator = "MLR")
summary(limbic_sem, fit.measures = TRUE)
# #diagnostics
# modificationindices(limbic_sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)
# resid(limbic_sem, type = "cor")

LIMBIC_scores <- predict(limbic_sem)


# Dorsal Attention --------------------------------------------------------

dorsal_reduced.thresh <- reduce_networks_forsem(toanalyze_binary, 3, prop.thresh = .15)

mod_network <- '
dorsal =~ V291 + V289 + V75 + V73  + V286
'

dorsal_sem <- sem(mod_network, dorsal_reduced.thresh, estimator = "MLR")
summary(dorsal_sem, fit.measures = TRUE)
# #diagnostics
# modificationindices(dorsal_sem) %>% dplyr::filter(mi > 2) %>% arrange(desc(mi)) %>% head(20)
# resid(dorsal_sem, type = "cor")

DORSAL_scores <- predict(limbic_sem)



