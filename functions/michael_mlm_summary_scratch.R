results_for_comparison <- mlm_res

rslope <- results_for_comparison[[1]]
cmod <- rslope[[1]]
summary(cmod) #centrality
car::Anova(cmod)
imod <- rslope[[2]]
summary(imod) #integration
car::Anova(imod)
wmod <- rslope[[3]]
summary(wmod) #withinmod
car::Anova(wmod)
