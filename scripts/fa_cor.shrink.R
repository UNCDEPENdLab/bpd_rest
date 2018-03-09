
# factor analysis from transformed metrics cor.shrink ---------------------

nodalmetrics_df_test <- get(load("/Users/nth7/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/cache/threshnodalmetrics_transformed_schaefer422_nosmooth_aroma_bp_nonaggr_cor.shrink_fc_binary_all.RData"))

fa.CFI<-function(x){
  nombre<-paste(x,"CFI",sep = ".")
  nombre<-
    ((x$null.chisq-x$null.dof)-(x$STATISTIC-x$dof))/(x$null.chisq-x$null.dof)
  return(nombre)
}

anyNA(nodalmetrics_df)

# head(dplyr::filter(nodalmetrics_df, density == 0.75))
# anyNA(nodalmetrics_df$part.coeff)

nodalmetrics_df[is.na(nodalmetrics_df)] <- 0

metrics.raw <- nodalmetrics_df %>% dplyr::select(id, node, density, one_of(reducemetrics)) %>%
  filter(density > den) %>% gather(key="variable", value="value", -id, -node, -density)
#recast for FA such that column names represent a metrics at a given density (e.g., 0.05_between.module.deg.zscore)
metrics.raw_fa <- dcast(metrics.raw, id + node ~ density + variable, value.var = "value")
toanalyze <- dplyr::select(metrics.raw_fa, id, node)


# individual centrality scores --------------------------------------------


# betweenness -------------------------------------------------------------
# btw.scores <- metrics.raw[which(metrics.raw$variable == "betweenness.node"),]  
# btw.fa.obj <- dcast(btw.scores, id+node~density+variable, value.var = "value" )
# faout <- fa(dplyr::select(btw.fa.obj, -id, -node), nfactors = 3, fm = "ml", rotate = "promax")
# print(faout$loadings, cutoff = 0.3)
# fasolution.btw <- data.frame(faout$scores)
# names(fasolution.btw) <- "betweenness"
# fasolution.btw.scores <- cbind(dplyr::select(btw.fa.obj, id, node), fasolution.btw)
# 
# toanalyze <- left_join(toanalyze, fasolution.btw.scores, by = c("id", "node"))

# degree ------------------------------------------------------------------
# ##betweenness centrality
# degree.scores <- metrics.raw[which(metrics.raw$variable == "degree"),]  
# degree.fa.obj <- dcast(degree.scores, id+node~density+variable, value.var = "value" )
# faout <- fa(dplyr::select(degree.fa.obj, -id, -node), nfactors = 1, fm = "ml", rotate = "promax")
# print(faout$loadings, cutoff = 0.3)
# fasolution.btw <- data.frame(faout$scores)
# names(fasolution.btw) <- "betweenness"
# fasolution.btw.scores <- cbind(dplyr::select(btw.fa.obj, id, node), fasolution.btw)
# 
# toanalyze <- left_join(toanalyze, fasolution.btw.scores, by = c("id", "node"))



merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)

toanalyze <- dplyr::left_join(toanalyze, merge.bpdage, by = "id")


# multiple graph metrics --------------------------------------------------


#####retain for analysis, 4 factors shows imporvements in model fit with low factor loadings on factor 4
#pcaout <- pca(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 3,  missing = TRUE, rotate = "promax")#, scores = "Bartlett")
faout <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 3,  missing = TRUE, rotate = "promax")#, scores = "Bartlett")
summary <- summary(faout)
CFI <- fa.CFI(faout)
fasolution <- data.frame(faout$scores)
names(fasolution) <- c("central","integration", "within.mod")
toanalyze_comb <- cbind(toanalyze, fasolution)

head(toanalyze_comb)


vss.out <- vss(dplyr::select(metrics.raw_fa, -id, -node), n = 3)

CFI <- rbind(CFI, SRMR = vss.out[["vss.stats"]]$SRMR[3])



if(!exists("subj_info")){subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)}
merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)

# browser()

toanalyze <- dplyr::left_join(toanalyze_comb, merge.bpdage, by = "id")


# testing -----------------------------------------------------------------


# faout_2 <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 2,  missing = TRUE, rotate = "varimax")
# print(faout_2$loadings, cutoff = .3)
# summary(faout_2)
# fa.CFI(faout_2)
# 
# faout_3 <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 3,  missing = TRUE, rotate = "promax")
# print(faout_3$loadings, cutoff = .3)
# summary(faout_3)
# fa.CFI(faout_3)
# 
# faout_4 <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 4,  missing = TRUE, rotate = "promax")
# print(faout_4$loadings, cutoff = .3)
# summary(faout_4)
# fa.CFI(faout_4)
# 
# 
# faout_5 <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 5,  missing = TRUE, rotate = "promax")
# print(faout_4$loadings, cutoff = .3)
# summary(faout_4)
# fa.CFI(faout_4)
# 
