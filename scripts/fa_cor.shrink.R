
# factor analysis from transformed metrics cor.shrink ---------------------

nodalmetrics_df <- get(load("/Users/natehall/Box Sync/bpd_rest/cache/threshnodalmetrics_transformed_MARCH_light_trans_schaefer422_nosmooth_aroma_bp_nonaggr_cor.shrink_fc_binary_all_update2018.RData"))

#untransformed
#nodalmetrics_df <- allmetrics.nodal.df



nodalmetrics_df[is.na(nodalmetrics_df)] <- 0
metrics.raw <- nodalmetrics_df %>% dplyr::select(id, node, wthresh, one_of(reducemetrics)) %>%
   gather(key="variable", value="value", -id, -node, -wthresh)

#recast for FA such that column names represent a metrics at a given density (e.g., 0.05_between.module.deg.zscore)
metrics.raw_fa <- dcast(metrics.raw, id + node ~ wthresh + variable, value.var = "value")
#convert non-finite values to 0
is.na(metrics.raw_fa)<-sapply(metrics.raw_fa, is.infinite)
metrics.raw_fa[is.na(metrics.raw_fa)]<-0

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
faout <- fa(dplyr::select(metrics.raw_fa, -id, -node), nfactors = 3, rotate = "oblimin", fm = "pa")#, scores = "Bartlett")
# summary <- summary(faout)
# CFI <- fa.CFI(faout)
fasolution <- data.frame(faout$scores)
names(fasolution) <- c("integration", "within.mod", "betweenness")
toanalyze_comb <- cbind(toanalyze, fasolution)

# head(toanalyze_comb)


# vss.out <- vss(dplyr::select(metrics.raw_fa, -id, -node), n = 3)

# CFI <- rbind(CFI, SRMR = vss.out[["vss.stats"]]$SRMR[3])



if(!exists("subj_info")){subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)}
merge.bpdage <- subj_info %>% dplyr::select(SPECC_ID, BPD, AgeAtScan) %>% dplyr::rename(id=SPECC_ID, Age=AgeAtScan)

# browser()

toanalyze <- dplyr::left_join(toanalyze_comb, merge.bpdage, by = "id")

toanalyze_melt <- melt(toanalyze, id.vars = c("id", "node", "BPD", "Age"))
head(toanalyze_melt)
ggplot(toanalyze_melt, aes(x=value)) + geom_histogram() + facet_wrap(~variable, scales = "free") 
