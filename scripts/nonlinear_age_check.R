
# Model comparisons for nonlinear age - related changes in network --------

if (!exists("basedir")){setwd("~/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/"); basedir <- getwd()}
####read in package dependencies and custom functions
#source("scripts/rs_initialize_pipeline.R") #setup globals, get adjmats, setup graphs, label communities, reduce across thresholds, calculate graph metrics


source("scripts/setup_globals.R")
#RIDGE
inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, preproc_pipeline = "nosmooth_aroma_bp_nonaggr") #leave blank for defaults, consult with function for deviations from default
#PEARSON
#inputs <- specify_inputs(thresh_weighted = "binary", fc_out_rm = FALSE, conn_method = "pearson", rs_desired_log = logspace(log10(.2), log10(.4), 20)) #leave blank for defaults, consult with function for deviations from default
for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements

subj_info <- get_subj_info(adjmats_base, parcellation, conn_method, preproc_pipeline, file_extension=".txt.gz", fd.scrub = TRUE, allowCache = TRUE)


require(rlang)
require(lme4)
# 
# toanalyze_c <- toanalyze_thresh %>% mutate(membership=factor(membership), BPD = factor(BPD))#, levels=1:7, labels=c("VIsual")))
# metrics <- colnames(toanalyze_c)
# metrics <- toanalyze_c %>% dplyr::select(-id, -node, -BPD, -Age, -membership)
# toanalyze_c$membership <- mapvalues(toanalyze_c$membership, from = c("1","2","3","4","5","6","7"), to = c("VIS", "SOMMOT", "DORSATTN", "SALVENTATTN", "LIMBIC", "FPN", "DMN"))

###load brain data
toanalyze_thresh <- get(load("/Users/nth7/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/cache/toanalyze.fa.thresh_schaefer422_nosmooth_aroma_bp_nonaggr_ridge.net_partial_fc_binary_all.RData"))
for(i in 1:length(toanalyze_thresh)) assign(names(toanalyze_thresh)[i], toanalyze_thresh[[i]])
path <- paste0(basedir, "/cache/membership_df.RData")
membership_df <- get(load("/Users/nth7/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/cache/membership_df.RData"))
toanalyze <- left_join(toanalyze, membership_df, by = "node")
toanalyze$id <- as.character(toanalyze$id)

metrics <- toanalyze %>% dplyr::select(-id,-node,-BPD,-Age,-membership)

metrics_results <- list()
anovas <- list()

pdf(file = "age_nonlinear_checkbynetwork_feb2018.pdf", width = 10, height = 7)
for(m in colnames(metrics)){
  for (net in unique(toanalyze$membership)){
    
    age_id <- dplyr::select(subj_info, c("SPECC_ID", "AgeAtScan")) 
    names(age_id) <- c("id", "Age")
    
    toanalyze_this.net <- toanalyze[which(toanalyze$membership == net),] %>% dplyr::select(id, node, Age, m)
    
    
    f_lin <- as.formula(paste(m, " ~ Age + (1 + membership|id) + (1|node)"))
    linear_lmer <- lmer(f_lin, data = toanalyze_this.net)
    
    f_asymp <- as.formula(paste(m, " ~ I(1/Age) + (1 + membership|id) + (1|node)"))
    asymptotic_lmer <- lmer(f_asymp, data = toanalyze_this.net)
    
    # f_quad <- as.formula(paste(m, " ~ poly(Age, 2) + (1|id) + (1|node)"))
    # quad_lmer <- lmer(f_quad, data = toanalyze_this.net)
    # 
    res <- anova(linear_lmer, asymptotic_lmer)
    
    toplot_this.net <- toanalyze[which(toanalyze$membership == net),] %>% group_by(membership, id) %>% 
      summarise(mean = mean(!!sym(m)), sd = sd(!!sym(m))) %>% left_join(age_id, by = "id")
    
    age_plot <- ggplot() + geom_point(data = toplot_this.net, aes(x = Age, y = mean)) + labs(title = net, y = m) +
      geom_smooth(data = toplot_this.net, aes(x=Age, y=mean), method="lm", formula = y ~ x, n= 40, se=TRUE, color="blue", size=1.5) +
      geom_smooth(data = toplot_this.net, aes(x=Age, y=mean), method="lm", formula = y ~ I(1/x), n= 40, se=TRUE, color="green", size=1.5) +
      geom_smooth(data = toplot_this.net, aes(x=Age, y=mean), method="lm", formula = y ~ poly(x,2), n= 40, se=TRUE, color="red", size=1.5)
    
    plot(age_plot)
    
    anovas[[net]] <- res
  }
  metrics_results[[m]] <- anovas 
}
dev.off()


a <- do.call(rbind, metrics_results)
b <- do.call(rbind, a)
b
# unlist(a)
# 
# data.frame(b)
# summary(anovas)
# str(anovas)
# str(anovas)
# 
# summary(res)
# res
# 
# stats::anova.lm(linear_lm, asymptotic_lm, quad_lm, test = "F")
# 
# require(nonnest2)
# 
# vuongtest(linear_lm, asymptotic_lm)


