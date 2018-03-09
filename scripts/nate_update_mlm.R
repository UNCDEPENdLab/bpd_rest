# Initialize --------------------------------------------------------------

setwd("~/Box Sync/bpd_rest/"); basedir <- getwd()
#initialize the graph analysis pipeline (includes sourcing pipeline inputs, creating graphs, thresholding, binarization, community assignment, and calculation and reduction of global and nodal graph metrics)

source("scripts/setup_globals.R")
library(MASS)
library(car)
library(lme4)
library(multcomp)
library(lattice)

#specify pipeline inputs.leave blank for defaults, consult with function for deviations from default#RIDGE
inputs <- specify_inputs(thresh_weighted = "binary",
                         conn_method = "cor.shrink",
                         fc_out_rm = FALSE, 
                         preproc_pipeline = "nosmooth_aroma_bp_nonaggr",
                         reducemetrics =  c( "degree", "page.rank", "part.coeff", "eigen.cent", "gateway.coeff.btw", "gateway.coeff.degree", "within.module.deg")) 

for(i in 1:length(inputs)) assign(names(inputs)[i], inputs[[i]]) #assign objects to names of input list elements


# Prep toanalyze ------------------------------------
# toanalyze_thresh <- get(load("/Users/nth7/Box Sync/DEPENd/Projects/RS_BPD_graph/bpd_rest/cache/toanalyze.fa.thresh_schaefer422_nosmooth_aroma_bp_nonaggr_ridge.net_partial_fc_binary_all.RData"))
# for(i in 1:length(toanalyze_thresh)) assign(names(toanalyze_thresh)[i], toanalyze_thresh[[i]])

load("toanalyze_cor.shrink_final_questionmark.RData")


toanalyze$Age_c <- toanalyze$Age - mean(toanalyze$Age)

membership_df <- get(load("~/Box Sync/bpd_rest/cache_protected/membership_df.RData"))
membership_df <- membership_df %>%
  mutate(membership=factor(membership, levels=1:7, labels=c("VIS", "SOMMOT", "DORSATTN", "SALVENTATTN", "LIMBIC", "FPN", "DMN")))

toanalyze <- left_join(toanalyze, membership_df, by="node") %>%
  mutate(BPD=factor(BPD, levels=0:1, labels=c("Control", "BPD")))



# qplot(log(toanalyze$within.mod + 2))
# 
# toanalyze$within.mod <- log(toanalyze$within.mod + 2)
# int.t<- Winsorize(toanalyze$integration, minval = quantile(toanalyze$integration, .01), maxval = quantile(toanalyze$integration, .9999))
# int.t <- int.t +4
# toanalyze$integration <- BoxCox_extract(int.t, seq(-3,25,.1))


# check normality ---------------------------------------------------------

# 
# pdf("normality_check_hists.pdf", height = 8, width =11)
# histogram( ~betweenness | node, toanalyze)
# histogram( ~integration | node, toanalyze)
# histogram( ~within.mod | node, toanalyze)
# dev.off()

# str(toanalyze)

# analyze separate networks (from MH) -------------------------------------

sink("sepnetworks_NH.txt")
pdf("sepnetworks_NH.pdf", width=12, height=6)
networks <- levels(toanalyze$membership)
for (n in networks) {
  cat("\n========\nNETWORK", n, "\n")
  df <- droplevels(filter(toanalyze, membership==n))
  df$Age_c <- df$Age - mean(df$Age, na.rm=TRUE)
  df$Agecat <- factor(df$Age < 18)
  df$Age_inv <- -1000*(1/df$Age - mean(1/df$Age, na.rm=TRUE)) #inverse age
  df$Age_quad <- I(df$Age^2)
  
  dfagg <- df %>% group_by(id) %>% dplyr::summarize(
    betweenness=median(betweenness), 
    integration=median(integration),
    within.mod=median(within.mod), BPD=head(BPD, n=1), Age=head(Age, n=1), Age_inv = head(Age_inv,1), Age_quad = head(Age_quad,1))
  
  dfcomb <- dfagg %>% gather(key="metric", value="value", betweenness, integration, within.mod)
  g <- ggplot(dfcomb, aes(x=Age, y=value, color=BPD)) + geom_point() + stat_smooth(method="lm", aes(fill = BPD), alpha=.2) + facet_wrap(~metric, scales="free_y", nrow=1) +
    ggtitle(paste("Network: ", n)) + theme_bw() + scale_color_brewer(palette = "Set1");plot(g)
  
  for (m in c("betweenness", "integration", "within.mod")) {
    cat("\n -- METRIC: ", m, "\n\n")

    f <- as.formula(paste(m, "~ 1  + BPD*Age_c + (1|node) + (1 |id)"))# + BPD*I(Age^2)"))
    #f <- as.formula(paste(m, "~ 1 + BPD*Agecat + (1|node) + (1|id)")) # + BPD*I(Age^2)
    mod <- lmer(f, df)

    cmat <- list(age=list(), group=list(), agexgroup=list())
    #membershipLevels <- levels(df$membership)
    coefs.rs <- fixef(mod)
    ageEff <- "Age_c"
    bpdEff <- "BPDBPD"
    bpd_age <- paste(bpdEff, ageEff, sep = ":")

      #setup empty contrast
      cont.rs <- rep(0, length(coefs.rs))
      names(cont.rs) <- names(coefs.rs)
      cont.rs[1] <- 1 #intercept

      ###age contrast
      age.cont <- cont.rs
      age.cont[ageEff] <- 1 #always include main effect of age for age tests

      conName <- paste(m, "Age", sep = "x")
      cmat$age[[conName]] <- age.cont

      ###bpd contrast
      bpd.cont <- cont.rs
      bpd.cont[bpdEff] <- 1

      conName <- paste(m, "BPD", sep = "x")
      cmat$group[[conName]] <- bpd.cont

      #####
      #bpd x age contrast
      agexbpd.cont <- cont.rs
      agexbpd.cont[bpd_age] <- 1 #BPD x age interaction, reference

      conName <- paste(m, "Age", "BPD", sep=" x ")
      cmat$agexgroup[[conName]] <- agexbpd.cont

    #print(summary(m1))
    cat("\nOmnibus\n")
    print(car::Anova(mod))

    cat("\nContrasts\n")
    cmat.mat <- do.call(rbind, lapply(cmat, function(el) do.call(rbind, el)))

    ##separate into list of 3
    cmat.matSeparate <- lapply(cmat, function(el) do.call(rbind, el))

    ###simultaneous inference
    results <- lapply(cmat.matSeparate, function(effect) {
      summary(glht(mod, linfct=effect), test=adjusted("single-step"))
    })

    print(results)
  }
  cat("\n===========\n")
}
sink()
dev.off()



# run for loop over metrics with and without random slope estimated-----------------------------------------------


sink("MLMresults.txt")
#pdf("MLMplots.pdf", width=11, height=8)
for (m in c("betweenness", "integration", "within.mod")) {
  df <- droplevels(toanalyze)
  cat("\n===========\n")
  cat("\n -- METRIC: ", m, "\n\n")
  
  f <- as.formula(paste(m, "~ 1  + BPD*Age_c*membership + (1|node) + (1 + membership|id)"))# + BPD*I(Age^2)"))
  #f <- as.formula(paste(m, "~ 1 + BPD*Agecat + (1|node) + (1|id)")) # + BPD*I(Age^2)
  mod <- lmer(f, df)
  
  cat("\nOmnibus\n")
  print(car::Anova(mod))
  
  cmat <- list(age=list(), group=list(), agexgroup=list())
  membershipLevels <- levels(df$membership)
  coefs.rs <- fixef(mod)
  ageEff <- "Age_c"
  bpdEff <- "BPDBPD"
  
  ##setup contrasts
  for(n in membershipLevels){
    #setup empty contrast
    cont.rs <- rep(0, length(coefs.rs))
    names(cont.rs) <- names(coefs.rs)
    
    thisnet <- paste0("membership", n)
    ##two-way interactions
    netbpd <- paste(bpdEff, thisnet, sep = ":")
    netage <- paste(ageEff, thisnet, sep = ":")
    bpd_age <- paste(bpdEff, ageEff, sep = ":")
    
    #three-way interaction:
    netbpdage <- paste(bpdEff, ageEff, thisnet, sep = ":")
    
    ###age contrast
    age.cont <- cont.rs
    age.cont[ageEff] <- 1 #always include main effect of age for age tests
    
    #age x network two-way interaction
    if (netage %in% names(coefs.rs)) age.cont[netage] <- 1
    
    conName <- paste(n, "Age", sep = "x")
    cmat$age[[conName]] <- age.cont
    
    ###bpd contrast
    bpd.cont <- cont.rs
    bpd.cont[bpdEff] <- 1
    
    #bpd x network two-way interaction
    if (netbpd %in% names(coefs.rs)) bpd.cont[netbpd] <- 1
    
    conName <- paste(n, "BPD", sep = "x")
    cmat$group[[conName]] <- bpd.cont
    
    #####
    #bpd x age contrast
    agexbpd.cont <- cont.rs
    agexbpd.cont[bpd_age] <- 1 #BPD x age interaction, reference
    
    #three-way interactions that include age x comt
    if (netbpdage %in% names(coefs.rs)) agexbpd.cont[netbpdage] <- 1
   
    
    conName <- paste(n, "Age", "BPD", sep=" x ")
    cmat$agexgroup[[conName]] <- agexbpd.cont
  }
 
  ##combine all contrasts 
  cmat.mat <- do.call(rbind, lapply(cmat, function(el) do.call(rbind, el)))
  
  ##separate into list of 3
  cmat.matSeparate <- lapply(cmat, function(el) do.call(rbind, el))
  
  ###simultaneous inference
  results <- lapply(cmat.matSeparate, function(effect) {
    summary(glht(mod, linfct=effect), test=adjusted("single-step"))
  })
  
  # results_uncorrected <- lapply(cmat.matSeparate, function(effect) {
  #   summary(glht(mod, linfct=effect), test=adjusted("none"))
  # })
  
  
  cat("\n===========\nResults when including random slopes for membership and using simultaneous inference procedure from Hothorn 2008\n")
  print(results)
  cat("\n===========\n")
  # cat("\n===========\nUncorrected\n")
  # print(results_uncorrected)
  # cat("\n===========\n")
  # 

# remove random slope ----------------------------------------------------

  
  f2 <- as.formula(paste(m, "~ 1  + BPD*Age_c*membership + (1|node) + (1|id)"))# + BPD*I(Age^2)"))
  #f <- as.formula(paste(m, "~ 1 + BPD*Agecat + (1|node) + (1|id)")) # + BPD*I(Age^2)
  mod2 <- lmer(f2, df)
  
  cat("\nOmnibus\n")
  print(car::Anova(mod))
  results <- lapply(cmat.matSeparate, function(effect) {
    summary(glht(mod2, linfct=effect), test=adjusted("single-step"))
  })
  
  # results_uncorrected <- lapply(cmat.matSeparate, function(effect) {
  #   summary(glht(mod2, linfct=effect), test=adjusted("single-step"))
  # })
  
  cat("\n===========\nResults when removing random slopes for membership and using simultaneous inference procedure from Hothorn 2008\n")
  print(results)
  cat("\n===========\n")
  # cat("\n===========\nUncorrected\n")
  # print(results_uncorrected)
  # cat("\n===========\n")
  
  
}
sink()


# nodal effects -----------------------------------------------------------
# m1 <- lmer(betweenness ~ 1 + node*Age_c*BPD  + (1|id), df)
# cmat <- pairs(lsmeans(m1,~BPD | node))
# forglht <- cmat@linfct
# summary(glht(m1, linfct=forglht), test=adjusted("single-step"))

sink("MLMresults_nodal_single_step.txt")
#pdf("MLMplots.pdf", width=11, height=8)
for (m in c("betweenness", "integration", "within.mod")) {
  df <- droplevels(toanalyze)
  
  cat("\n===========\n")
  cat("\n -- METRIC: ", m, "\n\n")
  
  f <- as.formula(paste(m, "~ 1 + node*Age_c*BPD +(1|id)"))# + BPD*I(Age^2)"))
  #f <- as.formula(paste(m, "~ 1 + BPD*Agecat + (1|node) + (1|id)")) # + BPD*I(Age^2)
  mod <- lmer(f, df)
  
  cat("\nOmnibus\n")
  print(car::Anova(mod))
  
  cmat <- list(age=list(), group=list(), agexgroup=list())
  nodeLevels <- levels(df$node)
  coefs.rs <- fixef(mod)
  ageEff <- "Age_c"
  bpdEff <- "BPDBPD"
  
  ##setup contrasts
  for(n in nodeLevels){
    #setup empty contrast
    cont.rs <- rep(0, length(coefs.rs))
    names(cont.rs) <- names(coefs.rs)
    
    thisnode <- paste0("node", n)
    ##two-way interactions
    nodebpd <- paste(thisnode,bpdEff, sep = ":")
    nodeage <- paste(thisnode, ageEff, sep = ":")
    bpd_age <- paste(ageEff, bpdEff, sep = ":")
    
    #three-way interaction:
    nodebpdage <- paste(thisnode, ageEff, bpdEff, sep = ":")
    
    ###age contrast
    age.cont <- cont.rs
    age.cont[ageEff] <- 1 #always include main effect of age for age tests
    
    #age x node two-way interaction
    if (nodeage %in% names(coefs.rs)) age.cont[nodeage] <- 1
    
    conName <- paste(n, "Age", sep = "x")
    cmat$age[[conName]] <- age.cont
    
    ###bpd contrast
    bpd.cont <- cont.rs
    bpd.cont[bpdEff] <- 1
    
    #bpd x node two-way interaction
    if (nodebpd %in% names(coefs.rs)) bpd.cont[nodebpd] <- 1
    
    conName <- paste(n, "BPD", sep = "x")
    cmat$group[[conName]] <- bpd.cont
    
    #####
    #bpd x age contrast
    agexbpd.cont <- cont.rs
    agexbpd.cont[bpd_age] <- 1 #BPD x age interaction, reference
    
    #three-way interactions 
    if (nodebpdage %in% names(coefs.rs)) agexbpd.cont[nodebpdage] <- 1
    
    
    conName <- paste(n, "Age", "BPD", sep=" x ")
    cmat$agexgroup[[conName]] <- agexbpd.cont
  }
  
  ##combine all contrasts 
  cmat.mat <- do.call(rbind, lapply(cmat, function(el) do.call(rbind, el)))
  
  ##separate into list of 3
  cmat.matSeparate <- lapply(cmat, function(el) do.call(rbind, el))
  
  ###simultaneous inference
  results <- lapply(cmat.matSeparate, function(effect) {
    summary(glht(mod, linfct=effect), test=adjusted("single-step"))
  })
  
  results_fdr <- lapply(cmat.matSeparate, function(effect) {
    summary(glht(mod, linfct=effect), test=adjusted("fdr"))
  })
  
  
  cat("\n===========\nResults for nodal analyses using simultaneous inference procedure from Hothorn 2008\n")
  print(results)
  cat("\n===========\n")
  # cat("\n===========\nFDR corrected:\n")
  # 
  # 
  # print(results_fdr)
  # 
  # cat("\n===========\n")

  
}
beepr::beep("coin")
sink()

# 
# toanalyze_nest <- toanalyze %>% gather(key="metric", value="value", betweenness, integration, within.mod) %>%
#   group_by(node, metric) %>% tidyr::nest()
# 
# fit_node <- function(data) {
#   data$Age <- data$Age - mean(data$Age, na.rm=TRUE)
#   mod <- lm(value ~ 1 + BPD*Age, data = data)
#   return(mod)
# }
# 
# results <- toanalyze_nest %>% dplyr::mutate(mod = purrr::map(data, fit_node))
# results <- results %>% dplyr::mutate(tidy = purrr::map(mod, broom::tidy))
# 
# effects <- results %>%
#   tidyr::unnest(tidy) %>%
#   dplyr::select(-std.error) %>%
#   filter(term != "(Intercept)")
# #gather(key="modelpart", value="value", estimate, statistic, p.value) %>%
# #spread(term, value) %>%
# #tidyr::spread(term, estimate)
# 
# ##global info about sig p values
# psig <- effects %>% filter(p.value < .05) %>% arrange(term, metric, node, p.value) %>%
#   inner_join(membership_df)
# 
# table(psig$membership, psig$term)
# 
# psigfdr <- effects %>% group_by(term) %>% mutate(padj = p.adjust(p.value, method="fdr")) %>% ungroup() %>%
#   filter(padj < .05) %>% inner_join(membership_df)
# 
# table(psigfdr$membership, psigfdr$term)
