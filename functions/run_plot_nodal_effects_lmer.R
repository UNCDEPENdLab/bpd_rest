run_plot_nodal_effects_lmer <- function (toanalyze, allowCache = TRUE, browse = FALSE, plot = TRUE){
  require(multcomp)
  stopifnot(file.exists(file.path(basedir, "results"))) 
  expectFile <- file.path(basedir, "results", paste0(parcellation, "_", conn_method),   paste0("all.sigttest.nodal.binary_lmer", data.reduce, ".RData"))
  if (file.exists(expectFile.ttest) && file.exists(expectFile.lm) && file.exists(expectFile.wle) && allowCache==TRUE) {
    message("Loading binary nodal metric comparisons from: ", expectFile)
    results_all <- load(expectFile)
  } else {
    metrics.toanalyze <- names(dplyr::select(toanalyze, -id,-node,-BPD,-Age))
    for(m in metrics.toanalyze){
      f <- as.formula(paste(m, "~ 1 + node*BPD*Age + (1|id)"))
      m1 <- lmer(f, toanalyze)
      cmat <- pairs(lsmeans(~BPD | node, m1))
      forglht <- cmat@linfct
      
      summary(glht(m1, linfct=forglht), test=adjusted("single-step")) 
    }
    
    
  }
}

toanalyze$BPD <- as.factor(toanalyze$BPD) 
str(toanalyze)
m_test <- lmer(central ~ 1 + node*BPD + (1 | id), toanalyze)
# cmat_test <- pairs(lsmeans(~BPD | node, m_test))
cmat_test <- pairs(lsmeans(m_test, ~BPD | node))
forglht <- cmat_test@linfct

test_glht <- summary(glht(m_test, linfct=forglht), test=adjusted("single-step")) 
beepr::beep()

mlm_res_nodal <- list()
for (i in fa.metrics){
  f <- as.formula(paste(i, "~ 1 + node*BPD*Age + (1|id)"))
  model <- lmer(f, toanalyze)
  mlm_res_nodal[[i]] <- model 
}

save(file = paste0(basedir, "/cache/mlm_res_nodal_02_21_18.RData"), mlm_res_nodal)
beepr::beep()

table_glht <- function(x) {
  #pq <- summary(x)$test
  pq <- x
  mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
  error <- attr(pq$pvalues, "error")
  pname <- switch(x$alternativ, less = paste("Pr(<", ifelse(x$df ==0, "z", "t"), ")", sep = ""), 
                  greater = paste("Pr(>", ifelse(x$df == 0, "z", "t"), ")", sep = ""), two.sided = paste("Pr(>|",ifelse(x$df == 0, "z", "t"), "|)", sep = ""))
  colnames(mtests) <- c("Estimate", "Std. Error", ifelse(x$df ==0, "z value", "t value"), pname)
  return(mtests)
  
}

test_table <- table_glht(test_glht)

test_glht$test

cmat_test

m <- metrics.toanalyze[1]
