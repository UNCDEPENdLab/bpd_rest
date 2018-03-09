#setwd(file.path(getMainDir(), "Miscellaneous", "Aarthi_ROIMultComp"))
setwd(paste0(basedir, "/multcomp_from_MH"))

library(lme4)
library(multcomp)
library(lattice)
rois <- read.csv("MATRIX_ALL.csv")

#looks normal enough to me!
histogram(~X | PFC + Subcortical, rois)

#store subject id as categorical
rois$Subject.ID <- factor(rois$Subject.ID)

#mixed correlation model
corrMixed <- lmer(X ~ 1 + PFC * Subcortical * COMT.Genotype * AgeInverseCentered + (1 | Subject.ID), rois, REML=TRUE)

fixef(corrMixed)



#compute cmat for each PFC * Subcort combination
cmat <- list(age=list(), genotype=list(), agexgeno=list())
subcortLevels <- levels(rois$Subcortical)
pfcLevels <- levels(rois$PFC)
coefs <- fixef(corrMixed)
ageEff <- "AgeInverseCentered"
comtEff <- "COMT.Genotype"

for (l in 1:length(subcortLevels)) {
  for (m in 1:length(pfcLevels)) {
    
    #setup empty contrast
    cont <- rep(0, length(coefs))
    names(cont) <- names(coefs)
    
    #setup possible effects of interest
    meSubcort <- paste("Subcortical", subcortLevels[l], sep="")
    mePFC <- paste("PFC", pfcLevels[m], sep="")
    
    #interactions must follow order of model specification above
    #two-way interactions (6 of them)
    pfc_subcort <- paste(mePFC, meSubcort, sep=":")
    pfc_comt <- paste(mePFC, comtEff, sep=":")
    pfc_age <- paste(mePFC, ageEff, sep=":")
    subcort_comt <- paste(meSubcort, comtEff, sep=":")
    subcort_age <- paste(meSubcort, ageEff, sep=":") 
    comt_age <- paste(comtEff, ageEff, sep=":")
    
    #three-way interactions (4 of them)
    pfc_subcort_comt <- paste(mePFC, meSubcort, comtEff, sep=":")
    pfc_subcort_age <- paste(mePFC, meSubcort, ageEff, sep=":")
    pfc_comt_age <- paste(mePFC, comtEff, ageEff, sep=":")
    subcort_comt_age <- paste(meSubcort, comtEff, ageEff, sep=":")
    
    #four-way interaction
    fourway <- paste(mePFC, meSubcort, comtEff, ageEff, sep=":")
    
    #####
    #age contrast
    #Generally: age + subcort:age + pfc:age + subcort:pfc:age
    age.cont <- cont
    age.cont[ageEff] <- 1 #always include main effect of age for age tests
    
    #two-way age x brain region interactions
    if (pfc_age %in% names(coefs)) age.cont[pfc_age] <- 1
    if (subcort_age %in% names(coefs)) age.cont[subcort_age] <- 1
    
    #three-way age x subcort x pfc
    if (pfc_subcort_age %in% names(coefs)) age.cont[pfc_subcort_age] <- 1
    
    conName <- paste(subcortLevels[l], pfcLevels[m], "Age", sep=" X ")
    cmat$age[[conName]] <- age.cont
    
    #####
    #genotype contrast
    #Generally: comt + subcort:comt + pfc:comt + subcort:pfc:comt
    genotype.cont <- cont
    genotype.cont[comtEff] <- 1 #always include main effect of genotype for genotype tests
    
    #two-way genotype x brain region interactions
    if (pfc_comt %in% names(coefs)) genotype.cont[pfc_comt] <- 1
    if (subcort_comt %in% names(coefs)) genotype.cont[subcort_comt] <- 1
    
    #three-way genotype x subcort x pfc
    if (pfc_subcort_comt %in% names(coefs)) genotype.cont[pfc_subcort_comt] <- 1
    
    conName <- paste(subcortLevels[l], pfcLevels[m], "COMT", sep=" X ")
    cmat$genotype[[conName]] <- genotype.cont
    
    #####
    #age x genotype
    #Contrasts test how the age x comt interaction deviates for this PFC x Subcort pair
    #Generally: age:comt + subcort:age:comt + pfc:age:comt + subcort:pfc:age:comt
    agexgeno.cont <- cont
    agexgeno.cont[comt_age] <- 1 #always include genotype x age interaction, which is the reference
    
    #three-way interactions that include age x comt
    if (pfc_comt_age %in% names(coefs)) agexgeno.cont[pfc_comt_age] <- 1
    if (subcort_comt_age %in% names(coefs)) agexgeno.cont[subcort_comt_age] <- 1
    
    #four-way interaction
    if (fourway %in% names(coefs)) agexgeno.cont[fourway] <- 1
    
    conName <- paste(subcortLevels[l], pfcLevels[m], "Age", "COMT", sep=" x ")
    cmat$agexgeno[[conName]] <- agexgeno.cont
    
    
  }
}

#check genotype contrast problem:


#print out contrast estimates for a couple of models.
for (i in 1:2) {
  cat("Age cont: ", names(cmat$age)[i], "\n")
  print(sum(coefs[which(cmat$age[[i]] > 0)]))
  
  cat("Comt cont: ", names(cmat$comt)[i], "\n")
  print(sum(coefs[which(cmat$genotype[[i]] > 0)]))
  
  cat("Age x geno cont: ", names(cmat$agexgeno)[i], "\n")
  print(sum(coefs[which(cmat$agexgeno[[i]] > 0)]))
}

#spot check three models 
laccum_acc <- subset(rois, Subcortical=="L_Accum" & PFC=="ACC")
summary(m1 <- lm(X ~ AgeInverseCentered * COMT.Genotype, data=laccum_acc))
coefs[which(cmat$agexgeno[["L_Accum x ACC x Age x COMT"]] > 0)]

lm(X ~ AgeInverseCentered, data=laccum_acc)
lm(X ~ COMT.Genotype, data=laccum_acc)

m1c <- coef(m1)

pred.low <- m1c["(Intercept)"] + m1c["COMT.Genotype"]

laccum_dacc <- subset(rois, Subcortical=="L_Accum" & PFC=="DACC")
summary(lm(X ~ AgeInverseCentered * COMT.Genotype, data=laccum_dacc))
coefs[which(cmat$agexgeno[["L_Accum x DACC x Age x COMT"]] > 0)]

lcaudate_rsfg <- subset(rois, Subcortical=="L_Caudate" & PFC=="RightSFG")
summary(lm(X ~ AgeInverseCentered * COMT.Genotype, data=lcaudate_rsfg))
coefs[which(cmat$agexgeno[["L_Caudate x RightSFG x Age x COMT"]] > 0)]
sum(coefs[which(cmat$agexgeno[["L_Caudate x RightSFG x Age x COMT"]] > 0)])

#get model-estimated correlations for each combination (useful for graphs)
cm <- lmerCellMeans(corrMixed, divide=c("AgeInverseCentered", "COMT.Genotype"))

####SIMULTANEOUS TESTS
#combine into a matrix
#this combines into all contrasts
cmat.mat <- do.call(rbind, lapply(cmat, function(el) do.call(rbind, el)))

#separated by age, genotype, and interaction
cmat.matSeparate <- lapply(cmat, function(el) do.call(rbind, el))


#no correction (per effect). Adapt below accordingly
lapply(cmat.matSeparate, function(effect) {
      print(summary(glht(corrMixed, linfct=effect), test=adjusted("none")))
    })

###
(nocorrect <- summary(glht(corrMixed, linfct=cmat.mat), test=adjusted("none")))
#print only significant effects
nocorrect$test$coefficients[which(summary(glht(corrMixed, linfct=cmat.mat), test=adjusted("none"))$test$pvalues < .05)]

#if you did the crazy bonferroni
summary(glht(corrMixed, linfct=cmat.mat), test=adjusted("bonferroni"))

#holm is no better
summary(glht(corrMixed, linfct=cmat.mat), test=adjusted("holm"))

#liberal FDR correction
summary(glht(corrMixed, linfct=cmat.mat), test=adjusted("fdr"))

#single-step correction (the Hothorn paper, and multcomp's default
summary(glht(corrMixed, linfct=cmat.mat), test=adjusted("single-step"))

#free combinations of parameters (Westfall 1999)
summary(glht(corrMixed, linfct=cmat.mat), test=adjusted("free"))

#logical constraint adjustments
summary(glht(corrMixed, linfct=cmat.mat), test=adjusted("Westfall"))

#logical constraint adjustments
summary(glht(corrMixed, linfct=cmat.mat), test=adjusted("Shaffer"))



###


cm <- lmerCellMeans(corrMixed, n.cont=c(AgeInverseCentered=30), fixat0="COMT.Genotype")


ggplot(cm, aes(x=AgeInverseCentered, y=X, color=factor(COMT.Genotype))) +
    facet_grid(PFC ~ Subcortical) + geom_line()

#helper function to get predictions of model-estimated correlations

lmerCellMeans <- function(lmerObj, divide=NULL, n.divide=3, divide.prefix=TRUE, n.cont=30, fixat0=NULL) { 
  #print cell means for lmer by expanding level combinations and multiplying against fixed effects  
  predNames <- attr(terms(lmerObj), "term.labels")
  whichME <- attr(terms(lmerObj), "order")
  predNames <- predNames[whichME==1]
  
  predData <- list()
  #divide into categorical and continuous predictors, determine which continuous predictors to discretize
  for (f in predNames) {
    if (f %in% fixat0) { predData[[f]] <- 0 #compute model prediction when this term is 0
    } else if (attr(terms(lmerObj), "dataClasses")[f] == "factor") {
      predData[[f]] <- levels(lmerObj@frame[[f]])  
    } else {
      if (f %in% divide) {
        #divide into -1 SD, M, + 1 SD; or -2SD, -1SD, M, +1SD, +2SD
        fsd <- sd(lmerObj@frame[[f]], na.rm=TRUE)
        fm <- mean(lmerObj@frame[[f]], na.rm=TRUE)
        predData[[f]] <- if (n.divide==3) { c(fm-fsd, fm, fm+fsd)
            } else { c(fm-fsd*2, fm-fsd, fm, fm+fsd, fm+fsd*2) }
      } else {
        if (!is.null(names(n.cont))) {
          #Named vector specifying number of points to predict for each IV
          if (is.na(n.cont[f])) stop("Cannot locate number of continuous pred points for: ", f)
          predData[[f]] <- seq(min(lmerObj@frame[[f]], na.rm=TRUE), max(lmerObj@frame[[f]], na.rm=TRUE), length=n.cont[f])          
        } else {
          #treat as truly continuous predictor and compute models estimates across the range of observed values
          predData[[f]] <- seq(min(lmerObj@frame[[f]], na.rm=TRUE), max(lmerObj@frame[[f]], na.rm=TRUE), length=n.cont)
        }
      }
    }
  }
  
  #dependent variable
  dvname <- as.character(terms(lmerObj)[[2]])
  
  #populate the model-predicted estimates with 0s prior to building model matrix
  predData[[dvname]] <- 0
  
  #Develop a grid 
  predData <- do.call(expand.grid, list(predData))
  
  mm <- model.matrix(terms(lmerObj),predData)
  
  
  predData[[dvname]] <- mm %*% fixef(lmerObj)
  
  pvar1 <- diag(mm %*% tcrossprod(vcov(lmerObj),mm))
  tvar1 <- pvar1+VarCorr(lmerObj)[[1]][1] #assumes that the first element in VarCorr is subject
  
  #confidence and prediction intervals
  predData <- data.frame(
      predData, se=sqrt(pvar1),
      plo = predData[[dvname]]-2*sqrt(pvar1),
      phi = predData[[dvname]]+2*sqrt(pvar1),
      tlo = predData[[dvname]]-2*sqrt(tvar1),
      thi = predData[[dvname]]+2*sqrt(tvar1)
  )
  
  for (f in divide) {
    if (n.divide==3) { flevels <- c("-1 SD", "M", "+1 SD")
    } else if (n.divide==5) { flevels <- c("-2SD", "-1 SD", "M", "+1 SD", "+2 SD") }
    if (divide.prefix) flevels <- paste(Hmisc::capitalize(f), flevels, sep=": ")
    predData[[f]] <- factor(predData[[f]], levels=sort(unique(predData[[f]])), labels=flevels)
  }
  
  return(predData)
}


#leftovers -- junk
#    agexgeno.cont[comtEff] <- 1
#    agexgeno.cont[ageEff] <- 1

#two-way (should include all)

#three-way (should include all)
#    if (pfc_subcort_age %in% names(coefs)) agexgeno.cont[pfc_subcort_age] <- 1
#    if (pfc_subcort_comt %in% names(coefs)) agexgeno.cont[pfc_subcort_comt] <- 1

#before including age and genotype, check for subcort and pfc me and two-way (will exist except for ref levels) 
#if ME subcort is a term, then include it
#    if (meSubcort %in% names(coefs)) cont[meSubcort] <- 1
#    if (mePFC %in% names(coefs)) cont[mePFC] <- 1
#    if (pfc_subcort %in% names(coefs)) cont[pfc_subcort] <- 1


#allSubcort <- grep(paste("Subcortical", subcortLevels[l], sep=""), names(coefs), value=TRUE)
#allPFC <- grep(paste("PFC", pfcLevels[m], sep=""), names(coefs), value=TRUE)

#    if (pfc_age %in% names(coefs)) agexgeno.cont[pfc_age] <- 1
#    if (subcort_age %in% names(coefs)) agexgeno.cont[subcort_age] <- 1
#    if (pfc_comt %in% names(coefs)) agexgeno.cont[pfc_comt] <- 1
#    if (subcort_comt %in% names(coefs)) agexgeno.cont[subcort_comt] <- 1

