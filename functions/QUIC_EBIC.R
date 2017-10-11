testts <- read.table(file.path("/Users/nth7/Downloads/roits.txt"))
covt <- cov(testts)
nsamples <- nrow(testts)
npath=100
gammaRng = seq(0, 1, 0.1)

library(QUIC)

#Code by Ben Cassidy ported from MATLAB
#helper function for calculating EBIC of a candidate inverse covariance matrix (concentration)
EBICcalc <- function(concentration, samplecov, nsamples, gamma) {
  require(pracma)
  #subfunction for computing loglikelihood
  logLik <- function(concentration, samplecov) { logDetMat(concentration) - sum(diag(concentration * samplecov)) }
  
  E = nnz(tril(concentration,-1));
  p = size(samplecov,1);
  
  L = logLik(concentration, samplecov)
  
  dimpen_samp = E*log(nsamples)/nsamples
  dimpen_EBIC = 4*E*gamma*log(p)/nsamples
  
  EBIC = - L + dimpen_samp + dimpen_EBIC;
  
  #returns [EBIC, L, dimpen_samp, dimpen_EBIC]
  return(list(EBIC=EBIC, L=L, dimpen_samp=dimpen_samp, dimpen_EBIC=dimpen_EBIC))
}

#log determinant of symmetric matrix A
logDetMat <- function(A) {
  # A MUST BE A SYMMETRIC MATRIX
  
  # BETTER way to calculate log(det(A)): goo.gl/VUe8S or goo.gl/bYnD3
  if (all(eigen(A)$values > 0)) { #positive definite, so use chol. method
    cholCov = chol(A)
    logdet = 2*sum(log(diag(cholCov)))
  } else {
    # otherwise can always use LU factorisation: http://goo.gl/MM5TG
    lufac = lu(A) #from Matrix package #matlab returns: [~, Ue, Pe]
    lumat = expand(lufac) #returns a list with L, U, P
    
    du = diag(lumat$U)
    ce = det(lumat$P) * prod(sign(du))
    logdet = Re(log(ce) + sum(log(abs(du))))
  }
  
  return(logdet)
}

# S is empirical covariance, cov(data); but we might only have cov already
# data is nsamples x p, so pxp covariance matrix
# gamma is EBIC tuning. 0 <= \gam <=1
# npath is the number of candidate solutions to estimate, along the regularization path.

QUIC_EBIC <- function (S, nsamples, npath, gammaRng) {
  #gammaRng is a vector of gamma values to test
  require(QUIC)
  require(Matrix)
  require(pracma) #has a bunch of matlab functions (e.g., logspace) to save us porting things
  maxReg = max(abs(S)) #upper regularization parameter
  
  # regpath = (maxReg:-(maxReg / npath):(maxReg / npath));   % fractions of L, along reg path
  regpath = sort( logspace(log10(maxReg*1e-5),log10(maxReg),npath) , decreasing=TRUE)
  # regpath = linspace(maxReg*1e-2, maxReg, npath)
  regpath <- regpath[1:15] #quick test
  #tol = 1e-15
  tol = 1e-6 #only for testing (fast)
  msg = 1 # or 2 
  #maxIter = 1e3
  maxIter = 1e3
  
  #Estimate QUIC along the regularization path defined above
  system.time(quic_est <- QUIC(S=S, rho=1.0, path=regpath, tol=tol, msg=msg, maxIter=maxIter))
  
  icov_est_QUIC <- quic_est$X
  WP <- quic_est$W
  optP <- quic_est$opt
  cputimeP <- quic_est$cputime
  iterP <- quic_est$iter
  dGapP <- quic_est$dGap
  
  #glasso comparison
  #[icov_est_GLASSO, cov_est, errflag, regpath, cmdout]=GLASSO_calc(S, sort(regpath));
  
  EBIC_QUIC = matrix(NA, nrow=npath, ncol=length(gammaRng))
  logLik = matrix(NA, nrow=npath, ncol=length(gammaRng))
  dimpenalty.samp = matrix(NA, nrow=npath, ncol=length(gammaRng))
  dimpenalty.EBIC = matrix(NA, nrow=npath, ncol=length(gammaRng))
  
  #compute EBIC for each regularization parameter rho (along path)
  for (gamIdx in 1:length(gammaRng)) {
    gamma = gammaRng[gamIdx]
    
    for (lp in 1:(dim(icov_est_QUIC)[3])) {
      concentration_QUIC = icov_est_QUIC[,,lp]
      eb <- EBICcalc(concentration_QUIC, S, nsamples, gamma)
      EBIC_QUIC[lp, gamIdx] <- eb$EBIC
      logLik[lp, gamIdx] <- eb$L
      dimpenalty.samp[lp, gamIdx] <- eb$dimpen_samp
      dimpenalty.EBIC[lp, gamIdx] <- eb$dimpen_EBIC
    }

  }
  
  #% [~, EBIC_opt_ind] = min(EBIC_QUIC(:,end)); %% USING EBIC == 1 ! 
    
  #convert partial covariance matrix to partial correlation metric
  parcor_QUIC = array(NA, dim=dim(icov_est_QUIC))
  for (lp in 1:(dim(icov_est_QUIC)[3])) {
    icov_est_buf = icov_est_QUIC[,,lp] #node x node
    icovbufDiag = 1/sqrt(diag(icov_est_buf))
    parcor_buf = - diag(icovbufDiag) %*% icov_est_buf %*% diag(icovbufDiag)
    #parcor_buf[eye(size(parcor_buf)) > 0] = 1
    diag(parcor_buf) = 1
    parcor_QUIC[,,lp] = parcor_buf
  }
  
  ret <- list()
  ret$oPARCOR <- parcor_QUIC  #oPARCOR.parcor = parcor_QUIC
  ret$oICOV  <- list(icov=icov_est_QUIC, cputimeP=cputimeP, optP=optP, iterP = iterP, dGapP <- dGapP)
  ret$oEBIC <- EBIC_QUIC
  ret$regpath <- regpath
  ret$maxReg <- maxReg
  ret$logLik <- logLik
  #% o.logLik = logLik;
  #% o.dimpenalty = dimpenalty;
  #o.o.WP = WP;
  
  return(ret)
  
  #% o.parcor.gamEq1_toUse = parcor_QUIC(:,:,EBIC_opt_ind); % this may get trumped afterwards, if we choose another EBIC val; but this is default.
  #% o.icov.gamEq1_toUse = icov_est_QUIC(:,:,EBIC_opt_ind);
  
  #ret$parcor.gamEq1_toUse <- parcor_QUIC[,,]
}

#also see https://cran.r-project.org/web/packages/stabs/vignettes/stabs_graphs.html
#for info about selection of rho in path mode
result <- QUIC_EBIC(covt, nsamples, npath, gammaRng)


