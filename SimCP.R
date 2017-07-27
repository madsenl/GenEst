# Function to simulate carcass persistence (CP) parameters. Input is
# nsims = number of simulations, survOut = a survreg object, and
# egDat = output from make_egDat().
# Output is a list containing pdb and pda. These are matrices with nsims rows
# and number of columns equal to the number of "cells" defined by the (at most)
# two categorical variables used in the model.
# Each column of pdb contains nsims values simulated from the distribution of 
# pdb for the corresponding cell. The columns of the pda matrix are identical
# since that parameter is the same for all cells. If the survival distribution
# is exponential, then the pda's are all 1's.
#
# A couple examples: The CP data set is from CPmod.R.
# survOut <- survival::survreg(surv ~ spec*vis, dist = "exp", data = CP)
# survOut <- survival::survreg(surv ~ spec*vis, dist = "weibull", data = CP)
# survOut <- survival::survreg(surv ~ spec*vis, dist = "lognormal", data = CP)
# 
# egDat <- make_egDat(vars=all.vars(formula(survOut))[-1],dat=CP)
# sims <- simCP(20,survOut,egDat)
# sims$pdbSim
# sims$pdaSim

simCP <- function(nsims,survOut,egDat) {
  Xcp <- model.matrix(as.formula(paste("~",as.character(formula(survOut))[3],collapse=NULL)),
                      data=egDat) # Extract formula and create model matrix.
  if (survOut$dist=="exponential") {
    # Exponential model's scale is 1.
    thetahat <- survOut$coefficients 
  } else {
    # For Weibull, loglogistic, and lognormal, the parameter vector consists of
    # the coefficients of the predictor variables followed by log(scale).
    thetahat <- c(survOut$coefficients,log(survOut$scale)) 
  }
  thetaSim <- mvtnorm::rmvnorm(nsims,mean=thetahat,sigma=survOut$var)
  # Models and reparametrizations.
  # exponential: a = 1; b = meanCP = exp(mod.e$coef)
  # Weibull: a = shape = 1/mod.w$scale; b = scale = exp(mod.w$coef*Xcp)
  # log-logistic: a = shape = 1/mod.ll$scale; b = scale = exp(mod.ll$*coef*Xcp)
  # lognormal: a = sdlog^2 = mod.ln$scale^2; b = meanlog = mod.lm$coef*Xcp
  pdaSim <- switch(survOut$dist,
                   exponential=matrix(1,nrow=nsims,ncol=NROW(Xcp)),
                   weibull=matrix(1/exp(thetaSim[,NCOL(thetaSim)]),nrow=nsims,ncol=NROW(Xcp)),
                   loglogistic=matrix(1/exp(thetaSim[,NCOL(thetaSim)]),nrow=nsims,ncol=NROW(Xcp)),
                   lognormal=matrix(1/thetaSim[,NCOL(thetaSim)]^2,nrow=nsims,ncol=NROW(Xcp)))
  pdbSim <- switch(survOut$dist, 
                   exponential=exp(thetaSim[,1:NCOL(Xcp)]%*%t(Xcp)),
                   weibull=exp(thetaSim[,1:NCOL(Xcp)]%*%t(Xcp)),
                   loglogistic=exp(thetaSim[,1:NCOL(Xcp)]%*%t(Xcp)),
                   lognormal=thetaSim[,1:NCOL(Xcp)]%*%t(Xcp))
  colnames(pdbSim) <- egDat$cellNames
  colnames(pdaSim) <- egDat$cellNames
  return(list(pdaSim=pdaSim,pdbSim=pdbSim))
}