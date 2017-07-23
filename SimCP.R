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
# survOut <- survival::survreg(surv ~ spec+vis, dist = "exp", data = CP)
# survOut <- survival::survreg(surv ~ spec*vis, dist = "weibull", data = CP)
# 
# egDat <- make_egDat(vars=all.vars(formula(survOut))[-1],dat=CP)
# sims <- simCP(20,survOut,egDat)
# sims$pdbSim
# sims$pdaSim

simCP <- function(nsims,survOut,egDat) {
  if (survOut$dist=="exponential") {
    betahat <- survOut$coefficients # Exponential model's scale is 1.
  } else {
    betahat <- c(survOut$coefficients,log(survOut$scale)) 
  }
  betaSim <- mvtnorm::rmvnorm(nsims,mean=betahat,sigma=survOut$var)
  Xcp <- model.matrix(as.formula(paste("~",as.character(formula(survOut))[3],collapse=NULL)),
                       data=egDat) # Extract formula and create model matrix.
  # eoa parameters (pda, pdb) <-> survival (1/scale, exp(coef[i]))
  pdbSim <- exp(betaSim[,1:NCOL(Xcp)]%*%t(Xcp))
  if (survOut$dist=="exponential") {
    pdaSim <- matrix(1,nrow=nsims,ncol=NROW(Xcp))
  } else {   
    pdaSim <- matrix(1/exp(betaSim[,NCOL(betaSim)]),nrow=nsims,ncol=NROW(Xcp))
  }  
  colnames(pdbSim) <- egDat$cellNames
  colnames(pdaSim) <- egDat$cellNames
  return(list(pdbSim=pdbSim,pdaSim=pdaSim))
}


