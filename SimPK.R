
make_egDat <- function(vars,dat) {
  varNames <- sort(vars)
  F1 <- levels(as.factor(dat[[varNames[1]]]))
  F2 <- levels(as.factor(dat[[varNames[2]]]))
  nF1 <- length(F1)
  nF2 <- length(F2)
  egDat <- data.frame(var1=gl(nF1,nF2,labels=F1),var2=gl(nF2,1,length=nF1*nF2,labels=F2))
  names(egDat) <- varNames
  egDat$cellNames <- apply(egDat,1,paste0,collapse="")
  return(egDat)
}

simPK <- function(nsims,fp,fk,betaphat,betakhat,hessian,egDat) {
  varHat <- solve(hessian)
  require(mvtnorm)
  betaSim <- rmvnorm(nsims,mean=c(betaphat,betakhat),sigma=varHat)
  Xp <- model.matrix(fp,egDat)
  Xk <- model.matrix(fk,egDat)
  pSim <- alogit(betaSim[,1:NCOL(Xp)]%*%t(Xp))
  kSim <- alogit(betaSim[,(NCOL(Xp)+1):(NCOL(Xp)+NCOL(Xk))]%*%t(Xk))
  colnames(pSim) <- egDat$cellNames
  colnames(kSim) <- egDat$cellNames
  return(list(pSim=pSim,kSim=kSim))
}


