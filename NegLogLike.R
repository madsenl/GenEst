logit <- function(x) log(x/(1-x))
alogit <- function(x) 1/(1+exp(-x))
require(matrixStats)

nLL <- function(theta,Xp,Xk,zeros,found){
  Betap <- matrix(theta[1:NCOL(Xp)],ncol=1)
  p <- alogit(Xp%*%Betap)
  Betak <- matrix(theta[(NCOL(Xp)+1):length(theta)],ncol=1)
  k <- alogit(Xk%*%Betak)

  # LL_1 is the component of the log-likelihood for the zeros.
  LL_1 <- rowSums2(log(1-sweep(outer(k[,1],0:19,"^")*zeros,1,p,"*")),na.rm=TRUE)  
  # LL_2 is the component of the log-likelihood for the ones.
  LL_2 <- ifelse(found>0,log(p*k^(found-1)),0)
  return(-sum(LL_1+LL_2))
}
