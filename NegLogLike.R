logit <- function(x) log(x/(1-x))
alogit <- function(x) 1/(1+exp(-x))
require(matrixStats)

nLL <- function(theta,Xk,Xp,zeros,found){
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
#
#  # Version 1. About 10 seconds for the optimization with n=1426.
# nLL <- function(theta,Xk,Xp,nzero,found){
#   n <- length(nzero)
#   Betak <- matrix(theta[1:NCOL(Xk)],ncol=1)
#   k <- alogit(Xk%*%Betak)
#   Betap <- matrix(theta[(NCOL(Xk)+1):length(theta)],ncol=1)
#   p <- alogit(Xp%*%Betap)
#   # LL_1 is the component of the log-likelihood for the zeros.
#   LL_1 <- unlist(lapply(1:n,
#                         FUN=function(x) ifelse(nzero[x]>0,sum(log(1-p[x,1]*k[x]^(0:(nzero[x]-1)))),0)))
#   # LL_2 is the component of the log-likelihood for the ones.
#   LL_2 <- ifelse(found>0,log(p*k^(found-1)),0)
#   return(-sum(LL_1+LL_2))
# }
# 

