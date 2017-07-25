# Simulation experiment for estimator sum_i X_i/ghat_i with i=1,2 and ghat_i dependent.
# Account for variance in ghat with parametric bootstrap. Also account for binomial
# variation in X_i with parametric bootstrap.

# two turbines
M<-c(50, 100)
g<-c(0.3, 0.1)

gRho <- matrix(c(1,0.99,0.99,1),2,2) # rank correlation matrix between ghat[1] and ghat[2]
sqrt_gRho <- chol(gRho) # Need Cholesky decomposition for simulation.

# if g is fixed and known, then sum(X[i]/g[i]) gives what kind of CI?
nsim<-1000
nominal<-c(0.5, 0.8, 0.9, 0.95, 0.99) # nominal coverage
# if g is estimated from binomial trials with
N <- c(100, 100) # field trial carcasses for two classes
# beta distribution parameters for posterior g derived from field trials in the two classes
Bab <- array(dim=c(nsim, 2, 2)) # simi, Bab, class
Bab[,1,1]<-rbinom(nsim, N[1], g[1]) + 0.5
Bab[,1,2]<-rbinom(nsim, N[2], g[2]) + 0.5
Bab[,2,1]<-N[1]-Bab[,1,1]+1
Bab[,2,2]<-N[2]-Bab[,1,2]+1

ngsim <- 1000 # number of ghats to simulate to account for uncertainty in ghat.
# for each X, build confidence interval by sampling from g and calculating X[1]/ghat[1] + X[2]/ghat[2]
X<-t(array(rbinom(nsim * 2, M, g), dim=c(2, nsim))) # carcasses found in searches
CI<-array(dim=c(nsim, 2, length(nominal))) # simulated bootstrap confidence intervals for various confidence levels
Mhat<-array(dim=c(nsim,2)) # point estimates of M (median and mean of bootstapped x/g's)
for (yi in 1:nsim){
  #ghat<-cbind(rbeta(ngsim, Bab[yi, 1, 1], Bab[yi, 2, 1]), rbeta(ngsim, Bab[yi, 1, 2], Bab[yi, 2, 2]))
  # The next two lines simulate dependent ghat[1],ghat[2]
  U <- pnorm(matrix(rnorm(2*ngsim),ngsim,2)%*%sqrt_gRho) 
  ghat <- cbind(qbeta(U[,1],Bab[yi,1,1],Bab[yi,2,1]),qbeta(U[,2],Bab[yi,1,2],Bab[yi,2,2]))
  Xtilde <- cbind(rbinom(ngsim,size=round(X[yi,1]/ghat[,1]),prob=ghat[,1]),rbinom(ngsim,size=round(X[yi,2]/ghat[,2]),prob=ghat[,2]))
  Mhat[yi, 1]<-quantile(Xtilde[,1]/ghat[,1] + Xtilde[,2]/ghat[,2], prob=.5)
  Mhat[yi, 2]<-mean(Xtilde[,1]/ghat[,1] + Xtilde[,2]/ghat[,2])
  for (cvi in 1:length(nominal)) CI[yi, , cvi]<-quantile(Xtilde[,1]/ghat[,1] + Xtilde[,2]/ghat[,2], c((1-nominal[cvi])/2, 1-(1-nominal[cvi])/2))
}
for (cvi in 1:length(nominal)){
  coverage[cvi]<-sum(sum(M)>=CI[,1,cvi] & sum(M)<=CI[,2,cvi])/nsim #coverage probability = 0.7436 for nominal 90% CI
}
cbind(nominal, coverage)

# Coverages when Mhat=sum(X/ghat) (as opposed to sum(Xtilde/ghat))
#     nominal coverage
#[1,]    0.50    0.343
#[2,]    0.80    0.609
#[3,]    0.90    0.740
#[4,]    0.95    0.823
#[5,]    0.99    0.910


