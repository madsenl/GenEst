# Function to estimate p and k using data in dat with models given in model
# matrices Xp and Xk and search data in Data. Expecting Data to have column
# names s1, s2, etc.

pkEst <- function(fp,fk,dat) {
  
  Xp <- model.matrix(fp,dat)
  Xk <- model.matrix(fk,dat)
  
  scols <- grep("s\\d{1}",names(dat),
                ignore.case=TRUE) # Identify columns containing search outcome data
  Data <- dat[,scols]
  Data <- dat[,scols] # Extract the search data and convert the 0s to 1s.
  Data[which(is.na(Data) | Data==-999,arr.ind=TRUE)] <- 3 # Convert the NAs and -999s to 3s
  Data <- as.data.frame(t(apply(Data,1,sort))) # Sort each row.
  Data[which(Data==3,arr.ind=TRUE)] <- NA # Convert the 3s back to NAs.
  names(Data) <- paste0("s",1:20)
  
  # Create the zeros array. It has the same dimensions as Data, but has a 1
  # wherever Data has a 0 and NAs elsewhere. nLL() takes it as input.
  zeros <- as.matrix(Data)
  zeros[which(zeros==1,arr.ind=TRUE)] <- NA
  zeros[which(zeros==0,arr.ind=TRUE)] <- 1
  
  # Create the found vector. It has length equal to NROWS(Data). Each element
  # gives the search occasion when the carcass is found and is 0 if the 
  # carcass is never found.
  found <- apply(Data,1,FUN=function(x) match(1,x)) # In which search is carcass found?
  found[is.na(found)] <- 0 # match() returns NA if string not found. Need these to be numeric.
   
  # Find reasonable starting values. First calculate empp, a vector of length
  # NROW(Data). The ith element is the observed proportion of carcasses found
  # in search 1 for the ith carcass's group, as defined by the model matrix
  # Xp. According to the model, the probability that a carcass is found during
  # the first search is p.
  ngroups <- NCOL(Xp)
  empp <- numeric(nrow(Xp))
  groups <- findGroups(Xp)
  for (i in 1:ngroups) {
    empp[which(groups==i)] <- with(Data,mean(s1[which(groups==i & s1>=0)],na.rm=TRUE))
  }
  empp[which(empp==0)] <- 0.1 # Cells with all 0's set to 0.1
  empp[which(empp==1)] <- 0.9 # Cells with all 1's set to 0.1
  
  # Then calculate emppk, which is like empp, except based on search 2. 
  # According to the model, the probability that a carcass is found during
  # search 2 is p*k.
  ### This doesn't work if emppk and empp don't correspond to the same groups.
#   emppk <- numeric(nrow(Xk))
#   groups <- findGroups(Xk)
#   for (i in 1:length(groups)) {
#     emppk[which(groups==i)] <- with(Data,mean(s2[which(groups==i & s2>=0)],na.rm=TRUE))
#   }
  # Then get least squares estimates for the betas. The first NCOL(Xp)
  # elements of theta are starting values for the betas in the p model. The
  # remaining NCOL(Xk) elements are the betas for the k model.
#   theta <- c(solve(t(Xp)%*%Xp)%*%t(Xp)%*%logit(empp),
#              solve(t(Xk)%*%Xk)%*%t(Xk)%*%logit(emppk/empp))
  theta <- c(solve(t(Xp)%*%Xp)%*%t(Xp)%*%logit(empp),
             logit(rep(0.7,times=NCOL(Xk))))
  
  # Perform the optimization.
  result <- optim(par=theta,fn=nLL,Xp=Xp,Xk=Xk,zeros=zeros,found=found,method="BFGS",hessian=TRUE)
  
  # Extract the parameters estimates for the p and k models.
  betaphat <- result$par[1:NCOL(Xp)]
  betakhat <- result$par[(NCOL(Xp)+1):length(theta)]
  hessian <- optimHess(par=c(betaphat,betakhat),fn=nLL,Xp=Xp,Xk=Xk,zeros=zeros,found=found)
  return(list(betaphat=betaphat,
              betakhat=betakhat,
              aic=2*result$value + 2*length(result$par),
              convergence=result$convergence,
              hessian=hessian))
}

# Utility functions for logit and its inverse:
logit <- function(x) log(x/(1-x))
alogit <- function(x) 1/(1+exp(-x))

# Find groups defined by a design matrix X.
findGroups <- function(X) {
  nparams <- NCOL(X)
  Xrows <- sapply(1:NROW(X),FUN=function(x) paste0(X[x,],collapse=""))
  Urows <- unique(Xrows) # Should be NCOL(X) values.
  match(Xrows,Urows)
}