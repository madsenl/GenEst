# Function to simulate nsims realizations of a vector of dependent ghats
# with beta(a,b) marginals. Beta parameters a and b should have one element
# for each class. grho is rank correlation between pairs of ghats.
# Returns g_hat which is nsims by length(a).
ghatSim <- function(nsims,a,b,grho) {
  nclasses <- length(a)
  Rho <- matrix(grho,nclasses,nclasses)
  diag(Rho) <- 1
  sqrtRho <- NA
  try(sqrtRho <- chol(Rho), TRUE)
  if (is.na(sqrtRho[1])) {
    stop("Correlation matrix is not positive definite.")
  }
  U <- pnorm(matrix(rnorm(nclasses*nsims),nsims,nclasses)%*%sqrtRho) # Rho^(1/2) * standard normals transformed to uniform(0,1) so that cor(t(Z)) gives Rho approximately.
  g_hat <- sapply(1:nclasses,FUN=function (x) qbeta(Z[x,],a[x],b[x])) 
  return(g_hat)
}
