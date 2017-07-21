# The first 13 lines are the same as the previous ForDan.R.
load("trialData.Rdata")

# Keep vis, size, and search data for all bats, all sites, 
# easy, medium, or hard class.
# dat <- subset(trialData,subset=vis %in% c("E","M","D") & size %in% c("s","m"),
#               select=c(vis,size,grep("s\\d{1}",names(trialData),ignore.case=TRUE))) 
dat <- subset(trialData,
              subset=vis %in% c("E","M","D") & 
                spec %in% c("LABO","LACI","LANO","MYLU","PISU"),
              select=c(vis,spec,grep("s\\d{1}",names(trialData),
                                     ignore.case=TRUE))) 


scols <- grep("s\\d{1}",names(dat),
              ignore.case=TRUE) # Identify columns containing search outcome data

# Discard any rows with NA for first search.
dat <- dat[which(!is.na(dat[,scols[1]])),]

# Estimate p and k.
source("pkEst1.R")
fp <- ~vis*spec
fk <- ~vis
result <- pkEst1(fp,fk,dat)

# Kluge: need Hessian matrix. I altered my version of pkEst() to return the
# hessian. There's a function "optimHess()" that just approximates the
# Hessian, but it needs the objective function. In pkEst1.R, the objective
# function isn't a standalone function.
source("pkEst.R")
source("NegLogLike.R")
hessian <- pkEst(fp,fk,dat)$hessian

# Now simulate p's and k's for each "cell."
source("SimPK.R")
# egDat makes a data set with one row for each cell defined by vars
# and the levels of those variables in dat.
egDat <- make_egDat(vars=c("vis","spec"),dat)
sims <- simPK(nsims=15,fp,fk,result$betaphat,result$betakhat,hessian,egDat)

sims$pSim
sims$kSim
