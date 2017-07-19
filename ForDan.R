load("trialData.Rdata")

# Keep vis, size, and search data for all bats, all sites, 
# easy, medium, or hard class.
dat <- subset(trialData,subset=vis %in% c("E","M","D") & size %in% c("s","m"),
              select=c(vis,size,grep("s\\d{1}",names(trialData),ignore.case=TRUE))) 

# Identify columns containing search outcome data
scols <- grep("s\\d{1}",names(dat),ignore.case=TRUE) 
sdata <- dat[,scols] # Extract the search data
sdata[which(is.na(sdata) | sdata==-999,arr.ind=TRUE)] <- 3 # Convert the NAs and -999s to 3s
sdata <- as.data.frame(t(apply(sdata,1,sort))) # Sort each row.
sdata[which(sdata==3,arr.ind=TRUE)] <- NA # Convert the 3s back to NAs.
names(sdata) <- paste0("s",1:20)
dat[,scols] <- sdata 
# The search data should now have no -999s.

# A couple of example calls to the fitting function:
pkEst(~vis*size,~vis,dat)
pkEst(~vis,~vis,dat) # This is the second model shown in PKEst.pdf.