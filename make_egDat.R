# This function creates a data set with columns given in vars. There is a row for
# every unique combination of the levels of vars found in dat. A final column
# gives a unique code which is the concatenation of the labels in the other columns.
# The rows define the unique cells given by combinations of categories in vars.
# If you want only one group, do: make_egDat(NULL,dat)
make_egDat <- function(vars,dat) {
  nvars <- length(vars) # Should be 0, 1, or 2 for now, but function will work for more.
  if (nvars==0) {
    return(data.frame(group="all",cellNames="all"))
  } else {
      if(any(is.na(match(vars,names(dat))))) {
        stop("At least one var not found in dat.")
      }
      varNames <- sort(vars)
      varLabels <- list()
      varNlevels <- list()
      for (i in 1:nvars) {
        varLabels[[i]] <- levels(as.factor(dat[[varNames[i]]]))
        varNlevels[[i]] <- length(varLabels[[i]])
      }
      reps <- cumprod(varNlevels)[nvars:1] # Reverse cumulative product
      egDat <- data.frame(var1=gl(varNlevels[[1]],1,length=reps[1],labels=varLabels[[1]]))
      if (nvars > 1) {
        for (j in 2:nvars) {
          egDat[[paste("var",j,collapse=NULL)]] <- gl(varNlevels[[j]],reps[j],length=reps[1],labels=varLabels[[j]])
        }
      }  
    }
  names(egDat) <- varNames
  egDat$cellNames <- apply(egDat,1,paste0,collapse="")
  return(egDat)
}