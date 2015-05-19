# CHNOSZ/util.program.R
# various programming-related functions

caller.name <- function(n=2) {
  # returns the name of the calling function n frames up
  # (n=2: the caller of the function that calls this one)
  # or character() if called interactively
  if(sys.nframe() < n) name <- character()
  else {
    sc <- sys.call(-n)[[1]]
    name <- try(as.character(sc),silent=TRUE)
    # also return character() if the value from sys.call is
    # the function itself (why does this sometimes happen,
    # e.g. when called from affinity()?)
    if(class(name)=="try-error") name <- character()
  }
  return(name)
}

palply <- function(X, FUN, ...) {
  # a wrapper function to run parLapply if length(X) > 10000
  # and package 'parallel' is available, otherwise run lapply
  if(length(X) > 10000 & "parallel" %in% (.packages())) {
    # the calculations - modified from ?parLapply
    ## Use option mc.cores to choose an appropriate cluster size.
    # or detectCores if that is NULL, and set max at 2 for now
    # (to be nice to CRAN etc.)
    nCores <- max(getOption("mc.cores", parallel::detectCores()), 2)
    # don't load methods package, to save startup time - ?makeCluster
    cl <- parallel::makeCluster(nCores, methods=FALSE)
    out <- parallel::parLapply(cl, X, FUN, ...)
    parallel::stopCluster(cl)
  } else out <- lapply(X, FUN, ...)
  return(out)
}
