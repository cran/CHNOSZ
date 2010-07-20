# CHNOSZ/util.list.R
# functions to work with lists
# Copyright (C) 2008-2009 Jeffrey M. Dick

which.pmax <- function (elts, na.rm = FALSE, pmin=FALSE) {
  # adapted from R's pmax. elts is a list of numeric vectors
  if(!is.numeric(elts[[1]])[1]) {
    if(is.list(elts[[1]])) elts[[1]] <- elts[[1]][[1]]
    else elts[[1]] <- as.numeric(elts[[1]])
  }
  mmm <- as.vector(elts[[1]])
  which.mmm <- rep(1,length(elts[[1]]))
  has.na <- FALSE
  for (i in 2:length(elts)) {
    if(!is.numeric(elts[[i]])[1]) {
      if(is.list(elts[[i]])) elts[[i]] <- elts[[i]][[1]]
      else elts[[i]] <- as.numeric(elts[[i]])
    }
    work <- cbind(mmm, as.vector(elts[[i]]))
    nas <- is.na(work)
    if (has.na || (has.na <- any(nas))) {
      work[, 1][nas[, 1]] <- work[, 2][nas[, 1]]
      work[, 2][nas[, 2]] <- work[, 1][nas[, 2]]
    }
    if(pmin) change <- work[, 1] > work[, 2]
    else change <- work[, 1] < work[, 2]
    change <- change & !is.na(change)
    work[, 1][change] <- work[, 2][change]
    which.mmm[change] <- i
    if (has.na && !na.rm) {
      work[, 1][nas[, 1] | nas[, 2]] <- NA
      which.mmm[nas[, 1] | nas[, 2]] <- NA
    }
    mmm <- work[, 1]
  }
  mostattributes(mmm) <- attributes(elts[[1]])
  which.mmm
}

lsub <- function(x,y) {
  # subtract elements of list y from those of
  # list x to give list z
  z <- x
  for(i in 1:length(x)) z[[i]] <- x[[i]] - y[[i]]
  return(z)
}

lsum <- function(x,y) {
  # sum up the respective elements of lists
  # x and y to give list z of the same length
  z <- x
  for(i in 1:length(x)) z[[i]] <- x[[i]] + y[[i]]
  return(z)
}

pprod <- function(x,y) {
  # multiply each element of vector y
  # by corresponding value in list x
  pfun <- function(i) x[[i]]*y[i]
  lapply(1:length(y),pfun)
}

psum <- function(x) {
  # sum all elements of a list
  s <- x[[1]]
  for(j in 2:length(x)) s <- s+x[[j]]
  return(s)
}

mylapply <- function(X,FUN,...) {
  # wrapper to run lapply or mclapply
  if(length(X) > 20) {
    if("multicore" %in% (.packages())) do.call("mclapply",list(X,FUN,...))
    else lapply(X,FUN,...)
  } else lapply(X,FUN,...)
}

