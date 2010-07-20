# CHNOSZ/util.character.R
# functions to work with character objects
# Copyright (C) 2006-2009 Jeffrey M. Dick

c2s <- function(x, sep=' ') {
  # make a string out of a character vector
  if(length(x) %in% c(0,1)) return(x)
  s <- x[1]
  for(i in 2:length(x)) s <- paste(s,x[i],sep=sep)
  return(s)
}

s2c <- function(x,sep=NULL,keep.sep=TRUE,n=NULL,move.sep=FALSE) {
  if(!is.null(sep)) {
    isep <- numeric()
    mysep <- character()
    for(i in 1:nchar(x)) {
    xnew <- x
      for(j in 1:length(sep)) {
        if(substr(x,i,i+nchar(sep[j])-1)==sep[j]) {
          if(!keep.sep) xnew <- c2s(c(substr(x,1,i-1),substr(x,i+nchar(sep[j]),nchar(x))),sep='')
          isep <- c(isep,i)
          mysep <- c(mysep,sep[j])
          break
        }
      }
    x <- xnew
    }
  } else isep <- 0:nchar(x)
  # the values of isep bracket the lengths of 
  # string we want to split
  if(length(isep)==0) {
    if(length(sys.calls())==1) cat(paste('s2c: no match for',sep,'in',x,'\n'))
    return(x)
  }
  if(!isep[1]==1 & !is.null(sep)) isep <- c(0,isep)
  if(isep[length(isep)] <= nchar(x) &!is.null(sep)) isep <- c(isep,nchar(x))
  v <- character()
  if(is.null(sep)) j <- 0 else {
    j <- -1
    isep[length(isep)] <- isep[length(isep)] + 1
  }
  # n: the number of elements we want to grab
  if(is.null(n)) n <- length(isep) else n <- n + 1
  for(i in 2:n) {
    v <- c(v,substr(x,isep[i-1]+1+j,isep[i]+j))
  }
  # alternatively move the separator to previous element
  if(move.sep & keep.sep & length(v) > 1) {
    for(i in 1:length(v)) {
      myv <- v[i]
      if(i!=length(v)) myv <- paste(myv,mysep[i],sep='')
      if(i!=1) myv <- substr(myv,start=length(mysep[i-1])+1,stop=nchar(myv))
      v[i] <- myv
    }
  }
  return(v)
}

can.be.numeric <- function(x) {
  # return FALSE if length of argument is zero
  if(length(x)==0) return(FALSE)
  if(length(x)>1) return(as.logical(sapply(x,can.be.numeric)))
  # don't warn about NAs in as.numeric
  oldopt <- options(warn=-1)
  cb <- FALSE
  if(!is.na(as.numeric(x))) cb <- TRUE
  if(x %in% c('.','+','-')) cb <- TRUE
  # let the user have their way
  options(oldopt)
  return(cb)
}

