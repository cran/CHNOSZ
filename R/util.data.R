# CHNOSZ/util.data.R
# Copyright (C) 2007-2008 Jeffrey M. Dick
# add or change entries in the thermodynamic database

mod.obigt <- function(species,...,missingvalues=NA) {
  # add or modify species in thermo$obigt
  args <- list(...)
  for(i in 1:length(args)) {
    args[[i]] <- rep(args[[i]],length(species))
  }
  # a function to write dates in a specific format
  mydate <- function() {
    t <- date()
    tt <- s2c(t,sep=' ',keep.sep=FALSE)
    tday <- tt[3]
    tmonth <- tt[2]
    tyear <- substr(tt[5],start=3,stop=4)
    return(paste(tday,tmonth,tyear,sep='.'))
  }
  inew <- numeric()
  for(i in 1:length(species)) {
    is <- NULL
    sp <- species[i]
    if(is.factor(sp)) sp <- as.character(sp)
    if(!is.numeric(sp)) {
      if('state' %in% names(args)) ii <- info(sp,args$state,quiet=TRUE,return.approx=FALSE)
      else ii <- info(sp,quiet=TRUE,return.approx=FALSE)
    } else ii <- sp
    # 20090203 use return.approx in info calls above, so
    # the test is if length of return is zero
    #if(!is.na(ii) & !is.list(ii)) {
    if(length(ii) != 0) {
      #is <- which(species[i]==thermo$obigt$name)
      is <- ii
      if('state' %in% names(args)) mystate <- args$state[i] else mystate <- thermo$opt$state
      if(mystate %in% thermo$obigt$state[is]) is <- is[match(mystate,thermo$obigt$state[is])]
      else {if('state' %in% names(args)) is <- NULL else is <- is[1]}
    }
    if(!is.null(is)) {
      # to modify a species
      newrow <- thermo$obigt[is,]
      #return(newrow)
      for(j in 1:ncol(newrow)) {
        cnames <- s2c(colnames(newrow)[j],sep='.',keep.sep=FALSE)
        if(any(tolower(cnames) %in% tolower(names(args)))) {
          # use a provided value
          newrow[,j] <- args[[which(tolower(names(args)) %in% tolower(cnames))]][i]
        } else {
          # use default values
          if(any(cnames %in% c('name','formula'))) next
          #else if(is.na(missingvalues)) newrow[,j] <- missingvalues
          #else if(!missingvalues=='') newrow[,j] <- missingvalues
          if(!missing(missingvalues) & !any(cnames %in% c('state','source1','source2'))) 
            newrow[,j] <- missingvalues
          if(any(cnames %in% 'source1')) newrow[,j] <- 'USER'
          if(any(cnames %in% 'source2')) newrow[,j] <- NA
        }
      }
      newrow$date <- mydate()
      r1 <- as.character(newrow)
      r2 <- as.character(thermo$obigt[is,])
      r2[is.na(r2)] <- 'NA'
      r1[is.na(r1)] <- 'NA'
      if(!identical(r1,r2)) {
        cat(paste('mod.obigt: updating ',newrow$name,' ',newrow$state,' (',is,').\n',sep=''))
      } else cat(paste('mod.obigt: no change for ',newrow$name,' ',newrow$state,' (',is,').\n',sep=''))
      thermo$obigt[is,] <<- newrow
      inew <- c(inew,is)
    } else {
      # add a species
      newrow <- thermo$obigt[1,]
      for(j in 1:ncol(newrow)) {
        cnames <- s2c(colnames(newrow)[j],sep='.',keep.sep=FALSE)
        if(any(cnames %in% names(args))) {
          # use a provided value
          newrow[,j] <- args[[which(names(args) %in% cnames)]][i]
        } else {
          # use default values
          if(any(cnames %in% c('name'))) newrow[,j] <- sp
          else if(any(cnames=='date')) newrow[,j] <- mydate()
          else if(is.na(missingvalues) & colnames(newrow)[j]=='state')
            newrow[,j] <- thermo$opt$state
          else if(!is.na(missingvalues)) newrow[,j] <- missingvalues
          #else if(!any(cnames=='formula')) newrow[,j] <- missingvalues
          else if(any(cnames=='source1')) newrow[,j] <- 'USER'
          else newrow[,j] <- NA
        }
      }
      if(is.na(newrow$formula)) warning(paste('mod.obigt: formula of ',newrow$name,
        ' ',newrow$state,' is NA.',sep=''),call.=FALSE)
      cat(paste('mod.obigt: adding ',newrow$name,' ',newrow$state,' (',nrow(thermo$obigt)+1,').\n',sep=''))
      thermo$obigt <<- rbind(thermo$obigt,newrow)
      inew <- c(inew,nrow(thermo$obigt))
    }
  }
  return(invisible(inew))
}

change <- function(name,...) {
  # a wrapper for mod.obigt and mod.buffer
  if(substr(name[1],1,1)=='_') {
    name <- substr(name,2,nchar(name))
    return(mod.buffer(name,...))
  } else {
    return(mod.obigt(species=name,...))
  }
}

add.obigt <- function(file="obigt.csv",force=FALSE) {
  # add/replace entries in thermo$obigt
  # from values saved in a file
  # only replace if force=TRUE
  to1 <- thermo$obigt
  id1 <- paste(to1$name,to1$state)
  to2 <- read.csv(file,as.is=TRUE)
  id2 <- paste(to2$name,to2$state)
  # check if the file is compatible with thermo$obigt
  tr <- try(rbind(to1,to2),silent=TRUE)
  if(identical(class(tr),'try-error')) stop(paste(file,"is not compatible with thermo$obigt data table."))
  # identify duplicated
  idup1 <- which(id1 %in% id2)
  idup2 <- which(id2 %in% id1)
  ndup <- length(idup2)
  nnew <- nrow(to2) - ndup
  iadd <- 1:nrow(to2)
  if(force) {
    # drop entries from original
    if(length(idup1) > 0) to1 <- to1[-idup1,]
  } else {
    if(length(idup2) > 0) iadd <- iadd[-idup2]
    ndup <- 0
  }
  inew <- numeric()
  if(length(iadd) > 0) {
    inew <- nrow(to1) + 1:length(iadd)
    to1 <- rbind(to1,to2[iadd,])
  }
  thermo$obigt <<- to1
  cat(paste("add.obigt: added",length(iadd),"of",nrow(to2),"species from",file,"(",ndup,"replacements,",nnew,"new )\n"))
  return(invisible(inew))
}

danger <- function() {
  cat("danger: loading supplemental thermodynamic data file\n")
  cat("danger: use with care; some properties are not thermodynamically\n")
  cat("danger: consistent with the default database\n")
  myfile <- system.file("data/OBIGT-2.csv",package="CHNOSZ")
  i <- add.obigt(myfile,force=TRUE)
  cat("danger: to undo, reinitialize data object with 'data(thermo)'\n")
  return(invisible(i))
}
