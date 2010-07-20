# CHNOSZ/revisit.R
# Copyright (C) 2009 Jeffrey M. Dick
# functions related to diversity calculations
# 20090415 jmd

revisit <- function(logact,which="cv",as.is=FALSE) {
  # given a list of logarithms of activities of species
  # (as matrices or vectors) calculate a diversity index
  # of the same dimensions
  # calculate the Shannon entropy "shannon"
  # or the standard deviation "sd"
  # or the coefficient of variation "cv"
  # ... of the activities, not their logarithms!
  # (unless as.is is TRUE)
  act <- logact
  # vectorize and get activities
  mydim <- dim(logact[[1]])
  for(i in 1:length(logact)) {
    logact[[i]] <- as.vector(logact[[i]])
    if(as.is) act[[i]] <- as.vector(logact[[i]])
    else act[[i]] <- as.vector(10^logact[[i]])
    if(i==1) acttotal <- act[[i]] else acttotal <- acttotal + act[[i]]
  }
  # make a place for the results
  H <- logact[[1]]
  H[] <- 0
  if(which=="shannon") {
    # do the calculation
    for(i in 1:length(act)) {
      dH <- -act[[i]]/acttotal*log(act[[i]]/acttotal)
      if(!is.na(dH)) H <- H + dH
    }
  } else if(which %in% c("sd","cv")) {
    # build a matrix; rows are species
    myad <- matrix(,length(logact),length(logact[[1]]))
    for(j in 1:length(logact)) myad[j,] <- act[[j]]
    # get the standard deviations
    H <- sd(myad)
    # for coefficient of variation, divide by the mean
    if(which=="cv") {
      myad <- as.data.frame(myad)
      H <- H / mean(myad)
    }
  }
  # replace dims
  dim(H) <- mydim
  return(H)
}

richness <- function(logact,logactmin=-4,mult=1,as.is=FALSE) {
  # given a list of logarithms of activities of species
  # (as matrices or vectors) calculate the richness 
  # make a place for the results
  H <- logact[[1]]
  H[] <- 0
  # loop over species
  # scale the values if requested
  mult <- rep(mult,length.out=length(logact))
  if(is.numeric(logactmin)) logactmin <- rep(logactmin,length.out=length(logact))
  else {
    if(as.is) logactmin <- sapply(1:length(logact), function(x) {mean(as.numeric(logact[[x]]))})
    else logactmin <- sapply(1:length(logact), function(x) {log10(mean(10^as.numeric(logact[[x]])))})
  }
  for(j in 1:length(logact)) {
    isthere <- logact[[j]] > logactmin[j]
    H[isthere] <- H[isthere] + mult[j] 
  }
  return(H)
}

where.extreme <- function(z,what="minimum") {
  # takes a matrix, returns the x,y coordinates of the extremum
  if(what %in% c("minimum","min")) iext <- which.min(z)
  else iext <- which.max(z)
  # in case there are multiple instances of the extremum
  # (esp. for richness calculations)
  iext <- which(z==z[iext])
  x.out <- y.out <- numeric()
  xres <- ncol(z)
  yres <- nrow(z)
  for(i in 1:length(iext)) {
    # column (x coord)
    x <- ceiling(iext[i]/xres)
    # and row (y coord)
    y <- iext[i] - floor(iext[i]/yres)*yres
    if(y==0) y <- yres  # there's a more eloquent way...
    x.out <- c(x.out,x)
    y.out <- c(y.out,y)
  }
  return(list(x=y.out,y=x.out))
}

extremes <- function(z,what="minimum") {
  # takes a matrix, returns the y as f(x) and x as f(y)
  # trajectories of the extreme
  y <- x <- numeric()
  xres <- ncol(z)
  yres <- nrow(z)
  if(what %in% c("minimum","min")) {
    for(i in 1:xres) y <- c(y,which.min(z[i,]))
    for(i in 1:yres) x <- c(x,which.min(z[,i]))
  } else {
    for(i in 1:xres) y <- c(y,which.max(z[i,]))
    for(i in 1:yres) x <- c(x,which.max(z[,i]))
  }
  return(list(x=x,y=y))
}


draw.diversity <- function(d,which="cv",logactmin=-4,
  col=par("fg"),as.is=FALSE,yline=2,ylim=NULL,ispecies=NULL,add=FALSE,
  cex=par("cex"),lwd=par("lwd"),mar=par("mar"),mult=1,side=1:4) {
  # plot diversity indices of relative abundances
  # 20090316 jmd
  if(is.null(ispecies)) ispecies <- 1:length(d$logact)
  # take a subset (or all) of the species
  logact <- d$logact[ispecies]
  # do diversity calculations
  if(which=="richness") H <- richness(logact,logactmin=logactmin,mult=mult,as.is=as.is)
  else H <- revisit(logact,which,as.is=as.is)
  # what are the values of the variables
  xname <- d$xname
  yname <- d$yname
  xres <- d$xlim[3]
  xlim <- d$xlim[1:2]
  #if(xname=="T") xlim <- outvert(xlim,"K")
  if(xname=="H+") {
    xname <- "pH"
    xlim <- -xlim 
  }
  xs <- seq(xlim[1],xlim[2],length.out=xres)
  # make the plot
  if(yname=="") {
    # a 1-D plot
    if(is.null(ylim)) {
      if(which=="richness") {
          ylim <- 0
          if(max(H) > ylim) ylim <- max(H) + 1
          ylim <- c(0,ylim)
      }
      else ylim <- range(H[[1]])
    }
    if(!add) thermo.plot.new(xlim=xlim,ylim=ylim,xlab=axis.label(xname),ylab=which,yline=yline,
      cex=cex,lwd=lwd,mar=mar,side=side)
    lines(xs,as.vector(H),col=col)
    return(invisible(list(xs=xs,H=H)))
  } else {
    # a 2-D plot
    yres <- d$ylim[3]
    ylim <- d$ylim[1:2]
    #if(yname=="T") ylim <- outvert(ylim,"K")
    if(yname=="H+") {
      yname <- "pH"
      ylim <- -ylim 
    }
    ys <- seq(ylim[1],ylim[2],length.out=yres)
    if(which %in% c("richness","shannon")) myext <- "maximum" else myext <- "minimum"
    # start the plot
    if(!add) thermo.plot.new(xlim=xlim,ylim=ylim,xlab=axis.label(xname),ylab=axis.label(yname),yline=yline,side=side)
    contour(xs,ys,H,add=TRUE)
    iext <- where.extreme(H,myext)
    # plot the location(s) of the extremum
    points(xs[iext$x],ys[iext$y],pch=8,cex=2)
    # show trajectories of the extrema
    iexts <- extremes(H,myext)
    # take out large jumps
    yext <- ys[iexts$y]
    yext.1 <- c(yext[2:length(yext)],yext[length(yext)])
    yext.2 <- c(yext[1],yext[1:length(yext)-1])
    yext[abs(yext.1-yext)/abs(diff(range(ys))) > 0.1] <- NA
    yext[abs(yext.2-yext)/abs(diff(range(ys))) > 0.1] <- NA
    lines(xs,yext,lty=3,col="blue")
    xext <- xs[iexts$x]
    xext.1 <- c(xext[2:length(xext)],xext[length(xext)])
    xext.2 <- c(xext[1],xext[1:length(xext)-1])
    xext[abs(xext.1-xext)/abs(diff(range(xs))) > 0.1] <- NA
    xext[abs(xext.2-xext)/abs(diff(range(xs))) > 0.1] <- NA
    lines(xext,ys,lty=3,col="seagreen")
    # what is the extreme value
    extval <- H[iext$y,iext$x]
    return(invisible(list(H=H,ixmin=iext$x,iymin=iext$y,
      xmin=xs[iext$x],ymin=ys[iext$y],extval=extval)))
  }
}
