# CHNOSZ/revisit.R
# 20090415 functions related to diversity calculations
# 20100929 merged draw.diversity and revisit

where.extreme <- function(z,target,do.sat=FALSE) {
  if(missing(target)) stop("no target specified")
  # are we interested in a maximum or minimum?
  if(tolower(target) %in% c("sd","sd.log","cv","cv.log","rmsd","cvrmsd")) myext <- "minimum" else myext <- "maximum"
  # do we care about the sign of the index?
  if(tolower(target) %in% c("sd","sd.log","cv","cv.log","rmsd","cvrmsd")) doabs <- TRUE else doabs <- FALSE
  # takes a matrix, returns the x,y coordinates of the extremum
  if(doabs) z <- abs(z)
  if(myext %in% c("minimum","min")) iext <- which.min(z)
  else iext <- which.max(z)
  ret.val <- iext
  # for a matrix, it gets more complicated esp.
  # if there are multiple instances of the extremum
  # (happens with saturated richnesses)
  if(do.sat & length(dim(z))==2) {
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
    ret.val <- list(ix=y.out,iy=x.out)
  }
  return(ret.val)
}

extremes <- function(z,target) {
  if(missing(target)) stop("no target specified")
  # are we interested in a maximum or minimum?
  if(tolower(target) %in% c("sd","sd.log","cv","cv.log","rmsd","cvrmsd")) myext <- "minimum" else myext <- "maximum"
  # do we care about the sign of the index?
  if(tolower(target) %in% c("sd","sd.log","cv","cv.log","rmsd","cvrmsd")) doabs <- TRUE else doabs <- FALSE
  # takes a matrix, returns the y as f(x) and x as f(y)
  # trajectories of the extreme
  if(doabs) z <- abs(z)
  y <- x <- numeric()
  xres <- ncol(z)
  yres <- nrow(z)
  if(myext %in% c("minimum","min")) {
    for(i in 1:xres) y <- c(y,which.min(z[i,]))
    for(i in 1:yres) x <- c(x,which.min(z[,i]))
  } else {
    for(i in 1:xres) y <- c(y,which.max(z[i,]))
    for(i in 1:yres) x <- c(x,which.max(z[,i]))
  }
  return(list(x=x,y=y))
}

revisit <- function(d,target="cv",loga.ref=NULL,
  do.plot=NULL,col=par("fg"),yline=2,ylim=NULL,ispecies=NULL,add=FALSE,
  cex=par("cex"),lwd=par("lwd"),mar=NULL,side=1:4,xlim=NULL,labcex=0.6,
  pch=1,legend="",legend.x=NULL,lpch=NULL,main=NULL,lograt.ref=NULL,plot.ext=TRUE) {
  # calculate and plot diversity indices of relative abundances
  # 20090316 jmd
  # d can be the output from diagram (enables plotting)
  # or simply a list of logarithms of activity 
  # (each list entry must have the same dimensions)
  # test if the entries have the same dimensions
  ud <- unique(lapply(1:length(d),function(x) dim(d[[x]])))
  if(length(ud)==1) {
    # d is list of logarithms of activity
    if(missing(do.plot)) do.plot <- FALSE
    if(do.plot) stop("can't make a plot if argument 'd' is not the output from diagram()")
    logact <- d
  } else {
     # d is the output from diagram()
    if(!"logact" %in% names(d)) {
      stop(paste("the list provided in 'd' doesn't look like it came from diagram()",
        "(or, for a 2-D diagram, did you specify mam=FALSE)?"))
    }
    if(missing(do.plot)) do.plot <- TRUE
    logact <- d$logact
  }
  # take a subset (or all) of the species
  if(is.null(ispecies)) ispecies <- 1:length(logact)
  logact <- logact[ispecies]
  # number of species
  ns <- length(logact)
  # the dimensions 
  mydim <- dim(logact[[1]])
  nd <- length(mydim)
  if(nd==1) if(mydim==1) nd <- 0
  cat(paste("revisit: calculating",target,"in",nd,"dimensions\n"))

  ## on to diversity calculations
  # given a list of logarithms of activities of species
  # (as vectors or matrices or higher dimensional arrays) 
  # calculate a diversity index of the same dimensions
  # available targets:
  # "shannon" shannon entropy
  # "sd" standard deviation
  # "cv" coefficient of variation
  # "sd.log", "cv.log" SD/CV for the logarithms of activity
  # "qqr" correlation coefficient on q-q plot (i.e., normality test)
  # "richness" species richness (loga.ref)
  # "cvrmsd" coefficient of variation of rmsd (loga.ref)
  # "spearman" spearman correlation coefficient (loga.ref)
  # "pearson" pearson correlation coefficient (log.ref)
  
  # vectorize the entries in the logact list
  for(i in 1:ns) {
    logact[[i]] <- as.vector(logact[[i]])
    # convert infinite values to NA
    logact[[i]][is.infinite(logact[[i]])] <- NA
  }
  # make a place for the results
  H <- logact[[1]]
  H[] <- 0
  # ratio-specific calculations
  if(!is.null(lograt.ref)) {
    if(!target %in% c("rmsd","cvrmsd","spearman","pearson"))
      stop(paste("target",target,"not available when comparing activity ratios"))
    else(cat("findit: calculating activity ratios\n"))
    # instead of logact we use calculated log activity ratio
    logact <- lograt(loga.ref,logact)
    loga.ref <- lograt.ref
  }
  # target-specific calculations
  target.lower <- tolower(target)
  if(target.lower %in% c("sd","cv")) {
    act <- logact
    for(i in 1:ns) act[[i]] <- 10^logact[[i]]
    # build a matrix; rows are species
    myad <- t(list2array(act))
    # get the standard deviations
    H <- sd(myad)
    # for coefficient of variation, divide by the mean
    if(target.lower=="cv") H <- H / colMeans(myad)
  } else if(target.lower %in% c("sd.log","cv.log")) {
    # build a matrix; rows are species
    myad <- t(list2array(logact))
    # get the standard deviations
    H <- sd(myad)
    # for coefficient of variation, divide by the mean
    if(target.lower=="cv.log") H <- H / colMeans(myad)
  } else if(target=="richness") {
    # given a list of logarithms of activities of species
    # (as matrices or vectors) calculate the richness 
    # make a place for the results
    # loop over species
    for(j in 1:length(logact)) {
      isthere <- logact[[j]] > loga.ref
      H[isthere] <- H[isthere] + 1
    }
  } else if(target=="shannon") {
    # for this calculation we need to loop twice;
    # first to convert logact to act and get acttotal
    act <- logact
    for(i in 1:ns) {
      # exponentiate the logarithmic values
      act[[i]] <- 10^logact[[i]]
      if(i==1) acttotal <- act[[i]] else acttotal <- acttotal + act[[i]]
    }
    # now do the calculation
    for(i in 1:ns) {
      dH <- -act[[i]]/acttotal*log(act[[i]]/acttotal)
      if(!any(is.na(dH))) H <- H + dH
      else(warning(paste("revisit: skipping species",i,"which gives NA")))
    }
  } else if(target=="qqr") {
    # normality test using correlation coefficient of a q-q plot 20100901
    actarr <- list2array(logact)
    qqrfun <- function(i,actarr) {
      y <- actarr[i,]
      # this is to catch errors from qqr (qqnorm)
      out <- try(qqr(y),silent=TRUE)
      if(class(out)=="try-error") out <- NA
      return(out)
    }
    H <- as.numeric(mylapply(1:length(H),qqrfun,actarr))
  } else if(target=="spearman") {
    # spearman rank correlation coefficient 20100912
    actarr <- list2array(logact)
    spearfun <- function(i,actarr) {
      y <- actarr[i,]
      out <- spearman(y,loga.ref)
      return(out)
    }
    H <- as.numeric(mylapply(1:length(H),spearfun,actarr))
  } else if(target=="pearson") {
    # pearson correlation coefficient
    actarr <- list2array(logact)
    pearfun <- function(i,actarr) {
      y <- actarr[i,]
      out <- cor(y,loga.ref)
      return(out)
    }
    H <- as.numeric(mylapply(1:length(H),pearfun,actarr))
  } else if(target=="rmsd") {
    # root mean squared deviation
    actarr <- list2array(logact)
    rmsdfun <- function(i,actarr) {
      y <- actarr[i,]
      out <- rmsd(loga.ref,y)
      return(out)
    }
    H <- as.numeric(mylapply(1:length(H),rmsdfun,actarr))
  } else if(target=="cvrmsd") {
    # coefficient of variation of the root mean squared deviation
    actarr <- list2array(logact)
    cvrmsdfun <- function(i,actarr) {
      y <- actarr[i,]
      out <- cvrmsd(loga.ref,y)
      return(out)
    }
    H <- as.numeric(mylapply(1:length(H),cvrmsdfun,actarr))
  } else stop(paste("specified target",target,"not available"))
  # replace dims
  dim(H) <- mydim
#print(H)
  ## now on to plotting
  # get information about the x-axis
  if(do.plot & nd > 0 & nd < 3) {
    xname <- d$xname
    yname <- d$yname
    xres <- d$xlim[3]
    xrange <- d$xlim[1:2]
    # special operations for pH
    if(xname=="H+") {
      xname <- "pH"
      xrange <- -xrange 
    }
    # the x-values
    xs <- seq(xrange[1],xrange[2],length.out=xres)
  }
  # make plots and return values
  if(nd==0) {
    # a 0-D plot
    if(do.plot) {
      plotted <- FALSE
      if(target=="qqr") {
        actarr <- list2array(logact)
        qqnorm(actarr,col=col,pch=pch,main="")
        qqline(actarr)
        plotted <- TRUE
      } else if(target %in% c("rmsd","cvrmsd","spearman","pearson")) {
        # plot the points
        xlab <- "loga.calc"
        ylab <- "loga.ref"
        if(!is.null(lograt.ref)) {
          xlab <- "lograt.calc"
          ylab <- "lograt.ref"
        }
        plot(actarr,loga.ref,xlab=xlab,ylab=ylab,pch=pch)
        # add a 1:1 line
        lines(range(loga.ref),range(loga.ref),col="grey")
        # add a lowess line
        ls <- loess.smooth(actarr,loga.ref)
        lines(ls$x,ls$y,col="red")
        plotted <- TRUE
      }
      if(plotted) {
        # add a title
        if(missing(main)) main <- paste(target,"=",round(H,3)) 
        title(main=main)
        if(!is.null(legend.x)) {
          if(is.null(lpch)) lpch <- unique(pch)
          legend(legend.x,legend=legend,pch=lpch)
        }
      }
    }
    ret.val <- list(H=H)
  } else if(nd==1) {
    # locate the extremum
    ix <- where.extreme(H,target)
    extval <- H[ix]
    # a 1-D plot
    if(do.plot) {
      if(is.null(ylim)) {
        if(target=="richness") {
            ylim <- 0
            if(max(H) > ylim) ylim <- max(H) + 1
            ylim <- c(0,ylim)
        }
        else ylim <- extendrange(H,f=0.075)
      }
      if(is.null(xlim)) xlim <- xrange
      if(!add) thermo.plot.new(xlim=xlim,ylim=ylim,xlab=axis.label(xname),ylab=target,yline=yline,
        cex=cex,lwd=lwd,mar=mar,side=side)
      # plot the values
      lines(xs,as.vector(H),col=col)
      x <- xs[ix]
      if(plot.ext) abline(v=x,lty=2)
      ret.val <- list(H=H,ix=ix,x=x,extval=extval)
    } else ret.val <- list(H=H,ix=ix,extval=extval)
  } else if(nd==2) {
    # a 2-D plot
    iext <- where.extreme(H,target,do.sat=TRUE)
    # what is the extreme value
    ix <- iext$ix
    iy <- iext$iy
    extval <- H[iy,ix]
    ret.val <- list(H=H,ix=ix,iy=iy,extval=extval) 
    if(do.plot) {
      yres <- d$ylim[3]
      yrange <- d$ylim[1:2]
      #if(yname=="T") yrange <- outvert(yrange,"K")
      if(yname=="H+") {
        yname <- "pH"
        yrange <- -yrange 
      }
      ys <- seq(yrange[1],yrange[2],length.out=yres)
      # start the plot
      if(is.null(xlim)) xlim <- xrange
      if(is.null(ylim)) ylim <- yrange
      if(!add) thermo.plot.new(xlim=xlim,ylim=ylim,xlab=axis.label(xname),ylab=axis.label(yname),
        yline=yline,side=side,cex=cex,mar=mar)
      contour(xs,ys,H,add=TRUE,labcex=labcex)
      # plot the location(s) of the extremum
      points(xs[iext$ix],ys[iext$iy],pch=8,cex=2)
      # show trajectories of the extrema
      iexts <- extremes(H,target)
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
      ret.val <- list(H=H,ix=ix,iy=iy,
        x=xs[ix],y=ys[iy],extval=extval)
    }
  } else {
    # we don't make plots for more than two dimensions
    # just return the values
    ret.val <- list(H=H)
  }
  # return the results
  return(invisible(ret.val))
}
