# CHNOSZ/findit.R
# find the minimum or maximum of a target 
# (eg cv, richness, shannon)
# as a function the specified chemical potentials
# 20100831 jmd

findit <- function(lims=list(),target="cv",n=NULL,iprotein=NULL,do.plot=TRUE,
  T=25,P="Psat",res=NULL,labcex=0.6,loga.ref=NULL,as.residue=FALSE,
  loga.tot=0) {
  # the lims list has the limits of the arguments to affinity()
  # we iteratively move toward a higher/lower value of the target
  # within these limits
  
  # fun stuff: when running on either side of 100 deg C,
  # set P=1 to force an extreme of some functions near that temperature
  
  # the number of dimensions
  nd <- length(lims)

  # get protein lengths if iprotein is given
  if(!is.null(iprotein)) pl <- protein.length(-iprotein)

  # how many gradations in each dimension we can timely 
  # consider goes down with the number of dimensions
  if(is.null(res)) res <- c(128,64,16,8,6,4,4)[nd]
  if(is.null(n)) n <- c(4,6,6,8,12,12,12)[nd]
  # the size fraction of the new intervals after each iteration
  rat <- 0.85
  # we can zoom in more quickly if the gradations are smaller
  if(res > 4) rat <- 0.8
  if(res > 8) rat <- 0.75
  if(res > 16) rat <- 2/(1+sqrt(5))

  # the initial values of the guesses (if midpoint==FALSE)
  basis <- thermo$basis

  # a hack so that we can use pH as a variable
  if("pH" %in% names(lims)) {
    iH <- match("H+",rownames(basis))
    if(length(iH) > 0) {
      rownames(basis)[iH] <- "pH"
      basis$logact[iH] <- -basis$logact[iH]
    } else(stop("pH is a requested variable but H+ is not in the basis"))
  }

  # function to get the current limits of a variable
  limfun <- function(lim,curr,i) {
    # in the ith loop we consider intervals
    # of rat^(i-1) of the ranges specified at start
    cat(paste("now",round(curr,4),"",""))
    range <- range(lim)
    if(i==1) int <- abs(diff(range)) 
    else int <- abs(diff(range)) * rat^(i-1)
    mylims <- c(curr-int/2,curr+int/2)
    # don't go beyond the starting limits
    if(any(mylims < min(range))) mylims <- mylims + (min(range) - min(mylims))
    if(any(mylims > max(range))) mylims <- mylims - (max(mylims) - max(range))
    # reverse the order if the axis is reversed
    if(diff(lim) < 0) mylims <- rev(mylims)
    cat(paste("limits",round(mylims[1],4),round(mylims[2],4),"\n"))
    return(mylims)
  }

  # place to put the output values
  teststat <- numeric()
  out <- vector("list",length(lims))
  names(out) <- names(lims)
  lolim <- out
  hilim <- out
  outlims <- lims

  # loop for the specified number of iterations
  # (todo: loop until an error threshhold is reached)
  for(i in 1:n) {
    cat(paste("\n###### findit: iteration",i,"of",n,"\n"))
    # to generate the argument list for affinity
    # with the variables specified in lims
    aargs <- list()
    for(j in 1:length(lims)) {
      if(names(lims)[j] %in% rownames(basis)) {
        ibasis <- match(names(lims)[j],rownames(basis))
        cat(paste("###",rownames(basis)[ibasis],""))
        # the starting search limits for this species
        lim <- lims[[j]]
        # center the search interval on the current values
        curr <- basis$logact[ibasis]
        # update the argument list with the intervals
        myarg <- list(c(limfun(lim,curr,i),res))
        names(myarg) <- rownames(basis)[ibasis]
        aargs <- c(aargs,myarg)
      } else if(names(lims[j]) %in% c("T","P")) {
        cat(paste("###",names(lims[j]),""))
        if(names(lims[j])=="T") {
          lim <- lims$T
          curr <- T
        } else {
          lim <- lims$P
          curr <- P
        }
        myarg <- list(c(limfun(lim,curr,i),res))
        names(myarg) <- names(lims[j])
        aargs <- c(aargs,myarg)
      } else warning(paste("findit: ignoring",names(lims[j]),"which is not a basis species, T or P"))
    }
    # we include single values of T and P in the arguments
    # if their limits weren't given
    if(!"T" %in% names(lims)) aargs <- c(aargs,list(T=T))
    if(!"P" %in% names(lims)) aargs <- c(aargs,list(P=P))
    # include iprotein if given
    if(!is.null(iprotein)) {
      aargs <- c(aargs,list(iprotein=iprotein))
    }

    # now calculate the affinities
    a <- do.call(affinity,aargs)
    # then calculate the values of the target function
    d <- diagram(a,do.plot=FALSE,mam=FALSE,as.residue=as.residue,logact=loga.tot)
    dd <- revisit(d$logact,target,loga.ref=loga.ref)$H
    # find the extreme value
    iext <- where.extreme(dd,target)
    # find the extreme value
    teststat <- c(teststat,dd[iext])
    # what are its coordinates
    dd1 <- dd
    dd1[] <- 1
    ai <- which(dd1==1,arr.ind=TRUE)
    if(nd==1) ai <- ai[iext]
    else ai <- ai[iext,]

    # loop to update the current parameters
    for(j in 1:length(lims)) {
      # the current limits for this variable
      # we ignore T and P if they are not variable:
      mylims <- aargs[[j]]
      # the increments used
      myinc <- seq(mylims[1],mylims[2],length.out=mylims[3])
      # the value of the variable at the extreme of the function
      myval <- myinc[ai[j]]
      # update the basis table, T or P
      if(names(lims)[j] %in% rownames(basis)) {
        ibasis <- match(names(lims)[j],rownames(basis))
        basis$logact[ibasis] <- myval
        basis(rownames(basis)[ibasis],myval)
      } else if(names(lims)[j]=="T") {
        T <- myval
      } else if(names(lims)[j]=="P") {
        P <- myval
      }
      # update the output list
      out[[j]] <- c(out[[j]],myval)
      lolim[[j]] <- c(lolim[[j]],mylims[1])
      hilim[[j]] <- c(hilim[[j]],mylims[2])
      outlims[[j]] <- mylims
    }

    # are we making a plot?
    if(do.plot) {
      # add our search lines and extreme points
      # to the plot
      if(nd==1) {
        if(i==1) revisit(d,target,loga.ref,xlim=lims[[1]])
        # on a 1-D diagram we add vertical lines to show our search limits
        abline(v=outlims[[1]][1:2])
        lines(myinc,dd)
        points(myval,dd[iext])
      } else if(nd==2) {
        if(i==1) revisit(d,target,loga.ref,xlim=lims[[1]],ylim=lims[[2]],labcex=labcex)
        else {
          # on a 2-D diagram we add a box around our search limits
          # and an updated map for this region
          ol1 <- outlims[[1]]
          ol2 <- outlims[[2]]
          rect(ol1[1],ol2[1],ol1[2],ol2[2],border=par("fg"),col="white")
          revisit(d,target,loga.ref,xlim=lims[[1]],ylim=lims[[2]],add=TRUE,labcex=labcex)
        }
        text(out[[1]],out[[2]])
        points(out[[1]],out[[2]],cex=2)
      } else {
        if(i==1) add <- FALSE else add <- TRUE
        # we only draw a 2-D diagram even if we're optimizing
        # for more than two variables
        # paint an empty box over the search area
        if(i > 1) {
          ol1 <- outlims[[1]]
          ol2 <- outlims[[2]]
          rect(ol1[1],ol2[1],ol1[2],ol2[2],border=par("fg"),col="white")
        }
        # we have to extract slices for dimensions greater than two
        # we extract the third dimension until only two remain
        for(j in 3:nd) {
          # ai[j] - which slice in this dimension has the extremum
          for(k in 1:length(d$logact)) d$logact[[k]] <- slice(d$logact[[k]],3,ai[j])
        }
        # now make the plot
        revisit(d,target,loga.ref,xlim=lims[[1]],ylim=lims[[2]],add=add,labcex=labcex)
        # indicate the location of the extremum
        text(out[[1]],out[[2]])
        points(out[[1]],out[[2]],cex=2)
      }
      # are we making a legend?
      # it looks something like
      #  n     qqr     O2     H2O   CO2   NH3
      #  1   0.993 -85.21  -6.378 2.032 -6.93
      # ...
      # 12   0.997 -85.26  -6.488 3.111 -6.50
    }
  }  # end main loop
  # build return list: values of chemical variables and test statistic
  teststat <- list(teststat)
  names(teststat) <- target
  value <- c(out,teststat)
  # build return list: limits at each step
  out <- list(value=value,lolim=lolim,hilim=hilim)
  class(out) <- "findit"
  return(out)
}

# to plot results from findit()
plot.findit <- function(x,which=NULL,mar=c(3.5,5,2,2),xlab="iteration",...) {
  # show the values of the test statistic and of the
  # chemical variables at each iteration
  l <- length(x$value)
  if(is.null(which)) which <- 1:l
  l <- length(which)
  opar <- par(mfrow=c(l,1),mar=mar) 
  for(i in which) {
    n <- length(x$value[[i]])
    ylab <- names(x$value)[i]
    if(ylab %in% c(rownames(thermo$basis),"T","P","pH","Eh")) ylab <- axis.label(ylab)
    # the values
    plot(1:n,x$value[[i]],xlab=xlab,ylab=ylab,...)
    lines(1:n,x$value[[i]])
    # the search intervals
    # (not available for the test statistic)
    if(i!=length(x$value)) {
      lines(1:n,x$lolim[[i]],lty=2)
      lines(1:n,x$hilim[[i]],lty=2)
    }
  }
  par(opar)
}

