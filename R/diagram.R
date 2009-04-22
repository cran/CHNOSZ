# CHNOSZ/diagram.R
# Copyright (C) 2006-2009 Jeffrey M. Dick
# plot predominance or activity diagrams 
# from affinities of formation reactions
# 20061023 jmd


diagram <- function(affinity,ispecies=NULL,balance=NULL,
  names=NA,color=NA,add=FALSE,dotted=0,cex=par('cex'),col=par('col'),pe=TRUE,pH=TRUE,
  ylim=c(-4,0),ylog=TRUE,title=NULL,cex.names=1,legend.x='topright',lty=NULL,
  col.names=par('fg'),cex.axis=par('cex'),logact=NA,property=NULL,lwd=par('lwd'),
  alpha=FALSE,mar=NULL,residue=FALSE,yline=par('mgp')[1]+1,xrange=NULL,
  ylab=NULL,xlab=NULL,do.plot=TRUE,as.residue=FALSE,mam=TRUE) {

  # store input values
  aval <- affinity$values
  asl <- affinity$species$logact
  # number of species possible
  nspecies <- length(aval)
  # number of dimensions (T, P or chemical potentials that are varied)
  nd <- 0
  if(length(affinity$xlim) > 1) nd <- nd + 1
  if(length(affinity$ylim) > 1) nd <- nd + 1
  # 'ispecies' which species to consider
  if(is.null(ispecies)) {
    ispecies <- 1:nspecies
    # take out species that have NA affinities
    if(nd > 0) {
      for(i in 1:nspecies) if(all(is.na(aval[[i]]))) ispecies[i] <- NA
      ispecies <- ispecies[!is.na(ispecies)]
    }
  }
  if(!identical(ispecies,1:nspecies)) {
    cat(paste('diagram: using',length(ispecies),'of',nspecies,'species.\n'))
    affinity$species <- affinity$species[ispecies,]
    aval <- aval[ispecies]
  }
  # number of species that are left
  nspecies <- length(aval)
  # handle line type/width arguments
  if(do.plot) {
    if(is.null(lty)) lty <- 1:nspecies
    if(is.null(lwd)) lwd <- rep(1,nspecies)
    if(length(lty)!=nspecies) lty <- rep(lty,length.out=nspecies)
    if(length(lwd)!=nspecies) lwd <- rep(lwd,length.out=nspecies)
  }
  # covert from activities of proton or electron to pH or pe
  if(affinity$xname=='H+' & pH) { affinity$xname <- 'pH'; affinity$xlim[1:2] <- -affinity$xlim[1:2] }
  if(affinity$xname=='e-' & pe) { affinity$xname <- 'pe'; affinity$xlim[1:2] <- -affinity$xlim[1:2] }
  if(affinity$yname=='H+' & pH) { affinity$yname <- 'pH'; affinity$ylim[1:2] <- -affinity$ylim[1:2] }
  if(affinity$yname=='e-' & pe) { affinity$yname <- 'pe'; affinity$ylim[1:2] <- -affinity$ylim[1:2] }
  # T/P conversions
  if(nuts('T')=='C') {
    if(affinity$xname=='T') { affinity$xlim[1:2] <- convert(affinity$xlim[1:2],nuts('T')) }
    if(affinity$yname=='T') { affinity$ylim[1:2] <- convert(affinity$ylim[1:2],nuts('T')) }
  }
  if(nuts('P')=='MPa') {
    if(affinity$xname=='P') { affinity$xlim[1:2] <- convert(affinity$xlim[1:2],nuts('P')) }
    if(affinity$yname=='P') { affinity$ylim[1:2] <- convert(affinity$ylim[1:2],nuts('P')) }
  }
  # the property we're plotting
  if(is.null(property)) {
    property <- affinity$property
    if(property=='A') property <- 'affinity'
  }
  # 'balance' how to balance reactions
  # (name of basis species, 'PBB', 1, or NULL for auto-select)
  # (or vector to set nbalance)
  # 'nbalance' vector of coefficients of the balanced component 
  isprotein <- grep('_',as.character(affinity$species$name))
  if(missing(balance) & length(isprotein)==nspecies)
    balance <- 'PBB'
  if(is.null(balance)) {
    ib <- which.balance(affinity$species)
    if(!is.null(ib)) {
      balance <- rownames(basis())[ib[1]]
      nbalance <- affinity$species[,ib[1]]
      cat(paste('diagram: balance method - per mole of ',balance,'.\n',sep=''))
    } else {
      balance <- 1
      nbalance <- rep(1,length(ispecies))
      cat(paste('diagram: balance method - per mole of species.\n'))
    }
  } else {
   if(is.character(balance[1])) {
      # is the balance is "PBB" for protein backbone
      if(identical(balance,'PBB')) {
        if(length(isprotein) != nspecies)
          stop('diagram: PBB was the requested balance, but ',
            affinity$species$name[-isprotein],' is not protein.')
        nbalance <- as.numeric(protein.length(as.numeric(names(aval))))
        cat(paste('diagram: balance method - per mole of protein backbone.\n'))
      } else {
      # is the balance the name of a basis species
        ib <- match(balance,rownames(basis()))
        if(!is.na(ib)) {
          nbalance <- affinity$species[,ib]
          if(TRUE %in% (nbalance==0)) {
            inotbalance <- which(nbalance==0)
            # take out species that can't be balanced
            if(length(inotbalance)==1) hh <- "has" else hh <- "have"
            cat(paste('diagram: removing ',c2s(affinity$species$name[inotbalance]),
              ', which ',hh,' no ',balance,'\n',sep=''))
            aval <- aval[-inotbalance]
            affinity$species <- affinity$species[-inotbalance,]
            nbalance <- nbalance[-inotbalance]
            nspecies <- nspecies - length(inotbalance)
          }
        } else {
          stop('requested balance (',balance,') is not a the name of a basis species or other value I understand.')
        }
      }
    } else if(is.numeric(balance[1])) {
      # user-defined balance
      nbalance <- rep(balance,length.out=nspecies)
      if(all(nbalance==1)) balance <- "species" else balance <- "user-defined"
    }
  }
  if(length(nbalance) < 100) cat(paste('diagram: balance coefficients are ',c2s(nbalance),'\n',sep=''))
  # 2D diagram mam and residue both FALSE can lead to long calculation times
  if(nd==2 & !mam & !residue) {
    warning("expect long calculation for 2D diagram with mam and residue both FALSE",immediate.=TRUE)
  }
  # rewrite the balance vector for residue reactions
  if(residue) {
    if(any(nbalance < 0)) stop("negative balance prohibits using residue reactions")
    # affinities divided by balance
    for(i in 1:nspecies) aval[[i]] <- aval[[i]]/nbalance[i]
    oldlogact <- affinity$species$logact
    affinity$species$logact <- affinity$species$logact + log10(nbalance)
    oldbalance <- nbalance
    nbalance <- rep(1,length(nbalance))
    cat(paste('diagram: using residue reactions\n'))
  }
  # get labels for the plot
  if(!is.null(names)) {
    if(is.na(names[1])) {
      # if names are missing ... make them up
      # properties of basis species or reactions?
      if(affinity$property %in% c('G.basis','logact.basis')) names <- rownames(affinity$basis)
      else {
        names <- as.character(affinity$species$name)
        # remove common organism label for proteins
        if(all(grep("_",names)>0)) {
          pname <- oname <- character()
          for(i in 1:length(names)) {
            mynames <- s2c(names[i],sep="_",keep.sep=FALSE)
            pname <- c(pname,mynames[1])
            oname <- c(oname,mynames[2])
          }
          if(length(unique(oname))==1) for(j in 1:length(names)) names[j] <- pname[j]
        }
        # append species label to distinguish ambiguous names
        #if(length(unique(names))!=length(names)) 
        #  names <- paste(names,' (',affinity$species$state,')',sep='')
        isdup <- names %in% names[duplicated(names)]
        if(any(isdup)) names[isdup] <- paste(names[isdup],
          " (",affinity$species$state[isdup],")",sep="")
      }
    }
  }
  # figure out colors
  if(missing(color)) {
    if(nd==2) {
      if( add | names(dev.cur())=='postscript' ) color <- NULL
      else color <- 'heat'
    } else color <- 'black'
  }
  if(!is.null(color)) {
    color <- rep(color,length.out=nspecies)
    if(is.character(color[1])) {
      if(color[1]=='rainbow') color <- rainbow(nspecies)
      else if(color[1]=='heat') color <- heat.colors(nspecies)
    }
  }
  ### end variable assignment

  ### various plotting functions
  ### color plot function
  plot.color <- function(xs,ys,out,color,nspecies) {
    # handle min/max reversal
    if(xs[1] > xs[length(xs)]) {
      tc <- out; t <- numeric()
      for(i in 1:length(xs)) {
        t[i] <- xs[length(xs)+1-i]
        tc[,i] <- out[,length(xs)+1-i]
      }
      out <- tc; xs <- t
    }
    if(ys[1] > ys[length(ys)]) {
      tc <- out; t <- numeric()
      for(i in 1:length(ys)) {
        t[i] <- ys[length(ys)+1-i]
        tc[i,] <- out[length(ys)+1-i,]
      }
      out <- tc; ys <- t
    }
    # the z values
    zs <- out
    for(i in 1:nrow(zs)) zs[i,] <- out[nrow(zs)+1-i,]
    zs <- t(zs)
    breaks <- c(0,1:nspecies) + 0.5
    image(x=xs,y=ys,z=zs,col=color,add=TRUE,breaks=breaks)
  }
  ### curve plot function
  plot.curve <- function(out,xlim,ylim,dotted,col,xrange) {
    # dotted tells us what fraction of points to ignore,
    # to make the line look dotted
    # what is the size of each grid point?
    # this is to be consistent with the way image() plots the grid,
    # so divide by (ncol-1) and (ncol-2), as half of each of the grid
    # boxes on the edges lie outside our limits
    xsize <- (xlim[length(xlim)]-xlim[1])/(ncol(out)-1)
    ysize <- (ylim[length(ylim)]-ylim[1])/(nrow(out)-1)
    # now, the coordinates of the upper left vertex of the upper left grid box
    # if xsize/2 and ysize/2 is used here, the results 
    # aren't aligned with the output of image()
    xleft <- xlim[1] - xsize
    yupper <- ylim[length(ylim)] + ysize
    # functions to plot boundaries around a given grid point
    vline <- function(ix1,ix2,iy) {
      # only plot dotted grid points if requested
      # (i.e. those points for which there is no remainder
      # when their coordinate is divided by 'dotted')
      xy.na <- list(x=NA,y=NA)
      if(!identical(dotted,0)) if(0 %in% (iy%%dotted)) return(xy.na)
      x1 <- xleft + ix1 * xsize
      x2 <- xleft + ix2 * xsize
      y  <- yupper - iy * ysize
      x <- (x1+x2)/2; y1 <- y - ysize/2; y2 <- y + ysize/2
      # xrange: don't plot predominance field boundaries outside this range
      if(!is.null(xrange)) {
        if(x < xrange[1] | x > xrange[2]) return(xy.na)
      }
      return(list(x=c(x,x),y=c(y1,y2)))
    }
    hline <- function(ix,iy1,iy2) {
      xy.na <- list(x=NA,y=NA)
      if(!identical(dotted,0)) if(0 %in% (ix%%dotted)) return(xy.na)
      x <- xleft + ix * xsize
      y1  <- yupper - iy1 * ysize
      y2  <- yupper - iy2 * ysize
      y <- (y1+y2)/2; x1 <- x - xsize/2; x2 <- x + xsize/2
      if(!is.null(xrange)) {
        if(x1 < xrange[1] | x1 > xrange[2]) return(xy.na)
        if(x2 < xrange[1] | x2 > xrange[2]) return(xy.na)
      }
      return(list(x=c(x1,x2),y=c(y,y)))
    }
    # plot curves to identify the boundaries
    # 20080910 gather them in myx and myy to call lines() only once
    myx <- myy <- NA
    for(i in 1:nrow(out)) {
      for(j in 1:ncol(out)) {
        xy <- xy.na <- list(x=NA,y=NA)
        xyfun <- function(xy1,xy2,xy3) list(x=c(xy1$x,xy2$x,xy3$x),y=c(xy1$y,xy2$y,xy3$y))
        # left, up (lu)
        if(!j==1) if(out[i,j]!=out[i,j-1]) xy <- xyfun(xy,vline(j,j-1,i),xy.na)
        if(!i==1) if(out[i,j]!=out[i-1,j]) xy <- xyfun(xy,hline(j,i,i-1),xy.na)
        # right, down
        if(!j==ncol(out)) if(out[i,j]!=out[i,j+1]) xy <- xyfun(xy,vline(j,j+1,i),xy.na)
        if(!i==nrow(out)) if(out[i,j]!=out[i+1,j]) xy <- xyfun(xy,hline(j,i,i+1),xy.na)
        lxy <- length(myx)
        # add values if we found a boundary
        if( !all(is.na(c(xy$x,xy$y))) ) {
          myx <- c(myx,xy$x)
          myy <- c(myy,xy$y)
        }
      }
    }
    lines(myx,myy,col=col,lwd=lwd[1])
  }
  ### label plot function
  # calculate coordinates for field labels
  plot.names <- function(out,xs,ys,names) {
    ll <- nspecies
    lx <- numeric(ll); ly <- numeric(ll); n <- numeric(ll)
    #for(i in 1:nspecies) {
    #ix <- 0; iy <- 0; n <- 0
    for(j in nrow(out):1) {
      for(k in 1:ncol(out)) {
        i <- out[j,k]
        lx[i] <- lx[i] + xs[k]
        ly[i] <- ly[i] + ys[nrow(out)+1-j]
        n[i] <- n[i] + 1
      }
    }
    lx <- lx[n!=0]
    ly <- ly[n!=0]
    is <- n!=0
    n <- n[n!=0]
    lx <- lx/n
    ly <- ly/n
    # plot field labels
    # the cex argument in this function specifies the character 
    # expansion of the labels relative to the current
    text(lx,ly,labels=names[is],cex=cex.names,col=col.names)
  }
  ### property description
  axis.title <- function(property) {
    if(property=='A') return('A/2.303RT')
    else if(property=='logact.basis') return('logQ*')
    else if(!property %in% c("logK","logQ")) return(paste("Delta",property))
    else return(property)
  }
  balance.title <- function(balance) {
    if(balance==1) return('per mole of product')
    else if(is.numeric(balance)) return(paste('per',balance,'moles of product'))
    else return(paste('per mole of',balance))
  }
  ### done with function definitions

  ### now on to the calculation and plots
  # do speciation calculations
  # unless using mam (maximum affinity method) for 2-D diagram
  if(property=='affinity' & (!mam | nd==0 | nd==1) ) {
    # compute the abundances of species
    # total logarithm of activity of the balanced component
    ib <- match(balance,colnames(affinity$species))
    if(!is.numeric(logact)) {
      thisa <- sum(10^affinity$species$logact * nbalance)
      if(thisa < 0) thisa <- -thisa
      logatotal <- log10(thisa)
      if(can.be.numeric(balance)) cat('diagram: log total activity of species is ',logatotal,'\n',sep='')
      else cat('diagram: log total activity of ',balance,' (from species) is ',logatotal,'\n',sep='')
    } else {
      logatotal <- logact
      cat('diagram: log total activity of ',balance,' (from argument) is ',logatotal,'\n',sep='')
    }
    # Astar: the affinities/2.303RT of formation of species
    # in their standard-state activities
    Astar <- aval
    if(residue) {
      for(j in 1:length(Astar)) Astar[[j]] <- Astar[[j]] + oldlogact[j]/oldbalance[j]
      A <- abundance(Astar,nbalance,logatotal)
    }
    else {
      for(j in 1:length(Astar)) Astar[[j]] <- Astar[[j]] + affinity$species$logact[j]
      A <- abundance.old(Astar,aval,nbalance,logatotal)
    }
    # if we rewrote the formation reactions per residue,
    # get back to activities of species
    if(residue & !as.residue) for(j in 1:length(A)) A[[j]] <- A[[j]] - log10(oldbalance[j])
  }
  ### 0-D properties of species or reactions for single set of conditions
  if(nd==0) {
    if(do.plot) {
      mgp <- par("mgp")
      mgp[1] <- yline
      if(property=="affinity") {
        # the logarithms of activities
        v <- numeric()
        for(i in 1:length(A)) v <- c(v,A[[i]])
        barplot(v,names.arg=names,ylab=axis.title(affinity$property),mgp=mgp)
        if(missing(title)) title <- "logarithm of activity"
        title(main=title)
      } else {
        # the value of a property like affinity, standard gibbs energy etc
        for(i in 1:nspecies) aval[[i]] <- aval[[i]]/nbalance[i]
        v <- as.numeric(aval)
        barplot(v,names.arg=names,ylab=axis.title(affinity$property),mgp=mgp)
        if(!is.null(title)) title(main=title)
        else title(main=paste(axis.title(affinity$property),balance.title(balance)))
      }
    }
    if(property=="affinity") return(invisible(list(logact=A)))
    else return(invisible(aval))
  }
  ### 1-D (property or speciation) diagram
  if(nd==1) {
    if(missing(color)) color <- rep(par('col'),nspecies)
    xvalues <- seq(affinity$xlim[1],affinity$xlim[2],length.out=affinity$xlim[3])
    xlab <- axis.label(affinity$xname,as.character(affinity$basis$state[match(affinity$xname,rownames(affinity$basis))]))
    if(property == 'affinity') {
      # alpha: plot degree of formation instead of logact
      # scale the activities to sum=1
      if(alpha) {
        ylog <- FALSE
        if(missing(ylim)) ylim <- c(0,1)
        for(i in 1:length(A[[1]])) {
          a <- numeric()
          for(j in 1:length(A)) a <- c(a,A[[j]][i])
          loga.sum <- log10(sum(10^a))
          loga <- 0; a <- loga - loga.sum
          for(j in 1:length(A)) A[[j]][i] <- A[[j]][i] + a
        }
      }
      if(do.plot) {
        if(ylog) {
          if(is.null(ylab)) ylab <- as.expression(quote(log~italic(a)))
          if(!add) thermo.plot.new(xlim=affinity$xlim[1:2],
            ylim=ylim,xlab=xlab,ylab=ylab,cex=cex,cex.axis=cex.axis,mar=mar,yline=yline)
          for(i in 1:length(A)) {
            lines(xvalues,(A[[i]]),col=color[i],lty=lty[i],lwd=lwd[i])
          }
        } else {
          # are we plotting alphas or logactivities?
          if(is.null(ylab)) {
            if(alpha) ylab <- as.expression(quote(alpha))
            else ylab <- as.expression(quote(a[i]))
          }
          if(missing(ylim)) ylim <- c(0,1)
          if(!add) thermo.plot.new(xlim=affinity$xlim[1:2],
            ylim=ylim,xlab=xlab,ylab=ylab,cex=cex,cex.axis=cex.axis,mar=mar,yline=yline)
          for(i in 1:length(A)) lines(xvalues,10^(A[[i]]),col=color[i],lty=lty[i],lwd=lwd[i])
        }
      }
    } else {
      A <- aval
      v <- numeric()
      # plot logarithms of activities or values of properties
      if(property=="affinity") for(i in 1:length(A)) v <- c(v,A[[i]][1,])
      else for(i in 1:nspecies) {
        aval[[i]] <- aval[[i]]/nbalance[i]
        v <- c(v,aval[[i]])
      }
      # determine the overall maximum and minimum values
      if(is.null(ylim)) ylim <- c(min(v),max(v))
      if(!add) thermo.plot.new(xlim=affinity$xlim[1:2],ylim=ylim,xlab=xlab,ylab=axis.title(affinity$property),cex=cex,mar=mar,yline=yline)
      for(i in 1:length(A)) 
        lines(xvalues,as.numeric(A[[i]]),col=color[i],lty=lty[i])
      if(!is.null(title)) title(main=title)
      else title(main=paste(axis.title(affinity$property),balance.title(balance)))
    }
    if(affinity$property=='logact.basis') names <- rownames(affinity$basis)
    if(do.plot & !add & !is.null(legend.x)) legend(x=legend.x,lty=lty,legend=names,col=color,bg=par('bg'),cex=cex.names,lwd=lwd)
    if(alpha) for(i in 1:length(A)) A[[i]] <- 10^A[[i]]
    # 20090324 return list with single element 'logact'
    if(property=="affinity") return(invisible(list(basis=affinity$basis,species=affinity$species,
      T=affinity$T,P=affinity$P,xname=affinity$xname,xlim=affinity$xlim,yname=affinity$yname,
      ylim=affinity$ylim,logact=A)))
    else return(invisible(aval))
  }
  ### 2-D predominance diagram aka equal activity diagram
  if(property=='affinity') {
    # we don't do only one species
    if(length(aval)==1) stop('refusing to make a predominance diagram for a single species')
    # predict predominant species
    if(mam) {
      # mam=TRUE: maximum affinity method
      if(residue) {
        # account for length 1: affinity of residues
        for(i in 1:length(aval))
          aval[[i]] <- aval[[i]] + oldlogact[i]/oldbalance[i]
        # account for length 2: affinity of species
        if(!as.residue)
          for(k in 1:length(aval)) 
            aval[[k]] <- aval[[k]] - log10(oldbalance[k])
      } else {
        # we're not interested in the residues but we still
        # normalize the affinities by the balanced component
        for(k in 1:length(aval)) aval[[k]] <- aval[[k]]/nbalance[[k]]
      }
      myvalues <- aval
    } else {
      # mam=FALSE: use abundances from speciation calculation
      myvalues <- A
    }
    # out: the index of the predominant species
    out <- myvalues[[1]]
    for(j in 1:nrow(out)) {
      values <- list()
      for(k in 1:nspecies) values[[k]] <- myvalues[[k]][j,]
      out[j,] <- which.pmax(values,na.rm=TRUE)
    }
    # reorder the rows in out
    o <- out
    for(i in 1:nrow(out)) o[i,] <- out[nrow(out)+1-i,]
    out <- o
    # the x, y and z values 
    xs <- seq(affinity$xlim[1],affinity$xlim[2],length.out=affinity$xlim[3])
    ys <- seq(affinity$ylim[1],affinity$ylim[2],length.out=affinity$ylim[3])
    # 20090217 added test for do.plot here
    if(do.plot) {
      if(!add) {
        xstate=as.character(affinity$basis$state[match(affinity$xname,rownames(affinity$basis))])
        ystate=as.character(affinity$basis$state[match(affinity$yname,rownames(affinity$basis))])
        if(is.null(xlab)) xlab <- axis.label(as.character(affinity$xname),xstate)
        if(is.null(ylab)) ylab <- axis.label(as.character(affinity$yname),ystate)
        thermo.plot.new(xlim=affinity$xlim[1:2],ylim=affinity$ylim[1:2],
          xlab=xlab,ylab=ylab,cex=cex,cex.axis=cex.axis,mar=mar,yline=yline)
      }
      # colors and curves
      if(!is.null(color)) plot.color(xs,ys,out,color,nspecies)
      if(!is.null(names)) plot.names(out,xs,ys,names)
      if(!is.null(dotted)) plot.curve(out,affinity$xlim[1:2],
        affinity$ylim[1:2],dotted,col,xrange=xrange)
      # return the formation reactions as they were balanced
      species <- affinity$species
      basis <- affinity$basis
      for(i in 1:nrow(basis)) species[,i] <- species[,i]/nbalance
    }
    # give back the results
    if(mam) return(invisible(list(basis=affinity$basis,species=species,T=affinity$T,P=affinity$P,
      xname=affinity$xname,xlim=affinity$xlim,yname=affinity$yname,ylim=affinity$ylim,out=out,aval=aval)))
    else return(invisible(list(basis=affinity$basis,species=species,T=affinity$T,P=affinity$P,xname=affinity$xname,
      xlim=affinity$xlim,yname=affinity$yname,ylim=affinity$ylim,out=out,aval=aval,logact=A)))
  } else stop(paste('diagram: 2-D plot of',property,'not available'))
}

abundance <- function(Astar,nbalance,thisloga) {
  # 20090217 new "abundance" function
  # return logactivity of species
  # works using Boltzmann distribution
  # A/At = e^(Astar/nbalance) / sum(e^(Astar/nbalance))
  # A/At = e^(Astar/nbalance) / sum(e^(Astar/nbalance))
  # where A is activity of the ith residue and
  # At is total activity of residues
  # advantages over abundance.old
  # 1) works on vectors (also matrices - 2D plots now possible)
  # 2) loops over species only - way faster
  # 3) always works (no root finding games)
  # disadvantage:
  # 1) only works for residue reactions
  # 2) can create NaN logacts if the Astars are huge/small

  # initialize output object
  A <- Astar
  # remember the dimensions of elements of Astar (could be vector,matrix)
  Astardim <- dim(Astar[[1]])
  # first loop: make vectors
  A <- mylapply(1:length(A),function(i) as.vector(A[[i]]))
  # second loop: get the exponentiated Astars (numerators)
  # need to convert /2.303RT to /RT
  #A[[i]] <- exp(log(10)*Astar[[i]]/nbalance[i])/nbalance[i]
  A <- mylapply(1:length(A),function(i) exp(log(10)*Astar[[i]]/nbalance[i]))
  # third loop: accumulate the denominator
  # initialize variable to hold the sum
  At <- A[[1]]; At[] <- 0
  for(i in 1:length(A)) At <- At + A[[i]]*nbalance[i]
  # fourth loop: calculate log abundances and replace the dimensions
  A <- mylapply(1:length(A),function(i) thisloga + log10(A[[i]]/At))
  # fifth loop: replace dimensions
  for(i in 1:length(A)) dim(A[[i]]) <- Astardim
  # we're done!
  return(A)
}

abundance.old <- function(Astar,av,nbalance,thisloga) {
  # 20090217 extracted from diagram and renamed to abundance.old
  # to turn the affinities/RT (A) of formation reactions into 
  # logactivities of species (logact(things)) at metastable equilibrium
  # 20080217 idea: for any reaction stuff = thing,
  # logQ = logact(thing) - logact(stuff),
  # A = logK - logQ = logK - logact(thing) + logact(stuff),
  # logact(thing) = Astar - A
  # where Astar = logK + logact(stuff)
  # ( Astar = A + logact(thing) )
  # and Abar = ( Astar - logact(thing) ) / n
  # ( or logact(thing) = Astar + Abar * n )
  # where thing has n of the balanced quantity
  # below, j indexes species and i indexes conditions
  # remember the dimensions (could be vector,matrix)
  Adim <- dim(Astar[[1]])
  # first loop: make vectors
  for(i in 1:length(Astar)) {
    Astar[[i]] <- as.vector(Astar[[i]])
    av[[i]] <- as.vector(av[[i]])
  }
  A <- Astar2 <- av

  # some function definitions
  # calculate logact(thing) from Abar and Astar
  activityfun <- function(Abar,j,i) (Astar[[j]][i] - Abar * nbalance[j])
  # difference between total activity of balanced quantity
  # computed from affinities and the mass-balanced quantity (thisloga)
  activityfun2 <- function(Abar,i) {
    act <- 0
    for(j in 1:length(A)) act <- act + (10^activityfun(Abar,j,i))*nbalance[j]
    if(act < 0) {
      act <- -act
      logact <- log10(act)
      diff <- thisloga - logact
    } else {
      logact <- log10(act)
      diff <- logact - thisloga
    }
    return(diff)
  }
  for(i in 1:length(A[[1]])) {
    # gather the min/max values of original affinities
    # to constrain our search interval
    Abar.max <- Abar.min <- NULL
    for(j in 1:length(A)) {
      thisAbar <- (A[[j]][i])/nbalance[j]
      thisAbarstar <- (Astar[[j]][i])/nbalance[j]
      thisAbarstar2 <- (Astar2[[j]][i])/nbalance[j]
      if(!is.infinite(activityfun2(thisAbar,i))) {
        if(is.null(Abar.max)) Abar.max <- thisAbar
          else if(thisAbar > Abar.max) Abar.max <-thisAbar
        if(is.null(Abar.min)) Abar.min <- thisAbar
          else if(thisAbar < Abar.min) Abar.min <-thisAbar
      }
      if(!is.infinite(activityfun2(thisAbarstar,i))) {
        if(is.null(Abar.max)) Abar.max <- thisAbarstar
          else if(thisAbarstar > Abar.max) Abar.max <-thisAbarstar
        if(is.null(Abar.min)) Abar.min <- thisAbarstar
          else if(thisAbarstar < Abar.min) Abar.min <-thisAbarstar
      }
    }
    # make sure A.min < A.max
    #if(Abar.min >= Abar.max) Abar.min <- Abar.max - 1
    # avoid the complication when using uniroot,
    # "f() values at end points not of opposite sign"
    fmin <- function(Abar.min) activityfun2(Abar.min,i)
    fmax <- function(Abar.max) activityfun2(Abar.max,i)
    maxiter <- 1000
    myiter <- 0
    while(sign(fmin(Abar.min))==sign(fmax(Abar.max))) {
      # spread out the values until we have opposite signs 
      diff <- Abar.max - Abar.min
      Abar.min <- Abar.min - diff/2
      Abar.max <- Abar.max + diff/2
      myiter <- myiter + 1
      if(myiter==maxiter) stop(paste('diagram: i tried it',maxiter,
        'times but can\'t make it work :< ... giving up on Abar'))
    }
    # how badly we want the right answer, might
    # have to be adjusted in certain cases
    Atol <- 1e-5
    # find the affinity that gives us the right amount of stuff
    Abar.new <- uniroot(activityfun2,interval=c(Abar.min,Abar.max),i=i,tol=Atol)$root
    # test: did we converge? this shouldn't be required,
    # as uniroot would spit out warnings or errors ... but it doesn't, 
    # even when the tolerance isn't reached by a factor of 100?!
    shouldbezero <- activityfun2(A=Abar.new,i=i)
    if(abs(shouldbezero) > Atol*100) 
      cat(paste('diagram: poor convergence in step ',i,' (remainder in logact of ',
        shouldbezero,').\n',sep=''))
    # and save the activities of the species
    for(j in 1:length(A)) A[[j]][i] <- activityfun(Abar.new,j,i)
  }
  # replace dimensions
  for(i in 1:length(A)) {
    dim(A[[i]]) <- Adim
    dim(av[[i]]) <- Adim
  }
  return(A)
}

