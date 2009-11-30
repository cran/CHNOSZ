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
  ylab=NULL,xlab=NULL,do.plot=TRUE,as.residue=FALSE,mam=TRUE,group=NULL,
  bg=par("bg"),side=1:4) {

  # store input values
  aval <- affinity$values
  asl <- affinity$species$logact
  # number of species possible
  nspecies.orig <- nspecies <- length(aval)
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
  # the property we're plotting
  if(is.null(property)) {
    property <- affinity$property
    if(property=='A') property <- 'affinity'
  }
  # consider a different number of species if we're grouping them together
  if(property=='affinity' & (!mam | nd==0 | nd==1) & !is.null(group)) {
    nspecies.orig <- ngroup <- length(group) 
    ispecies <- 1:ngroup
  } else {
    ngroup <- nspecies
  }
  # handle line type/width/color arguments
  if(do.plot) {
    if(is.null(lty)) lty <- 1:ngroup
    if(is.null(lwd)) lwd <- rep(1,ngroup)
    if(length(lty)!=ngroup) lty <- rep(lty,length.out=ngroup)
    if(length(lwd)!=ngroup) lwd <- rep(lwd,length.out=ngroup)
    if(missing(col.names)) col.names <- col
    if(length(col.names) != ngroup) col.names <- rep(col.names,length.out=ngroup)
    # figure out colors
    if(missing(color)) {
      if(nd==2) {
        if( add | names(dev.cur())=='postscript' ) color <- NULL
        else color <- 'heat'
      } else color <- 'black'
    }
    if(!is.null(color)) {
      color <- rep(color,length.out=ngroup)
      if(is.character(color[1])) {
        # 20091027: take colors as a subset of the original number of species
        if(color[1]=='rainbow') color <- rainbow(nspecies.orig)[ispecies]
        else if(color[1]=='heat') color <- heat.colors(nspecies.orig)[ispecies]
      }
    }
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
  # 'balance' how to balance reactions
  # (name of basis species, 'PBB', 1, or NULL for auto-select)
  # (or vector to set nbalance)
  # 'nbalance' vector of coefficients of the balanced component 
  isprotein <- grep('_',as.character(affinity$species$name))
  if(missing(balance) & length(isprotein)==nspecies)
    balance <- 'PBB'
  # 20091119 make residue=TRUE the default for protein reactions
  if(missing(residue) & length(isprotein)==nspecies) 
    residue <- TRUE
  if(is.null(balance)) {
    ib <- which.balance(affinity$species)
    if(!is.null(ib)) {
      balance <- rownames(basis())[ib[1]]
      nbalance <- affinity$species[,ib[1]]
      cat(paste('diagram: immobile component is',balance,'\n'))
    } else {
      balance <- 1
      nbalance <- rep(1,nspecies)
      cat(paste('diagram: immobile components is formulas as written\n'))
    }
  } else {
   if(is.character(balance[1])) {
      # is the balance is "PBB" for protein backbone
      if(identical(balance,'PBB')) {
        if(length(isprotein) != nspecies)
          stop('diagram: PBB was the requested balance, but ',
            affinity$species$name[-isprotein],' is not protein.')
        nbalance <- as.numeric(protein.length(as.numeric(names(aval))))
        cat(paste('diagram: immobile component is protein backbone group\n'))
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
  if(length(nbalance) < 100) cat(paste('diagram: conservation coefficients are ',c2s(nbalance),'\n',sep=''))
  # 2D diagram mam and residue both FALSE can lead to long calculation times
  if(nd==2 & !mam & !residue) {
    warning("expect long calculation for 2D diagram with mam and residue both FALSE",immediate.=TRUE)
  }
  # rewrite the balance vector for residue reactions
  if(residue) {
    if(any(nbalance < 0)) stop("negative balance prohibits using residue equivalents")
    # affinities divided by balance
    for(i in 1:nspecies) aval[[i]] <- aval[[i]]/nbalance[i]
    oldlogact <- affinity$species$logact
    affinity$species$logact <- affinity$species$logact + log10(nbalance)
    oldbalance <- nbalance
    nbalance <- rep(1,length(nbalance))
    cat(paste('diagram: using residue equivalents\n'))
  }
  # get labels for the plot
  if(!is.null(names)) {
    if(is.na(names[1])) {
      # if names are missing ... make them up
      # properties of basis species or reactions?
      if(affinity$property %in% c('G.basis','logact.basis')) names <- rownames(affinity$basis)
      else {
        if(!is.null(group)) names <- names(group)
        else names <- as.character(affinity$species$name)
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
  # 20091116 replaced plot.curve with plot.line; different
  # name, same functionality, *much* faster
  plot.line <- function(out,xlim,ylim,dotted,col,xrange) {
    # plot boundary lines between predominance fields
    vline <- function(out,ix) {
      ny <- nrow(out)
      xs <- rep(ix,ny*2+1)
      if(0 %in% (ix%%dotted)) 
        return(list(xs=xs,ys=rep(NA,length(xs))))
      ys <- c(rep(ny:1,each=2),0)
      y1 <- out[,ix]
      y2 <- out[,ix+1]
      iy <- which(y1==y2)
      ys[iy*2] <- NA
      return(list(xs=xs,ys=ys))
    }
    hline <- function(out,iy) {
      nx <- nrow(out)
      ny <- ncol(out)
      ys <- rep(ny-iy,nx*2+1)
      if(0 %in% (iy%%dotted)) 
        return(list(xs=rep(NA,length(ys)),ys=ys))
      xs <- c(0,rep(1:nx,each=2))
      x1 <- out[iy,]
      x2 <- out[iy+1,]
      ix <- which(x1==x2)
      xs[ix*2] <- NA
      return(list(xs=xs,ys=ys))
    }
    clipfun <- function(z,zlim) {
      if(zlim[2] > zlim[1]) {
        z[z>zlim[2]] <- zlim[2]
        z[z<zlim[1]] <- zlim[1]
      } else {
        z[z>zlim[1]] <- zlim[1]
        z[z<zlim[2]] <- zlim[2]
      }
      return(z)
    }
    rx <- (xlim[2] - xlim[1]) / (ncol(out) - 1)
    ry <- (ylim[2] - ylim[1]) / (nrow(out) - 1)
    # vertical lines
    xs <- ys <- NA
    for(ix in 1:(ncol(out)-1)) {
      vl <- vline(out,ix)
      xs <- c(xs,vl$xs,NA)
      ys <- c(ys,vl$ys,NA)
    }
    xs <- xlim[1] + (xs - 0.5) * rx
    ys <- ylim[1] + (ys - 0.5) * ry
    ys <- clipfun(ys,ylim)
    lines(xs,ys,col=col)
    # horizontal lines
    xs <- ys <-NA
    for(iy in 1:(nrow(out)-1)) {
      hl <- hline(out,iy)
      xs <- c(xs,hl$xs,NA)
      ys <- c(ys,hl$ys,NA)
    }
    xs <- xlim[1] + (xs - 0.5) * rx
    ys <- ylim[1] + (ys - 0.5) * ry
    xs <- clipfun(xs,xlim)
    if(!is.null(xrange)) xs <- clipfun(xs,xrange)
    lines(xs,ys,col=col)
  }
  ### label plot function
  # calculate coordinates for field labels
  plot.names <- function(out,xs,ys,names) {
    ll <- ngroup
    lx <- numeric(ll); ly <- numeric(ll); n <- numeric(ll)
    for(j in nrow(out):1) {
      # 20091116 for speed, loop over ngroup instead of k (columns)
      for(i in 1:ngroup) {
        k <- which(out[j,]==i)
        if(length(k)==0) next
        lx[i] <- lx[i] + sum(xs[k])
        ly[i] <- ly[i] + length(k)*ys[nrow(out)+1-j]
        n[i] <- n[i] + length(k)
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
    text(lx,ly,labels=names[is],cex=cex.names,col=col.names[is])
  }
  ### property description
  axis.title <- function(property,suffix="") {
    if(property=='A') return(as.expression(substitute(italic(bold("A"))/2.303*italic(R)*italic(T)~~x,list(x=suffix))))
    else if(property=="affinity") return(as.expression(substitute(log*italic(a)~~x,list(x=suffix))))
    else if(property=='logact.basis') return(paste('logQ*',suffix))
    else if(!property %in% c("logK","logQ")) return(paste("Delta",property,suffix))
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
      A <- abundance.new(Astar,nbalance,logatotal)
    }
    else {
      for(j in 1:length(Astar)) Astar[[j]] <- Astar[[j]] + affinity$species$logact[j]
      A <- abundance.old(Astar,aval,nbalance,logatotal)
    }
    # if we rewrote the formation reactions per residue,
    # get back to activities of species
    if(residue & !as.residue) for(j in 1:length(A)) A[[j]] <- A[[j]] - log10(oldbalance[j])
    # groups the species together into meta-species 20090524
    if(!is.null(group)) {
      # store our dimensions
      mydim <- dim(A[[1]])
      # set up the output
      B <- A[1:length(group)]
      for(i in 1:length(B)) {
        B[[i]] <- as.numeric(A[[i]])
        B[[i]][] <- 0
      }
      # add up the activities
      for(i in 1:length(group)) {
        for(j in 1:length(group[[i]])) {
          B[[i]] <- B[[i]] + 10^as.numeric(A[[group[[i]][j]]])
        }
      }
      # return to logarithms and replace dimensions
      for(i in 1:length(B)) {
        # do we want to divide by the number of representative for 
        # each group? probably not 20091017
        #B[[i]] <- log10(B[[i]]/length(group[[i]]))
        B[[i]] <- log10(B[[i]])
        dim(B[[i]]) <- mydim
      }
      A <- B
    }
  }
  ### 0-D properties of species or reactions for single set of conditions
  if(nd==0) {
    if(do.plot) {
      mgp <- par("mgp")
      mgp[1] <- yline
      if(property=="affinity") {
        # the logarithms of activities
        v <- numeric()
        # alpha: plot degree of formation instead of logact
        # scale the activities to sum=1  ... 20091017
        if(alpha) {
          a <- numeric()
          for(j in 1:length(A)) a <- c(a,A[[j]])
          loga.sum <- log10(sum(10^a))
          loga <- 0; a <- loga - loga.sum
          for(j in 1:length(A)) A[[j]] <- 10^(A[[j]] + a)
          ylab <- "alpha"
          if(missing(title)) title <- "alpha"
        } else {
          ylab <- axis.title(property)
          if(missing(title)) title <- "logarithm of activity"
        }
        for(i in 1:length(A)) v <- c(v,A[[i]])
        barplot(v,names.arg=names,ylab=ylab,mgp=mgp,cex.names=cex.names)
        if(missing(title)) title <- title
        title(main=title)
      } else {
        # the value of a property like affinity, standard gibbs energy etc
        for(i in 1:nspecies) aval[[i]] <- aval[[i]]/nbalance[i]
        v <- as.numeric(aval)
        barplot(v,names.arg=names,ylab=axis.title(affinity$property),mgp=mgp,cex.names=cex.names)
        if(!is.null(title)) title(main=title)
        else title(main=as.expression(axis.title(affinity$property,balance.title(balance))))
      }
    }
    if(property=="affinity") return(invisible(list(logact=A)))
    else return(invisible(aval))
  }
  ### 1-D (property or speciation) diagram
  if(nd==1) {
    xvalues <- seq(affinity$xlim[1],affinity$xlim[2],length.out=affinity$xlim[3])
    if(is.null(xlab)) xlab <- axis.label(affinity$xname,as.character(affinity$basis$state[match(affinity$xname,rownames(affinity$basis))]))
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
        if(missing(color)) color <- rep(par('col'),ngroup)
        if(ylog) {
          if(is.null(ylab)) ylab <- as.expression(quote(log~italic(a)))
          if(!add) thermo.plot.new(xlim=affinity$xlim[1:2],
            ylim=ylim,xlab=xlab,ylab=ylab,cex=cex,cex.axis=cex.axis,mar=mar,
            yline=yline,side=side)
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
            ylim=ylim,xlab=xlab,ylab=ylab,cex=cex,cex.axis=cex.axis,mar=mar,
            yline=yline,side=side)
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
      if(missing(ylim)) ylim <- c(min(v),max(v))
      if(!add) thermo.plot.new(xlim=affinity$xlim[1:2],ylim=ylim,xlab=xlab,
        ylab=axis.title(affinity$property),cex=cex,mar=mar,yline=yline,side=side)
      for(i in 1:length(A)) 
        lines(xvalues,as.numeric(aval[[i]]),col=color[i],lty=lty[i])
      if(!is.null(title)) title(main=title)
      else title(main=axis.title(affinity$property,balance.title(balance)))
    }
    if(affinity$property=='logact.basis') names <- rownames(affinity$basis)
    # 20090826: use bg argument
    if(do.plot & !add & !is.null(legend.x)) legend(x=legend.x,lty=lty,legend=names,col=color,bg=bg,cex=cex.names,lwd=lwd)
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
      for(k in 1:ngroup) values[[k]] <- myvalues[[k]][j,]
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
          xlab=xlab,ylab=ylab,cex=cex,cex.axis=cex.axis,mar=mar,yline=yline,side=side)
      }
      # colors and curves
      if(!is.null(color)) plot.color(xs,ys,out,color,ngroup)
      if(!is.null(names)) plot.names(out,xs,ys,names)
      if(!is.null(dotted)) plot.line(out,affinity$xlim[1:2],
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

abundance.new <- function(Astar,nbalance,thisloga) {
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
        shouldbezero,')\n',sep=''))
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


##
## everything below is related to diversity calculations
## 20090415 jmd
##

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
  return(list(x=x.out,y=y.out))
}

extremes <- function(z,what="minimum") {
  # takes a matrix, returns the y as f(x) and x as f(y)
  # trajectories of the extreme
  y <- x <- numeric()
  xres <- ncol(z)
  yres <- nrow(z)
  if(what %in% c("minimum","min")) {
    for(i in 1:xres) y <- c(y,which.min(z[,i]))
    for(i in 1:yres) x <- c(x,which.min(z[i,]))
  } else {
    for(i in 1:xres) y <- c(y,which.max(z[,i]))
    for(i in 1:yres) x <- c(x,which.max(z[i,]))
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
    contour(xs,ys,t(H),add=TRUE)
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

strip <- function(affinity,ispecies=NULL,col=NULL,ns=NULL,
  xticks=NULL,ymin=-0.2,xpad=1) {
  # make strip chart(s) showing the degrees of formation
  # of species as color bars of varying widths
  # extracted from bison/plot.R 20091102 jmd
  # figure out defaults
  a <- affinity
  xlim <- a$xlim[1:2]
  xlab <- axis.label(a$xname)
  if(is.null(ispecies)) ispecies <- list(1:nrow(a$species))
  if(!is.list(ispecies)) ispecies <- list(ispecies)
  if(!is.null(ns) & !is.list(ns)) ns <- list(ns)
  if(is.null(col)) col <- rainbow(length(ispecies[[1]]))
  # how many strip charts on this plot
  # determined by the length of the ispecies list
  nstrip <- length(ispecies)
  # start up the plot
  plot.xlim <- c(xlim[1]-xpad,xlim[2]+xpad)
  ymax <- nstrip+0.3
  thermo.plot.new(xlim=plot.xlim,ylim=c(ymin,ymax),xlab=xlab,ylab="",
      side=c(1,3),mar=par('mar'),do.box=FALSE)
  if(!is.null(xticks)) {
    # mark the positions of the sites on the x-axis
    for(i in 1:5) lines(rep(xticks[i],2),c(ymin,ymin+0.1),lwd=6,col=col[i])
    for(i in 1:5) lines(rep(xticks[i],2),c(ymax,ymax-0.1),lwd=6,col=col[i])
  }
  for(j in 1:nstrip) {
    # get the degrees of formation
    d <- diagram(a,residue=TRUE,do.plot=FALSE,alpha=TRUE,ispecies=ispecies[[j]])
    # depict the relative stabilities of the proteins as color bars
    # make vertical color bars sizes proportional to abundances
    xs <- seq(xlim[1],xlim[2],length.out=length(d$logact[[1]]))
    # total height of the stack
    ly <- 0.5  
    # where to start plotting on the y-axis
    y0 <- j-ly
    for(i in 1:length(d$logact[[1]])) {
      # create a vector of abundances
      loga <- numeric()
      for(k in 1:length(ispecies[[j]])) loga <- c(loga,d$logact[[k]][i])
      loga.order <- order(loga,decreasing=FALSE)
      # plot the bars, least abundant first 
      dy <- y0    # keep track of the height of the stack
      for(k in 1:length(ispecies[[j]])) {
        y1 <- dy
        y2 <- y1 + loga[loga.order[k]] * ly
        dy <- y2
        # note: lwd should be lowered here for very high-resolution plots
        lines(c(xs[i],xs[i]),c(y1,y2),col=col[loga.order[k]],lwd=3,lend="butt")
      }
    }
    # label the color bar
    text(xlim[1],j-ly*1.4,names(ispecies[j]),adj=0,cex=0.7)
    # add inset plot showing the relative numbers of species
    if(!is.null(ns)) {
      ys1 <- y0 - 0.85*ly
      ys2 <- y0 - 0.15*ly
      xmax <- xlim[2]
      xmin <- xlim[1] + 17/22*diff(xlim)
      xss <- seq(xmin,xmax,length.out=length(ispecies[[j]]))
      yss <- numeric()
      for(i in 1:length(ispecies[[j]])) yss <- c(yss,ys1+(ys2-ys1)*(ns[[j]][i]/max(ns[[j]])))
      points(xss,yss,pch=20,cex=0.5)
      lines(xss,yss)
    }
  }
}

