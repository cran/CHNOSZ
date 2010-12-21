# CHNOSZ/diagram.R
# plot predominance or activity diagrams 
# from affinities of formation reactions
# 20061023 jmd

diagram <- function(affinity,what="logact",ispecies=NULL,balance=NULL,
  names=NA,color=NA,add=FALSE,dotted=0,cex=par('cex'),col=par('col'),
  pe=TRUE,pH=TRUE,ylog=TRUE,main=NULL,cex.names=1,legend.x='topright',
  lty=NULL,col.names=par('fg'),cex.axis=par('cex'),logact=NA,lwd=par('lwd'),
  alpha=FALSE,mar=NULL,residue=FALSE,yline=par('mgp')[1]+1,xrange=NULL,
  ylab=NULL,xlab=NULL,do.plot=TRUE,as.residue=FALSE,mam=TRUE,group=NULL,
  bg=par("bg"),side=1:4,xlim=NULL,ylim=NULL) {

  # store input values
  aval <- affinity$values
  asl <- affinity$species$logact
  # number of species possible
  nspecies.orig <- nspecies <- length(aval)
  # number of dimensions (T, P or chemical potentials that are varied)
  mydim <- dim(aval[[1]])
  nd <- length(mydim)
  if(nd==1) if(mydim==1) nd <- 0
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
    if(length(ispecies)==0) {
      cn <- caller.name()
      stop(paste("the length of ispecies is zero ( from",cn,")"))
    }
    cat(paste('diagram: using',length(ispecies),'of',nspecies,'species.\n'))
    affinity$species <- affinity$species[ispecies,]
    aval <- aval[ispecies]
  }
  # number of species that are left
  nspecies <- length(aval)
  # what property we're plotting (jmd 20101014)
  do.loga.equil <- FALSE
  if(affinity$property=="A") {
    if(what %in% rownames(affinity$basis)) {
      # to calculate the loga of basis species at equilibrium
      do.loga.equil <- TRUE
      residue <- FALSE
      balance <- 1
    }
  } else {
    # it might be the values returned by affinity()
    if(missing(what)) what <- affinity$property
    else if(what != affinity$property) stop(paste("can't plot",what,"since the property in 'a' is",affinity$property))
  }

  # consider a different number of species if we're grouping them together
  if(what=='affinity' & (!mam | nd==0 | nd==1) & !is.null(group)) {
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
  # handle legend argument (turn off legend
  # in 1D diagrams if too many species)
  if(nd==1 & ngroup > 10 & missing(legend.x)) legend.x <- NULL
  # covert from activities of proton or electron to pH or pe
  if(affinity$xname=='H+' & pH) { affinity$xname <- 'pH'; affinity$xlim[1:2] <- -affinity$xlim[1:2]; affinity$xvals <- -affinity$xvals }
  if(affinity$xname=='e-' & pe) { affinity$xname <- 'pe'; affinity$xlim[1:2] <- -affinity$xlim[1:2]; affinity$xvals <- -affinity$xvals }
  if(affinity$yname=='H+' & pH) { affinity$yname <- 'pH'; affinity$ylim[1:2] <- -affinity$ylim[1:2]; affinity$yvals <- -affinity$yvals }
  if(affinity$yname=='e-' & pe) { affinity$yname <- 'pe'; affinity$ylim[1:2] <- -affinity$ylim[1:2]; affinity$yvals <- -affinity$yvals }
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
  axis.title <- function(what,suffix="") {
    if(what=="A") return(as.expression(substitute(italic(bold("A"))/2.303*italic(R)*italic(T)~~x,list(x=suffix))))
    else if(what=="logact") return(as.expression(substitute(log*italic(a)~~x,list(x=suffix))))
    else if(what=='logact.basis') return(paste('logQ*',suffix))
    else if(!what %in% c("logK","logQ")) return(paste("Delta",what,suffix))
    else return(what)
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
  if(what=="logact" & (!mam | nd==0 | nd==1) ) {
    # compute the activities of species
    # logarithm of total activity of the balanced component
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
      A <- logact.mb(Astar,nbalance,logatotal)
    }
    else {
      for(j in 1:length(Astar)) Astar[[j]] <- Astar[[j]] + affinity$species$logact[j]
      A <- logact.react(Astar,aval,nbalance,logatotal)
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
  } else if(do.loga.equil) {
    # calculate the logarithm of activity of a basis species
    # at equilibrium
    ibasis <- match(what,rownames(affinity$basis))
    # the reference logact
    AV <- aval
    loga.basis <- as.numeric(affinity$basis$logact[ibasis])
    if(is.na(loga.basis)) stop(paste("the logact for basis species",what,"is not numeric"))
    nu.basis <- affinity$species[,ibasis]
    # the logact where affinity = 0
    A <- lapply(1:length(AV),function(x) {loga.basis - AV[[x]]/nu.basis[x]})
  }

  ### 0-D properties of species or reactions for single set of conditions
  if(nd==0) {
    if(do.plot) {
      mgp <- par("mgp")
      mgp[1] <- yline
      if(what=="logact") {
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
          if(missing(main)) main <- "alpha"
        } else {
          ylab <- axis.title(what)
          if(missing(main)) main <- "logarithm of activity"
        }
        for(i in 1:length(A)) v <- c(v,A[[i]])
        barplot(v,names.arg=names,ylab=ylab,mgp=mgp,cex.names=cex.names,col=col)
        if(missing(main)) main <- main
        title(main=main)
      } else if(do.loga.equil) {
        # the equilibrium logact of a basis species
        v <- as.numeric(A)
        barplot(v,names.arg=names,ylab=axis.label(what),mgp=mgp,cex.names=cex.names,col=col)
        if(!is.null(main)) title(main=main)
      } else {
        # the value of a property like affinity, standard gibbs energy etc
        for(i in 1:nspecies) aval[[i]] <- aval[[i]]/nbalance[i]
        v <- as.numeric(aval)
        barplot(v,names.arg=names,ylab=axis.title(what),mgp=mgp,cex.names=cex.names,col=col)
        if(!is.null(main)) title(main=main)
        else title(main=as.expression(axis.title(what,balance.title(balance))))
      }
    }
    if(what=="logact" | do.loga.equil) return(invisible(list(basis=affinity$basis,species=affinity$species,
      T=affinity$T,P=affinity$P,xname=affinity$xname,xlim=affinity$xlim,yname=affinity$yname,
      ylim=affinity$ylim,logact=A)))
    else return(invisible(list(aval)))
  }
  ### 1-D (property or speciation) diagram
  if(nd==1) {
    xvalues <- affinity$xvals
    if(is.null(xlab)) xlab <- axis.label(affinity$xname,as.character(affinity$basis$state[match(affinity$xname,rownames(affinity$basis))]))
    if(what == 'logact') {
      # alpha: plot degree of formation instead of logact
      if(alpha) {
        ylog <- FALSE
        # scale the activities to sum=1
        for(i in 1:length(A[[1]])) {
          a <- numeric()
          for(j in 1:length(A)) a <- c(a,A[[j]][i])
          loga.sum <- log10(sum(10^a))
          loga <- 0; a <- loga - loga.sum
          for(j in 1:length(A)) A[[j]][i] <- A[[j]][i] + a
        }
      }
      aval <- A
      if(ylog) {
        if(is.null(ylab)) ylab <- as.expression(quote(log~italic(a)))
      } else {
        # are we plotting alphas or logactivities?
        if(is.null(ylab)) {
          if(alpha) ylab <- as.expression(quote(alpha))
          else ylab <- as.expression(quote(a[i]))
        }
        for(i in 1:length(A)) aval[[i]] <- 10^A[[i]]
      }
    } else if(do.loga.equil) {
      # plot equilibrium logarithms of activities of basis species
      aval <- A
      if(is.null(ylab)) ylab <- axis.label(what)
    } else {
      # plot metastable equilibrium logarithms of activities of species or values of properties
      A <- aval
      if(is.null(ylab)) ylab <- axis.title(what)
      for(i in 1:nspecies) aval[[i]] <- aval[[i]]/nbalance[i]
    }
    # now make the plot
    if(do.plot) {
      if(!add) {
        if(is.null(xlim)) xlim <- affinity$xlim[1:2]
        # automatically get range for y-axis
        # use only those points that are in the x-range
        if(is.null(ylim)) {
          isx <- xvalues >= min(xlim) & xvalues <= max(xlim)
          xfun <- function(x) x[isx]
          myval <- sapply(aval,xfun)
          ylim <- extendrange(myval)
        }
        thermo.plot.new(xlim=xlim,ylim=ylim,xlab=xlab,
          ylab=ylab,cex=cex,mar=mar,yline=yline,side=side)
      }
      if(length(col) < ngroup) col <- rep(col,ngroup)
      for(i in 1:length(A)) lines(xvalues,as.numeric(aval[[i]]),col=col[i],lty=lty[i],lwd=lwd[i])
      # label the plot
      if(what=="logact.basis") names <- rownames(affinity$basis)
      # 20090826: use bg argument
      if(do.plot & !add & !is.null(legend.x)) legend(x=legend.x,lty=lty,legend=names,col=col,bg=bg,cex=cex.names,lwd=lwd)
      # add a title
      if(!missing(main)) title(main=main)
    }
    if(alpha) for(i in 1:length(A)) A[[i]] <- 10^A[[i]]
    # 20090324 return list with single element 'logact'
    out <- list(basis=affinity$basis,species=affinity$species,
      T=affinity$T,P=affinity$P,xname=affinity$xname,xlim=affinity$xlim,yname=affinity$yname,
      ylim=affinity$ylim)
    if(what=="logact" | do.loga.equil) out <- c(out,list(logact=A))
    else out <- c(out,list(aval=aval))
    return(invisible(out))
  }

  ### 2-D predominance diagram aka equal activity diagram
  if(what=="logact") {
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
      # mam=FALSE: use activities from speciation calculation
      myvalues <- A
    }
    # only make a 2-D plot
    if(nd==2) {
      # out: the index of the predominant species
      out <- myvalues[[1]]
      for(j in 1:nrow(out)) {
        values <- list()
        for(k in 1:ngroup) values[[k]] <- myvalues[[k]][j,]
        out[j,] <- which.pmax(values,na.rm=TRUE)
      }
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
        # put out in the right order for image() etc
        o <- out
        for(i in 1:ncol(out)) o[,i] <- out[,ncol(out)+1-i]
        out <- t(o)
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
    } else {
      #cat(paste('diagram: 2-D plot of',property,'not available\n'))
      if(mam) return(invisible(list(basis=affinity$basis,species=species,T=affinity$T,P=affinity$P,
        xname=affinity$xname,xlim=affinity$xlim,yname=affinity$yname,ylim=affinity$ylim,aval=aval)))
      else return(invisible(list(basis=affinity$basis,species=species,T=affinity$T,P=affinity$P,xname=affinity$xname,
        xlim=affinity$xlim,yname=affinity$yname,ylim=affinity$ylim,aval=aval,logact=A)))
    }
  } else {
    # if we're not calculating predominances we can only make a contour plot for properties
    # of a single reaction
    if(!do.loga.equil) A <- aval
    if(length(A)!=1) warning(paste("can't make contour plot of",what,"for more than one reaction. suggestion: select a single one with 'ispecies'"))
    else if(do.plot & nd==2) {
      if(!add) {
        xstate=as.character(affinity$basis$state[match(affinity$xname,rownames(affinity$basis))])
        ystate=as.character(affinity$basis$state[match(affinity$yname,rownames(affinity$basis))])
        if(is.null(xlab)) xlab <- axis.label(as.character(affinity$xname),xstate)
        if(is.null(ylab)) ylab <- axis.label(as.character(affinity$yname),ystate)
        thermo.plot.new(xlim=affinity$xlim[1:2],ylim=affinity$ylim[1:2],
          xlab=xlab,ylab=ylab,cex=cex,cex.axis=cex.axis,mar=mar,yline=yline,side=side)
        if(missing(main)) main <- axis.label(what)
        title(main=main)
      }
      xs <- seq(affinity$xlim[1],affinity$xlim[2],length.out=affinity$xlim[3])
      ys <- seq(affinity$ylim[1],affinity$ylim[2],length.out=affinity$ylim[3])
      contour(xs,ys,A[[1]],add=TRUE,col=col,lty=lty,lwd=lwd,labcex=cex)
    }
    return(invisible(list(basis=affinity$basis,species=species,T=affinity$T,P=affinity$P,xname=affinity$xname,
      xlim=affinity$xlim,yname=affinity$yname,ylim=affinity$ylim,aval=aval,logact=A)))
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
    # make vertical color bars sizes proportional to activities
    xs <- seq(xlim[1],xlim[2],length.out=length(d$logact[[1]]))
    # total height of the stack
    ly <- 0.5  
    # where to start plotting on the y-axis
    y0 <- j-ly
    for(i in 1:length(d$logact[[1]])) {
      # create a vector of activities
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

