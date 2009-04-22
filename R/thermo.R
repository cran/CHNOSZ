# CHNOSZ/thermo.R
# Copyright (C) 2006-2009 Jeffrey M. Dick
# utility functions for the CHNOSZ package
# speciate/thermo.R 20051021 jmd

eos.args <- function(eos='',property=NULL,T=NULL,P=NULL) {

  # the available properties for supcrt, probably
  props <- c('G','H','S','Cp','V','kT','E')
  if(eos=='water') {
    # things we also get with water
    #props <- c(colnames(thermo$water)[4:length(colnames(thermo$water))])
    props <- c(props,'A','U','Cv','Psat','rho','Q','X','Y','epsilon','w')
    # they keep on coming: things we also get with SUPCRT92
    if(length(agrep(tolower(thermo$opt$water),'supcrt9',max=0.3))>0)
      props <- c(props,'Z','visc','tcond','tdiff','Prndtl','visck','albe','daldT')
    else 
      props <- c(props,'P','N','UBorn','de.dT','de.dP')
  }

  # so others can find out what we're about
  if(is.null(property)) return(props)

  # default: all properties and contributions
  if(is.null(property)) property <- props
  # the lowercase equivalent of the argument
  prop <- tolower(property)
  # and returns its lower- and upper- case equivalents
  # (prop and Prop) and the available properties
  Prop <- props[match(prop,tolower(props))]
  #contrib <- tolower(contrib)
  # error: not a property
  notprop <- ! prop %in% tolower(props)
  if(TRUE %in% notprop) 
    stop('thermo.args: properties ',c2s(prop[notprop]),' not in ',c2s(props),'.\n')
  # return arguments
  return(list(props=props,prop=prop,Prop=Prop))

}

TP.args <- function(T=NULL,P=NULL) {
  if(!is.null(P)) {
    if(P[1]=='Psat') {
      P <- water('Psat',T,P=NULL)
      P <- P[,1]
      # water.SUPCRT92 issues its own warnings about 
      # exceeding Psat's temperature limit
      if(length(agrep(tolower(thermo$opt$water),'supcrt9',max=0.3))==0)
        if(length(which(is.na(P)))>0) 
          warning('TP.args: NAs in Psat (likely T > Tc where Tc = 647.096 K)',call.=FALSE)
    }
  }
  if(length(P) < length(T)) P <- rep(P, length.out=length(T))
  else if(length(T) < length(P)) T <- rep(T, length.out=length(P))
  # something we do here so the SUPCRT water calculations work
  T[T==273.15] <- 273.16
  return(list(T=T,P=P))
}

state.args <- function(state=NULL) {
  if(is.null(state) | is.numeric(state[1])) return(state)
  # normalize state arguments
  for(i in 1:length(state)) {
    if(tolower(state[i])=='a') state[i] <- 'aq'
    if(tolower(state[i])=='c') state[i] <- 'cr'
    if(tolower(state[i])=='g') state[i] <- 'gas'
    if(tolower(state[i])=='l') state[i] <- 'liq'
  }
  return(state)
}

aminoacids <- function(seq=NULL,nchar=1) {
  # if no sequence is given just return the
  # abbreviations of the amino acids
  # these are in the same order as in thermo$protein
  aa <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
  aa.3 <- c('Ala','Cys','Asp','Glu','Phe','Gly','His','Ile','Lys','Leu',
            'Met','Asn','Pro','Gln','Arg','Ser','Thr','Val','Trp','Tyr')
  aa.NA <- c("alanine","cysteine","aspartic acid","glutamic acid","phenylalanine","glycine",
    "histidine","isoleucine","lysine","leucine","methionine","asparagine","proline","glutamine",
    "arginine","serine","threonine","valine","tryptophan","tyrosine")
  aa.Z <- c("alanine","cysteine","aspartate","glutamate","phenylalanine","glycine",
    "histidine","isoleucine","lysinium","leucine","methionine","asparagine","proline","glutamine",
    "argininium","serine","threonine","valine","tryptophan","tyrosine")
  if(is.null(seq)) {
    if(is.na(nchar)) return(aa.NA)
    if(nchar=="Z") return(aa.Z)
    if(as.numeric(nchar)==1) return(aa)
    if(as.numeric(nchar)==3) return(aa.3)
  }
  # sequences are given as elements of the list seq
  # to count the number of each amino acids in a sequence
  count.aa <- function(aa,seq) sum(seq==aa)
  count.i <- function(i,seq) as.numeric(lapply(aa,count.aa,strsplit(toupper(seq[i]),"")[[1]]))
  # count amino acids in each sequence
  a <- t(as.data.frame(mylapply(1:length(seq),count.i,seq),optional=TRUE))
  # clean up row/column names
  colnames(a) <- aa
  rownames(a) <- 1:nrow(a)
  return(a)
}


GHS <- function(species=NULL,DG=NA,DH=NA,S=NA,T=thermo$opt$Tr) {
  # convert among G, H, S, Se
  # calculate Se (entropy of elements)
  if(!is.null(species)) Se <- element(species,'entropy')[,1] else Se <- 0
  if(NA%in%(DG) & NA%in%(DH) & NA%in%(S)) return(Se)
  if(NA%in%(DG) & NA%in%(DH) | NA%in%(DH) & NA%in%(S) | NA%in%(DG) & NA%in%(S))
    stop('please supply zero, two or three of G,H,S.')
  if(NA%in%(DH)) return(DG+T*(S-Se)) else
  if(NA%in%(S)) return((DH-DG)/T+Se) 
  # G, finally
  return(DH-T*(S-Se))
}

element <- function(compound,property=c('mass','entropy')) {
  # calculate the mass and/or entropy of elements in compound(s)
  if(!is.character(compound[1])) stop('please specify a character argument')
  for(j in 1:length(compound)) {
    m <- makeup(compound[j])
    ie <- match(rownames(m),thermo$element$element)
    if('mass' %in% property)
      mass <- data.frame(mass=sum(thermo$element$mass[ie]*m$count))
    if('entropy' %in% property)
      entropy <- data.frame(entropy=sum(thermo$element$s[ie]/thermo$element$n[ie]*m$count))
    if(length(property)>1) {
      p <- get(property[1])
      for(i in 2:length(property)) p <- cbind(p,get(property[i]))
    } else p <- get(property)
    if(j>1) prop <- rbind(prop,p) else prop <- p
  }
  return(prop)
}

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

nuts <- function(units=NULL) {
  # change preferred units or list the current ones
  
  # show the current units, if none are specified
  if(missing(units)) {
    cat(paste('nuts: temperature in ',thermo$opt$T.units,'.\n',sep=''))
    cat(paste('nuts: energy in ',thermo$opt$E.units,'.\n',sep=''))
    cat(paste('nuts: pressure in ',thermo$opt$P.units,'.\n',sep=''))
    return(invisible())
  }
  
  for(i in 1:length(units)) {
    # argument handling
    if(is.character(units[i])) {
      units[i] <- tolower(units[i])
      if(!units[i] %in% c('t','e','p','c','k','cal','j','bar','mpa')) stop('units must be one of T, E, P, C, K, cal, j, bar, MPa.')
    } else stop('units needs to be character.')
    # return the current units of a particular property
    if(units[i]=='t') return(thermo$opt$T.units)
    if(units[i]=='e') return(thermo$opt$E.units)
    if(units[i]=='p') return(thermo$opt$P.units)

    # tests and update
    if(units[i] %in% c('c','k')) {
      if(units[i]=='c') thermo$opt$T.units <<- 'C'
      if(units[i]=='k') thermo$opt$T.units <<- 'K'
      cat(paste('nuts: temperature in ',thermo$opt$T.units,'.\n',sep=''))
    }
    if(units[i] %in% c('j','cal')) {
      if(units[i]=='j') thermo$opt$E.units <<- 'J'
      if(units[i]=='cal') thermo$opt$E.units <<- 'cal'
      cat(paste('nuts: energy in ',thermo$opt$E.units,'.\n',sep=''))
    }
    if(units[i] %in% c('bar','mpa')) {
      if(units[i]=='bar') thermo$opt$P.units <<- 'bar'
      if(units[i]=='mpa') thermo$opt$P.units <<- 'MPa'
      cat(paste('nuts: pressure in ',thermo$opt$P.units,'.\n',sep=''))
    }
  }
}

outvert <- function(value,units) {
  # converts the given value from the given units to
  # those specified in thermo$opt
  units <- tolower(units)
  if(units %in% c('c','k')) {
    if(units=='c' & thermo$opt$T.units=='K') return(convert(value,'k'))
    if(units=='k' & thermo$opt$T.units=='C') return(convert(value,'c'))
  }
  if(units %in% c('j','cal')) {
    if(units=='j' & thermo$opt$E.units=='Cal') return(convert(value,'cal'))
    if(units=='cal' & thermo$opt$E.units=='J') return(convert(value,'j'))
  }
  if(units %in% c('bar','mpa')) {
    if(units=='mpa' & thermo$opt$P.units=='bar') return(convert(value,'bar'))
    if(units=='bar' & thermo$opt$P.units=='MPa') return(convert(value,'mpa'))
  }
  return(value)
}

envert <- function(value,units) {
  # convert values to the specified units
  # from those given in thermo$opt
  if(!is.numeric(value[1])) return(value)
  units <- tolower(units)
  if(units %in% c('c','k','t.units')) {
    if(units=='c' & thermo$opt$T.units=='K') return(convert(value,'c'))
#    if(units=='t.units' & thermo$opt$T.units=='C') return(convert(value,'c'))
    if(units=='k' & thermo$opt$T.units=='C') return(convert(value,'k'))
  }
  if(units %in% c('j','cal','e.units')) {
    if(units=='j' & thermo$opt$T.units=='Cal') return(convert(value,'j'))
#    if(units=='e.units' & thermo$opt$T.units=='J') return(convert(value,'j'))
    if(units=='cal' & thermo$opt$T.units=='J') return(convert(value,'cal'))
  }
  if(units %in% c('bar','mpa','p.units')) {
    if(units=='mpa' & thermo$opt$P.units=='bar') return(convert(value,'mpa'))
#    if(units=='p.units' & thermo$opt$P.units=='MPa') return(convert(value,'mpa'))
    if(units=='bar' & thermo$opt$P.units=='MPa') return(convert(value,'bar'))
  }
  return(value)
}

convert <- function(value,units,T=thermo$opt$Tr,P=thermo$opt$Pr,pH=7,logaH2O=0) {
  # converts value(s) to the specified units

  if(is.null(value)) return(NULL)
  ### argument handling
  if(!is.character(units)) stop(paste('convert: please specify',
    'a character argument for the destination units.\n',
    'possibilities include (G or logK) (C or K) (J or cal) (cm3bar or calories) (Eh or pe)\n',
    'or their lowercase equivalents.\n'),call.=FALSE)
  Units <- units # for the possible message to user
  units <- tolower(units)

  # make everything same length
  makesame <- function(x) {
    maks <- max(as.numeric(sapply(x,length)))
    for(i in 1:length(x)) x[[i]] <- rep(x[[i]],length.out=maks)
    return(x)
  }
  if(!is.matrix(value)) {
    args <- makesame(list(value=value,units=units,T=T,pH=pH,logaH2O=logaH2O))
    # except we only accept one type of destination units
    value <- args$value; units <- args$units[1]; T <- args$T; pH <- args$pH; logaH2O <- args$logaH2O
  }
  # tests and calculations for the specified units
  if(units %in% c('c','k')) {
    CK <- 273.15
    if(units=='k') value <- value + CK
    if(units=='c') value <- value - CK 
  }
  else if(units[1] %in% c('j','cal')) {
    Jcal <- 4.184
    if(units=='j') value <- value * Jcal
    if(units=='cal') value <- value / Jcal
  }
  else if(units %in% c('g','logk')) {
    if(units=='logk') value <- value / (-log(10) * thermo$opt$R * T)
    if(units=='g') value <- value * (-log(10) * thermo$opt$R * T)
  }
  else if(units %in% c('cm3bar','calories')) {
    if(units=='cm3bar') value <- convert(value,'J') * 10
    if(units=='calories') value <- convert(value,'cal') / 10
  }
  else if(units %in% c('eh','pe')) {
    R <- 0.00831470
    F <- 96.4935
    if(units=='pe') value <- value * F / ( log(10) * R * T )
    if(units=='eh') value <- value * ( log(10) * R * T ) / F
  }
  else if(units %in% c('bar','mpa')) {
    barmpa <- 10
    if(units=='mpa') value <- value / barmpa
    if(units=='bar') value <- value * barmpa
  }
  else if(units %in% c('e0','logfo2')) {
    # convert between Eh and logfO2
    supcrt.out <- subcrt(c('H2O','O2','H+','e-'),c(-1,0.5,2,2),T=T,P=P,convert=FALSE)
    if(units=='logfo2') value <- 2*(supcrt.out$out$logK + logaH2O + 2*pH + 2*(convert(value,'pe',T=T)))
    if(units=='e0') value <- convert(( -supcrt.out$out$logK - 2*pH + value/2 - logaH2O )/2, 'Eh',T=T)
  }
  else cat(paste('convert: no conversion to ',Units,' found.\n',sep=''))
  return(value)
}

thermo.axis <- function(lab='x-axis',side=1,line=1.5,cex=par('cex'),lwd=par('lwd'),T=NULL,col=par('col')) {
  # if T isn't NULL, looks like we want make a second
  # oxidation scale corresponding to one already plotted.
  # e.g.,  Eh-pe, Eh-logfO2, or logfO2-Eh
  if(!is.null(T)) {
    usr <- par('usr')
    if(side %in% c(1,3)) lim <- usr[1:2] else lim <- usr[3:4]
    if(length(grep('pe',lab)) > 0) {
      lim <- convert(lim,'pe',T=T)
    } else if(length(grep('O2',lab)) > 0) {
      lim <- convert(lim,'logfO2',T=T)
    } else if(length(grep('Eh',lab)) > 0) {
      lim <- convert(lim,'E0',T=T)
    }
    if(side %in% c(1,3)) usr[1:2] <- lim else usr[3:4] <- lim
    opar <- par(usr=usr)
  }
  if(!is.null(lwd)) {
    ## plot major tick marks and numeric labels
    do.label <- TRUE
    if(missing(cex) & side %in% c(3,4) & is.null(T)) do.label <- FALSE
    at <- axis(side,labels=do.label,tick=TRUE,lwd=lwd,col=col,col.axis=col) 
    ## plot minor tick marks
    # the distance between major tick marks
    da <- abs(diff(at[1:2]))
    # distance between minor tick marks
    di <- da / 4
    if(da %% 2 | !(da %% 10)) di <- da / 5
    # number of minor tick marks
    if(side %in% c(1,3)) {
      ii <- c(1,2) 
      myasp <- par('xaxp')
    } else {
      ii <- c(3,4)
      myasp <- par('yaxp')
    }
    myusr <- par('usr')[ii]
    daxis <- abs(diff(myusr))
    nt <- daxis / di + 1
    ## if nt isn't an integer, it probably
    ## means the axis limits don't correspond
    ## to major tick marks (expect problems)
    ##at <- seq(myusr[1],myusr[2],length.out=nt)
    # start from (bottom/left) of axis?
    bl <- 1
    #if(myasp[2]==myusr[2]) bl <- 2
    # is forward direction (top/right)?
    tr <- 1
    if(xor(myusr[2] < myusr[1] , bl==2)) tr <- -1
    #at <- myusr[bl] + tr * di * seq(0:(nt-1))
    # well all of that doesn't work in a lot of cases,
    # where none of the axis limits correspond to
    # major tick marks. perhaps the following will work
    at <- myusr[1] + tr * di * (0:(nt-1))
    # apply an offset
    axt <- axTicks(side)[1]
    daxt <- (axt - myusr[1])/di
    daxt <- (daxt-round(daxt))*di
    at <- at + daxt
    tcl <- par('tcl') * 0.5
    axis(side,labels=FALSE,tick=TRUE,lwd=lwd,col=col,col.axis=col,at=at,tcl=tcl)
  }

  # rotate labels on side axes
  if(side %in% c(2,4)) las <- 0 else las <- 1
  if(!is.null(lab)) mtext(lab,side=side,line=line,cex=cex,las=las)
  # reset limits if we were plotting a second axis
  if(!is.null(T)) par(opar)
}

thermo.plot.new <- function(xlim,ylim,xlab,ylab,cex=par('cex'),mar=NULL,lwd=par('lwd'),ticks=c(1,2,3,4),mgp=c(1.2,0.3,0),cex.axis=par('cex'),col=par('col'),yline=NULL,axs='i') {
  # 20090324 mar handling: NULL - a default setting; NA - par's setting
  # 20090413 changed mar of top side from 2 to 2.5
  if(is.null(mar)) mar <- c(3,3.5,2.5,1) else if(is.na(mar[1])) mar <- par('mar')
  par(mar=mar,mgp=mgp,tcl=0.3,las=1,xaxs=axs,yaxs=axs,cex=cex,lwd=lwd)
  plot.new()
  plot.window(xlim=xlim,ylim=ylim)
  box()
  # labels
  thermo.axis(xlab,side=1,line=mgp[1],cex=cex.axis,lwd=NULL)
  if(is.null(yline)) yline <- mgp[1]
  thermo.axis(ylab,side=2,line=yline,cex=cex.axis,lwd=NULL)
  # (optional) tick marks
  if(1 %in% ticks) thermo.axis(NULL,side=1,lwd=lwd,col=par('col'))
  if(2 %in% ticks) thermo.axis(NULL,side=2,lwd=lwd,col=par('col'))
  if(3 %in% ticks) thermo.axis(NULL,side=3,lwd=lwd,col=par('col'))
  if(4 %in% ticks) thermo.axis(NULL,side=4,lwd=lwd,col=par('col'))
}

water.lines <- function(xaxis='pH',yaxis='Eh',T=298.15,P='Psat',which=c('oxidation','reduction'),logaH2O=0,lty=2,col=par('fg'),xpoints=NULL) {
  # draw water stability limits
  # if we're on an Eh-pH diagram, or logfO2-pH diagram,
  # or logfO2-T or Eh-T
  # calculate them exactly (nicer looking lines), otherwise 
  # (TODO) add them using affinity() and diagram()
  
  # get the x and y limits from the current plot
  pu <- par('usr')
  xlim <- pu[1:2]
  ylim <- pu[3:4]
  # exact lines
  # warning: Eh calculations are reliable only at a single T
  if(xaxis=='pH' & (yaxis=='Eh' | yaxis=='O2' | yaxis=="pe")) {
    if('reduction' %in% which) {
      logfH2 <- 0
      logK <- subcrt(c('H2O','O2','H2'),c(-1,0.5,1),T=T,P=P,convert=FALSE)$out$logK 
      logfO2 <- 2 * logK - logfH2 + 2 * logaH2O
      if(yaxis=='O2') abline(h=logfO2,lty=lty,col=col) 
      else if(yaxis=="Eh") lines(xlim,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,col=col)
      else if(yaxis=="pe") lines(xlim,convert(convert(logfO2,'E0',T=T,P=P,pH=xlim),"pe",T=T),lty=lty,col=col)
    }
    if('oxidation' %in% which) {
      logfO2 <- 0
      if(yaxis=='O2') abline(h=logfO2,lty=lty,col=col) 
      else if(yaxis=="Eh") lines(xlim,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,col=col)
      else if(yaxis=="pe") lines(xlim,convert(convert(logfO2,'E0',T=T,P=P,pH=xlim),"pe",T=T),lty=lty,col=col)
    }
  } else if(xaxis %in% c('T','P') & yaxis %in% c('Eh','O2') ) {
    #if(xaxis=='T') if(is.null(xpoints)) xpoints <- T
    # 20090212 get T values from plot limits
    # TODO: make this work for T on y-axis too
    if(xaxis=='T') {
      if(missing(T)) {
        xpoints <- seq(xlim[1],xlim[2],length.out=100)
        T <- envert(xpoints,"K")
      }
    }
    if(xaxis=='P') if(is.null(xpoints)) xpoints <- P
    if('oxidation' %in% which) {
      logfO2 <- rep(0,length(xpoints))
      if(yaxis=='Eh') lines(xpoints,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,col=col)
      else lines(xpoints,logfO2,lty=lty,col=col)
    }
    if('reduction' %in% which) {
      logfH2 <- 0
      logK <- subcrt(c('H2O','O2','H2'),c(-1,0.5,1),T=T,P=P,convert=FALSE)$out$logK 
      logfO2 <- 2 * logK - logfH2 + 2 * logaH2O
      if(yaxis=='Eh') lines(xpoints,convert(logfO2,'E0',T=T,P=P,pH=xlim),lty=lty,col=col)
      else lines(xpoints,logfO2,lty=lty,col=col)
    }
  } else {
    # inexact lines
    #
  }
}

label.plot <- function(x,xfrac=0.95,yfrac=0.9,cex=1,paren=TRUE,adj=1) {
  # make a text label e.g., "(a)" in the corner of a plot
  # xfrac, yfrac: fraction of axis where to put label (default top right)
  # paren: put a parenthesis around the text, and italicize it?
  if(paren) x <- as.expression(substitute(group('(',italic(a),')'),list(a=x)))
  pu <- par('usr')
  text( pu[1]+xfrac*(pu[2]-pu[1]), pu[3]+yfrac*(pu[4]-pu[3]), labels=x, cex=cex , adj=adj)
}

axis.label <- function(x,opt=NULL,do.state=TRUE,oldstyle=FALSE,do.upper=FALSE,mol='mol') {
  # make axis labels
  lab <- x
  if(missing(opt)) do.opt <- TRUE else do.opt <- FALSE
  if(!is.null(opt)) if(is.na(opt)) do.opt <- TRUE
  if(lab %in% c('T','P','Eh','pH','pe','logK','IS')) {
    if(lab=='Eh') lab <- paste(lab,'(volt)')
    else if(lab=='T') {
      if(do.opt) T.units <- nuts('T') else T.units <- opt
      if(T.units=='C') lab <- as.expression(quote(list(italic(T),degree*C)))
      else lab <- as.expression(quote(list(italic(T),K)))
    } 
    else if(lab=='P') {
      if(do.opt) P.units <- nuts('P') else P.units <- opt
      if(P.units=='bar') lab <- as.expression(quote(list(italic(P),bar)))
      else lab <- as.expression(quote(list(italic(P),MPa)))
    } 
    else if(lab=='logK') lab <- as.expression(quote(log~italic(K)))
    else if(lab=='IS') lab <- as.expression(quote(list(IS,mol~~kg^{-1})))

    return(lab)
  } else {
    if(is.null(thermo$basis)) rn <- '' else rn <- rownames(basis())
    if(lab %in% rn) {
      if(do.opt) {
        # the default state
        #opt <- thermo$opt$state
        # 20090215: the state as this basis species is defined 
        opt <- as.character(thermo$basis$state)[match(lab,rn)]
      }
      state <- opt
      if(oldstyle) {
        # append (log a) or (log f)
        if(state %in% c('gas')) llab <- '(log f)' else llab <- '(log a)'
        newlab <- paste(lab,llab,sep=' ')
        return(newlab)
      } else {
        # make the subscripts
        labform <- makeup(lab)
        newlab <- ''
        for(i in 1:nrow(labform)) {
            #newlab <- paste(newlab,rownames(labform)[i],sep='')
            #if(labform$count[i]!=1) newlab <- paste(newlab,'[',labform$count[i],']',sep='')
            #if(i!=nrow(labform)) newlab <- paste(newlab,'*',sep='')
            if(rownames(labform)[i] != 'Z') {
              newlab <- substitute(paste(a,b),list(a=newlab,b=rownames(labform)[i]))
              # subscripts under subscripts are getting too small
              #if(labform$count[i]!=1) newlab <- substitute(a[b],list(a=newlab,b=labform$count[i]))
              if(labform$count[i]!=1) newlab <- substitute(a*b,list(a=newlab,b=labform$count[i]))
            } else {
              # for charged species, don't show "Z" but do show e.g. "+2"
              lc <- labform$count[i]
              if(lc > 0) lc <- paste("+",as.character(lc),sep="")
              newlab <- substitute(paste(a,b),list(a=newlab,b=lc))
            }
        }
        llab <- 'f'; slab <- 'g'
        if(state %in% c('aq','cr','liq')) {llab <- 'a'; slab <- state}
        #if(do.state) newlab <- substitute(a[group('(',italic(b),')')],list(a=newlab,b=slab))
        if(do.state) newlab <- substitute(a*group('(',italic(b),')'),list(a=newlab,b=slab))
        llab <- substitute(log*italic(a),list(a=llab))
        newlab <- as.expression(substitute(a[b],list(a=llab,b=newlab)))
        return(newlab)
      }
    } else {
      # a way to make expressions for various properties
      # e.g. axis.label('DG0r','k')
      clab <- s2c(lab)
      wlab <- ''
      doital <- TRUE; dosub <- FALSE
      for(i in 1:length(clab)) {
        blab <- clab[i]
        #if(dosub) blab <- substitute(phantom0[a],list(a=blab))
        # D for Delta
        if(i==1 & blab=='D') wlab <- quote(Delta)
        else if(i > 1 & blab=='0') {
          wlab <- substitute(a*degree,list(a=wlab))
          dosub <- TRUE
        }
        else if(i > 1 & (can.be.numeric(blab) | blab=='P' | dosub)) 
          wlab <- substitute(a[italic(b)],list(a=wlab,b=blab))
        else if(i > 1 & blab==',') {
          wlab <- substitute(a*b,list(a=wlab,b=blab))
          doital <- FALSE
        }
        else {
          if(blab=='A') blab <- substitute(bold(a),list(a=blab))
          if(i > 1) {
            if(blab=='p') wlab <- substitute(a[italic(b)],list(a=wlab,b=toupper(blab)))
            else if(doital) wlab <- substitute(a*italic(b),list(a=wlab,b=blab))
            else wlab <- substitute(a*b,list(a=wlab,b=blab))
          } else wlab <- substitute(italic(b),list(b=blab))
        }
      }
      # now to the nuts
      if(do.state) {
        mytoupper <- function(lab) {
          if(do.upper & lab!='g') return(toupper(lab))
          else return(lab)
        }
        if(clab[1]=='D') clab <- clab[-1]
        ulab <- mytoupper(nuts('E'))
        if(clab[1] %in% c('C','S')) mylab <- substitute(a~~K^-1,list(a=ulab))
        else if(clab[1] == c('V')) mylab <- substitute(a^3,list(a=mytoupper('cm')))
        else if(clab[1] == 'E') mylab <- substitute(a^3~~K^-1,list(a=mytoupper('cm')))
        else mylab <- ulab
        if(!is.null(opt)) {
          if(can.be.numeric(opt)) mylab <- substitute(10^a~~b,list(a=opt,b=mylab))
          else {
            opt <- mytoupper(opt)
            mylab <- substitute(a*b,list(a=opt,b=mylab))
          }
        }
        mylab <- substitute(a~~b^-1,list(a=mylab,b=mytoupper(mol)))
        wlab <- substitute(list(a,b),list(a=wlab,b=mylab))
      }
      return(as.expression(wlab))
    }
  }
}


thermo.postscript <- function(file,family='Helvetica',width=8,height=6,horizontal=FALSE) {
  postscript(onefile=FALSE,horizontal=horizontal,paper='special',width=width,height=height,file=file,family=family)
}

# a function to expand formulas (makeups) to any
# number of elements, i.e. 
# expand.formula(c('C','H','O'),makeup('H2O')) returns c(0,2,1)
expand.formula <- function(elements,makeup) {
  # if the makeup is missing, return an dataframe with
  # zero rows and ncol=number of elements
  if(missing(makeup)) {
    comp0 <- as.data.frame(matrix(rep(NA,length(elements)),nrow=1))[0,,drop=FALSE]
    # with column names of elements
    #if(length(elements) > 1) colnames(comp0) <- elements
    colnames(comp0) <- elements
    return(comp0)
  }
  # the elements that are present
  ematch <- match(elements,rownames(makeup))
  # the exapnding
  ecount <- makeup[ematch,]
  # strip NAs
  ecount[is.na(ecount)] <- 0
  # 
  return(ecount)
}

protein.length <- function(protein) {
  # character, name of protein
  # positive numeric, index of protein in thermo$species
  # negative numeric, index of protein in thermo$protein
  # 20090331 added negative numeric option and multicore support

  # to get the length of the ith protein in thermo$protein
  iprotein2length <- function(iprotein) sum(thermo$protein[iprotein,6:25])
  # to find the index of the named protein
  name2iprotein <- function(name) {
    if(length(grep('_',name))==0) {
      warning(name,' is not the name of a protein.')
      return(NA)
    }
    po <- s2c(as.character(name),sep='_',keep.sep=FALSE)
    ip <- as.character(thermo$protein$protein)==po[1] & 
      as.character(thermo$protein$organism)==po[2]
    ip <- ip[!is.na(ip)]
    ip <- which(ip)[1]
    return(ip)
  }
  # to find the protein index of the numbered species
  ispecies2iprotein <- function(ispecies) {
    name2iprotein(as.character(thermo$species$name)[match(ispecies,thermo$species$ispecies)])
  }
  # to get protein length for iprotein, name or ispecies
  x2length <- function(x) {
    if(can.be.numeric(x)) {
      if(x > 0) iprotein2length(ispecies2iprotein(as.numeric(x)))
      else iprotein2length(-as.numeric(x))
    } else iprotein2length(name2iprotein(x))
  }
  # do it!
  return(as.numeric(mylapply(protein,x2length)))
  
  
}



dPdTtr <- function(x) {
  # calculate dP/dT for a phase transition
  # (argument is index of the lower-T phase)
  t1 <- subcrt(x,P=0,T=thermo$obigt$z.T[x],convert=FALSE,do.phases=FALSE)
  t2 <- subcrt(x+1,P=0,T=thermo$obigt$z.T[x],convert=FALSE,do.phases=FALSE)
  # if these aren't the same mineral all we can say is zero
  if(as.character(t1$species$name) != as.character(t2$species$name)) return(0)
  # we really hope the G's are the same ...
  #if(abs(t2$out[[1]]$G - t2$out[[1]]$G) > 0.1) warning('dP.dT: inconsistent values of G for different phases of ',x,call.=FALSE)
  dP.dT <- convert( ( t2$out[[1]]$S - t1$out[[1]]$S ) / ( t2$out[[1]]$V - t1$out[[1]]$V ), 'cm3bar' )
  return(dP.dT)
}

Ttr <- function(x,P=1,dPdT=NULL) {
  # calculate a phase transition temperature
  # taking account of P (from dP/dT)
  T <- thermo$obigt$z.T[x]
  if(is.null(dPdT)) dPdT <- dPdTtr(x)
  return(T + P / dPdT)
}

basis.comp <- function(basis) {
  #if(any(!is.numeric(basis))) stop('argument must be numeric index of species')
  #else if(any(basis < 0)) stop('numeric argument must be positive')
  # retrieve basis stoichiometries of species
  elements <- colnames(basis())
  comp0 <- expand.formula(rownames(basis()))
  for(i in 1:length(basis)) {
    if(identical(basis,0)) makeup <- makeup(thermo$obigt$formula[1]) 
    else if(is.numeric(basis[i])) makeup <- makeup(as.character(thermo$obigt$formula[basis[i]])) 
    else makeup <- makeup(basis[i])
    comprow <- expand.formula(elements,makeup)
    # solve for the basis coefficients
    comprow <- solve(t(basis()),comprow)
    comprow <- data.frame(matrix(comprow,byrow=TRUE,nrow=1))
    colnames(comprow) <- rownames(basis())
    missingelements <- which(!rownames(makeup) %in% elements)
    if(length(missingelements)>0) {
      if(length(missingelements)==1) p <- '' else p <- 's'
      cat(paste('basis.comp: missing element',p,' ',
        c2s(rownames(makeup)[missingelements],sep=' '),
        ' from species ',basis[i],'.\n',sep=''))
      rownames(comprow) <- 'missing'
    }
    # strip small numbers
    comprow <- makeup(comprow,FALSE)
    if(length(basis)==1) return(comprow)
    else if(i==1) comp <- rbind(comp0,comprow)
    else comp <- rbind(comp,comprow)
  }
  return(comp)
} 
    

ZC <- function(x) {
  # calculate nominal carbon oxidation state of chemical formulas
  # FIXME: add more elements, warn about missing ones
  # defactorize
  if(is.factor(x[1])) x <- as.character(x)
  # make sure all the elements we want are listed (even zeros)
  ZC <- numeric()
  for(i in 1:length(x)) {
    m <- as.data.frame(t(makeup(c(x[i],'C0H0N0O0S0Z0'))))
    ZC <- c(ZC,(-1*m$H+3*m$N+2*m$O+2*m$S+m$Z)/m$C)
  }
  return(ZC)
}

MP90.cp <- function(T,protein) {
  # T (temperature, degrees C), protein (name of protein)
  # returns heat capacity of protein (kj/mol)
  # using algorithm of makhatadze and privalov, 1990.
  TMP <- c(5,25,50,75,100,125)
  A.cp <- splinefun(TMP,c(175.7,166.7,156.2,144.7,134.6,124.1))
  C.cp <- splinefun(TMP,c(225.4,237.6,250.8,260.7,268.2,276.1))
  D.cp <- splinefun(TMP,c( 72.8, 89.0,106.2,124.5,140.7,154.3))
  E.cp <- splinefun(TMP,c(168.3,179.0,192.0,203.7,211.4,217.8))
  F.cp <- splinefun(TMP,c(395.7,383.0,370.3,358.4,348.3,339.6))
  G.cp <- splinefun(TMP,c( 82.3, 78.0, 71.7, 66.4, 59.7, 53.9))
  H.cp <- splinefun(TMP,c(205.7,179.6,177.2,179.6,187.1,196.8))
  I.cp <- splinefun(TMP,c(406.8,402.3,397.1,390.8,386.0,380.8))
  K.cp <- splinefun(TMP,c(328.8,332.5,334.0,337.5,339.4,343.6))
  L.cp <- splinefun(TMP,c(385.9,381.7,377.8,372.9,369.4,365.5))
  M.cp <- splinefun(TMP,c(197.1,175.9,158.1,150.3,148.1,143.9))
  N.cp <- splinefun(TMP,c( 72.9, 88.8,109.8,125.2,140.5,154.2))
  P.cp <- splinefun(TMP,c(214.6,177.7,152.3,142.8,135.6,130.1))
  Q.cp <- splinefun(TMP,c(168.0,180.2,193.4,203.3,210.8,218.7))
  R.cp <- splinefun(TMP,c(204.6,273.4,305.8,315.1,318.7,318.5))
  S.cp <- splinefun(TMP,c( 75.6, 81.2, 85.7, 91.4, 97.3,102.1))
  T.cp <- splinefun(TMP,c(194.2,184.5,182.2,186.5,199.0,216.2))
  V.cp <- splinefun(TMP,c(324.6,314.4,305.0,294.7,285.7,269.6))
  W.cp <- splinefun(TMP,c(471.2,458.5,445.8,433.9,423.8,415.1))
  Y.cp <- splinefun(TMP,c(310.6,301.7,295.2,294.5,300.1,304.0))
  AA.cp <- splinefun(TMP,c(-158.3,-90.4,-21.5,-32.3,-92.4,-150.0))
  UPBB.cp <- splinefun(TMP,c(3.7,15.2,26.2,29.8,33.7,33.7))
  cnew <- numeric()
  for(i in 1:length(T)) {
    Ti <- T[i]
    cp <- c(A.cp(Ti),C.cp(Ti),D.cp(Ti),E.cp(Ti),F.cp(Ti),
            G.cp(Ti),H.cp(Ti),I.cp(Ti),K.cp(Ti),L.cp(Ti),
            M.cp(Ti),N.cp(Ti),P.cp(Ti),Q.cp(Ti),R.cp(Ti),
            S.cp(Ti),T.cp(Ti),V.cp(Ti),W.cp(Ti),Y.cp(Ti))
    # get the protein composition
    tt <- protein(protein)[,6:25]
    cnew <- c(cnew, sum(cp * as.numeric(tt)) + sum(as.numeric(tt)) * UPBB.cp(Ti))
  }
  return(cnew)
}

describe <- function(x=NULL,T=NULL,P=NULL,use.name=FALSE,as.reaction=NULL,digits=1) {
  # summarize (write) a reaction or the thermodynamic potentials
  # 'x' can be a reaction dataframe (output by subcrt)
  # or the basis dataframe (thermo$basis)
  do.reaction <- FALSE
  if(!is.null(x)) {
    if('coeff' %in% colnames(x)) do.reaction <- TRUE
    if(do.reaction) {
      use.name <- rep(use.name,length.out=nrow(x))
      un.r <- un.p <- logical()
      # we should group the products and reactants
      x.reactants <- x.products <- x[0,]
      for(i in 1:nrow(x)) {
        if(x$coeff[i] < 0) {
          x.reactants <- rbind(x.reactants,x[i,])
          un.r <- c(un.r,use.name[i])
        }
        else if(x$coeff[i] > 0) {
          x.products <- rbind(x.products,x[i,])
          un.p <- c(un.p,use.name[i])
        }
      }
      use.name <- c(un.r,un.p)
      if(is.null(as.reaction)) {
        as.reaction <- FALSE
        # if there are reactants and products we say it's a reaction
        if(nrow(x.reactants) > 0 & nrow(x.products) > 0) as.reaction <- TRUE
        if(as.reaction) x <- rbind(x.reactants,x.products)
      }
    }
  }
  mystuff <- character()
  # this function takes T, P arguments in K, bar,
  # which are converted as necessary ...
  if(is.null(T)) if(missing(T)) if(!do.reaction) T <- thermo$opt$Tr
  if(is.null(P)) if(missing(P)) if(!do.reaction) P <- thermo$opt$Pr
  T.units <- nuts('T')
  if(T.units=='C') T <- convert(T,'C')
  if(!is.null(T)) mystuff <- paste(T,' ',T.units,', ',sep='')
  P.units <- nuts('P')
  if(P.units=='MPa') P <- convert(P,'MPa')
  if(!is.null(P)) mystuff <- paste(mystuff,P,' ',P.units,', ',sep='')
  # put together the basis text
  if(!is.null(x)) {
    for(i in 1:nrow(x)) {
      if(do.reaction) {
        coefftext <- abs(x$coeff[i])
        if(as.reaction) if(coefftext==1) coefftext <- ''
        if(!as.reaction) {
          #if(x$coeff[i]==1) coefftext <- '' 
          if(i==1) coefftext <- x$coeff[i]
        }
        if(use.name[i]) speciestext <- x$name[i] else speciestext <- x$formula[i]
        if(i==nrow(x)) optext <- '' else {
          optext <- ' + '
          if(sign(x$coeff[i])==-1 & sign(x$coeff[i+1])==1) optext <- ' = '
          if(!as.reaction) if(sign(x$coeff[i+1])==-1) optext <- ' - ' else optext <- ' + '
        } 
        mystuff <- paste(mystuff,coefftext,speciestext,optext,sep='')
      } else {
        if(x$state[i]=='gas') activity <- 'f' else activity <- 'a'
        if(can.be.numeric(x$logact[i])) logacttext <- round(as.numeric(x$logact[i]),digits)
        else logacttext <- x$logact[i]
        mystuff <- paste(mystuff,activity,rownames(x)[i],'(',logacttext,'), ',sep='')
      }
    }
  }
  # backspace over the last bit of punctuation
  if(!do.reaction) mystuff <- substr(mystuff,1,nchar(mystuff)-2)
  return(mystuff)
}

mod.buffer <- function(name,species=NULL,state=thermo$opt$state,logact=-3) {
  # 20071102 add or change a buffer system
  if(is.null(species)) {
    iname <- which(name==thermo$buffer$name)
    if(length(iname)>0) species <- thermo$buffer$species[iname]
    else species <- character()
  }
  ls <- length(species)
  if(ls < length(name) | ls < length(state) | ls < length(logact))
    stop('species must be at least as long as the other arguments')
  if(length(name)!=ls) name <- rep(name,length.out=ls)
  add <- TRUE
  if(TRUE %in% (name %in% thermo$buffer$name)) {
    add <- FALSE
    imod <- which(thermo$buffer$name %in% name & thermo$buffer$species %in% species)
    if(length(imod)>0) {
      if(state[1]=='') {
        thermo$buffer <<- thermo$buffer[-imod,]
        cat(paste('mod.buffer: removed ',c2s(species),' in ',
          c2s(unique(name)),' buffer.\n',sep=''))
      } else {
        if(missing(state)) state <- thermo$buffer$state[imod]
        if(missing(logact)) logact <- thermo$buffer$logact[imod]
        if(length(state)!=ls) state <- rep(state,length.out=ls)
        if(length(logact)!=ls) logact <- rep(logact,length.out=ls)
        state.old <- thermo$buffer$state[imod]
        logact.old <- thermo$buffer$logact[imod]
        thermo$buffer$state[imod] <<- state
        thermo$buffer$logact[imod] <<- logact
        if(identical(state.old,state) & identical(logact.old,logact)) {
          cat(paste('mod.buffer: nothing changed for ',
            c2s(species),' in ',c2s(unique(name)),' buffer.\n',sep=''))
        } else {
          cat(paste('mod.buffer: changed state and/or logact of ',
            c2s(species),' in ',c2s(unique(name)),' buffer.\n',sep=''))
        }
      }
    } else {
      add <- TRUE
    }
  } 
  if(add) {
    if(state[1]=='') state <- rep(thermo$opt$state,length.out=ls)
    t <- data.frame(name=name,species=species,state=state,logact=logact)
    thermo$buffer <<- rbind(thermo$buffer,t)
    cat(paste('mod.buffer: added ',c2s(unique(name)),'.\n',sep=''))
  }
  return(invisible(thermo$buffer[thermo$buffer$name %in% name,]))
}

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

examples <- function() {
  # run all the examples in CHNOSZ
  .ptime <- proc.time()
  example(CHNOSZ,ask=FALSE)
  example(thermo,ask=FALSE)
  example(examples,ask=FALSE)  # i.e., "utilities"
  example(info,ask=FALSE)
  example(hkf,ask=FALSE)  # i.e., "eos"
  example(water,ask=FALSE)
  example(subcrt,ask=FALSE)
  example(nuts,ask=FALSE)
  example(makeup,ask=FALSE)
  example(basis,ask=FALSE)
  example(species,ask=FALSE)
  example(affinity,ask=FALSE)
  example(diagram,ask=FALSE)
  example(buffer,ask=FALSE)
  example(protein,ask=FALSE)
  example(ionize,ask=FALSE)
  example(get.protein,ask=FALSE)
  example(diversity,ask=FALSE)
  example(transfer,ask=FALSE)
  cat("Time elapsed: ", proc.time() - .ptime, "\n")
}

CHNOSZ <- function() {
  # just runs the CHNOSZ.Rd examples
  example(CHNOSZ,ask=FALSE)
}

nonideal <- function(species,proptable,IS,T) {
  # generate nonideal contributions to thermodynamic properties
  # number of species, same length as proptable list
  # T in Kelvin, same length as nrows of proptables
  # a function that does a lot of the work
  loggamma2 <- function(Z,I,T,prop='log') {
    # extended Debye-Huckel equation ('log')
    # and its partial derivatives ('G','H','S','Cp')
    # T in Kelvin
    B <- 1.6
    # equation for A from Clarke and Glew, 1980
    #A <- expression(-16.39023 + 261.3371/T + 3.3689633*log(T)- 1.437167*(T/100) + 0.111995*(T/100)^2)
    # equation for alpha from Alberty, 2003 p. 48
    A <- alpha <- expression(1.10708 - 1.54508E-3 * T + 5.95584E-6 * T^2)
    # from examples for deriv to take first and higher-order derivatives
    DD <- function(expr,name, order = 1) {
      if(order < 1) stop("'order' must be >= 1")
      if(order == 1) D(expr,name)
      else DD(D(expr, name), name, order - 1)
    }
    # Alberty, 2003 Eq. 3.6-1
    loggamma <- function(a,Z,I,B) { - a * Z^2 * I^(1/2) / (1 + B * I^(1/2)) }
    # XXX check the following equations 20080208 jmd @@@
    if(prop=='log') return(loggamma(eval(A),Z,I,B))
    else if(prop=='G') return(thermo$opt$R * T * loggamma(eval(A),Z,I,B))
    else if(prop=='H') return(thermo$opt$R * T^2 * loggamma(eval(DD(A,'T',1)),Z,I,B))
    else if(prop=='S') return(- thermo$opt$R * T * loggamma(eval(DD(A,'T',1)),Z,I,B))
    else if(prop=='Cp') return(thermo$opt$R * T^2 *loggamma(eval(DD(A,'T',2)),Z,I,B))
  }
  if(!is.numeric(species[[1]])) species <- info(species,'aq')
  proptable <- as.list(proptable)
  # which gamma function to use
  #mlg <- get(paste('loggamma',which,sep=''))
  ndid <- 0
  for(i in 1:length(species)) {
    myprops <- proptable[[i]]
    Z <- as.numeric(makeup(species[i],'Z'))
    # don't do anything for neutral species
    if(Z==0) next
    # this would keep unit activity coefficient of the proton and electron
    #if(species[i] %in% c(info('H+',quiet=TRUE),info('e-',quiet=TRUE))) next
    didit <- FALSE
    for(j in 1:ncol(myprops)) {
      #if(colnames(myprops)[j]=='G') myprops[,j] <- myprops[,j] + thermo$opt$R * T * mlg(Z,IS,convert(T,'C'))
      cname <- colnames(myprops)[j]
      if(cname %in% c('G','H','S','Cp')) {
        myprops[,j] <- myprops[,j] + loggamma2(Z,IS,T,cname)
        didit <- TRUE
      }
    }
    # append a loggam column if we did any nonideal calculations of thermodynamic properties
    if(didit) myprops <- cbind(myprops,loggam=loggamma2(Z,IS,T))
    proptable[[i]] <- myprops
    if(didit) ndid <- ndid + 1
  }
  if(ndid > 0) cat(paste('nonideal:',ndid,'species were nonideal.\n'))
  return(proptable)
}

add.obigt <- function(file="obigt.csv") {
  # restore our entries in thermo$obigt
  # from values saved in a file
  to1 <- thermo$obigt
  to2 <- read.csv(file)
  j <- which(!to2$name %in% to1$name)
  if(length(j)>0) to1 <- rbind(to1,to2[j,])
  thermo$obigt <<- to1
  cat(paste("add.obigt: added",length(j),"species from",file,"\n"))
}

# which.balance is used by transfer(),
# keep it as a global function
which.balance <- function(species) {
  # find the first basis species that
  # is present in all species of interest
  # ... it can be used to balance the system
  nbasis <- function(species) return(ncol(species)-4)
  ib <- NULL
  nb <- 1
  nbs <- nbasis(species)
  for(i in 1:nbs) {
    coeff <- species[,i]
    if(length(coeff)==length(coeff[coeff!=0])) {
      ib <- c(ib,i)
      nb <- nb + 1
    } else ib <- c(ib,NA)
  }
  return(ib[!is.na(ib)])
}

grep.file <- function(file,x="",y=NULL,ignore.case=TRUE,startswith=">") {
  # return the line numbers of the
  # file that contain
  # the search term x
  # and optionally don't contain y
  lines <- readLines(file)
  ix <- grep(x,lines,ignore.case=ignore.case)
  if(!is.null(y)) {
    iy <- grep(y,lines,ignore.case=ignore.case)
    ix <- ix[!ix %in% iy] 
  }
  if(!is.null(startswith)) {
    ihead <- which(substr(lines,1,1)==startswith)
    ix <- ix[ix %in% ihead]
  }
  return(ix)
}

mylapply <- function(X,FUN,...) {
  # wrapper to run lapply or mclapply
  if(length(X) > 20) {
    if("multicore" %in% (.packages())) do.call("mclapply",list(X,FUN,...))
    else lapply(X,FUN,...)
  } else lapply(X,FUN,...)
}

.First.lib <- function(lib,pkg) {
  # version figuring adapted from package mgcv
  pkghelp <- library(help=CHNOSZ)$info[[1]]
  # things are different for older versions of R
  if(length(pkghelp)==1) pkghelp <- library(help=CHNOSZ)$info[[2]]
  version <- pkghelp[pmatch("Version",pkghelp)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  date <- pkghelp[pmatch("Date",pkghelp)]
  um <- strsplit(date," ")[[1]]
  date <- um[nchar(um)>0][2]
  # identify the program and version
  cat(paste('CHNOSZ version ',version,' (',date,')\n',sep=''))
  # load the data object
  data(thermo)
  # load data files in user's directory
  if(file.exists('protein.csv')) {
    tr <- try(rbind(read.csv('protein.csv',as.is=TRUE),thermo$protein),silent=TRUE)
    if(identical(class(tr),'try-error')) cat("thermo: protein.csv in current directory is not compatible with thermo$protein data table.\n")
    else {
      add.protein("protein.csv")
    }
  }
  if(file.exists('obigt.csv')) {
    tr <- try(rbind(read.csv('obigt.csv',as.is=TRUE),thermo$obigt),silent=TRUE)
    if(identical(class(tr),'try-error')) cat("thermo: obigt.csv in current directory is not compatible with thermo$obigt data table.\n")
    else {
      add.obigt("obigt.csv")
    }
  }
}


