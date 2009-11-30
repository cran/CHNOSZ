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
    stop('thermo.args: properties ',c2s(prop[notprop]),' not in ',c2s(props),'\n')
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

nucleicacids <- function(seq=NULL,type="DNA",comp=NULL,comp2=NULL) {
  # count bases or compute the formula, e.g.
  # n <- nucleicacids(list("AGCT","TTTT"))  # a dataframe of counts
  # f <- nucleicacids(n)  # a series of formulas
  # 20090926 jmd
  if(is.null(seq)) stop("please provide a sequence")
  if(type=="DNA") {
    na <- c("A","C","G","T")
    na.NA <- c("adenine","cytosine","guanine","thymine")
  } else if(type=="RNA") {
    na <- c("U","G","C","A")
    na.NA <- c("uracil","guanine","cytosine","adenine")
  } else stop(paste("invalid type:",type))
  if(is.data.frame(seq)) {
    # compute the chemical formula of bases
    if(!all(na %in% colnames(seq))) {
      nabases <- c2s(na[which(!na %in% colnames(seq))],sep=" ")
      stop(paste("requested type is",type,"but",nabases,"is/are not in the colnames of supplied dataframe"))
    }
    f.base <- thermo$obigt$formula[info(na.NA)]
    f0 <- makeup("C0H0N0O0")
    for(i in 1:4) {
      m.base <- makeup(f0,makeup(f.base[i]))
      if(i==1) m <- m.base else m <- cbind(m,m.base)
    }
    m <- as.matrix(m)
    f.out <- character()
    for(i in 1:nrow(seq)) {
      s <- t(seq[i,])
      c <- m %*% s
      c <- as.data.frame(c,row.names=row.names(c))
      colnames(c) <- "count"
      f <- makeup(makeup(c,""),"")
      f.out <- c(f.out,f)
    }
    return(f.out)
  } else {
    # count the numbers of nucleic acid bases in a sequence
    # sequences are given as elements of the list seq
    # to count the number of each amino acids in a sequence
    count.na <- function(na,seq) sum(seq==na)
    count.i <- function(i,seq) as.numeric(lapply(na,count.na,strsplit(toupper(seq[i]),"")[[1]]))
    # count bases in each sequence
    n <- t(as.data.frame(mylapply(1:length(seq),count.i,seq),optional=TRUE))
    n <- as.data.frame(n)
    # clean up row/column names
    colnames(n) <- na
    rownames(n) <- 1:nrow(n)
    # return the complement if requested e.g.
    # nucleicacids(x,type,"DNA")  # DNA complement
    # nucleicacids(x,type,"RNA")  # RNA complement
    # nucleicacids(x,type,"DNA","RNA")  # DNA, then RNA complement
    if(!is.null(comp)) {
      if(comp=="DNA") colnames(n) <- c("T","G","C","A")
      else if(comp=="RNA") colnames(n) <- c("U","G","C","A")
      else stop(paste("invalid complement request:",comp))
    }
    if(!is.null(comp2)) {
      if(comp2=="DNA") colnames(n) <- c("A","C","G","T")
      else if(comp2=="RNA") colnames(n) <- c("A","C","G","U")
      else stop(paste("invalid complement request:",comp))
    }
    return(n)
  }
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
    cat(paste('nuts: temperature in',thermo$opt$T.units,'\n'))
    cat(paste('nuts: energy in',thermo$opt$E.units,'\n'))
    cat(paste('nuts: pressure in',thermo$opt$P.units,'\n'))
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
      cat(paste('nuts: temperature in',thermo$opt$T.units,'\n'))
    }
    if(units[i] %in% c('j','cal')) {
      if(units[i]=='j') thermo$opt$E.units <<- 'J'
      if(units[i]=='cal') thermo$opt$E.units <<- 'cal'
      cat(paste('nuts: energy in',thermo$opt$E.units,'\n'))
    }
    if(units[i] %in% c('bar','mpa')) {
      if(units[i]=='bar') thermo$opt$P.units <<- 'bar'
      if(units[i]=='mpa') thermo$opt$P.units <<- 'MPa'
      cat(paste('nuts: pressure in',thermo$opt$P.units,'\n'))
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

thermo.plot.new <- function(xlim,ylim,xlab,ylab,cex=par('cex'),mar=NULL,lwd=par('lwd'),side=c(1,2,3,4),
  mgp=c(1.2,0.3,0),cex.axis=par('cex'),col=par('col'),yline=NULL,axs='i',do.box=TRUE,ticks=NULL) {
  # start a new plot with some customized settings
  # 20091108 changed argument name from 'ticks' to 'side' but
  # keep 'ticks' for backward compatibility
  if(!is.null(ticks)) side <- ticks 
  # 20090324 mar handling: NULL - a default setting; NA - par's setting
  # 20090413 changed mar of top side from 2 to 2.5
  if(is.null(mar)) mar <- c(3,3.5,2.5,1) else if(is.na(mar[1])) mar <- par('mar')
  par(mar=mar,mgp=mgp,tcl=0.3,las=1,xaxs=axs,yaxs=axs,cex=cex,lwd=lwd)
  plot.new()
  plot.window(xlim=xlim,ylim=ylim)
  if(do.box) box()
  # labels
  thermo.axis(xlab,side=1,line=mgp[1],cex=cex.axis,lwd=NULL)
  if(is.null(yline)) yline <- mgp[1]
  thermo.axis(ylab,side=2,line=yline,cex=cex.axis,lwd=NULL)
  # (optional) tick marks
  if(1 %in% side) thermo.axis(NULL,side=1,lwd=lwd,col=par('col'))
  if(2 %in% side) thermo.axis(NULL,side=2,lwd=lwd,col=par('col'))
  if(3 %in% side) thermo.axis(NULL,side=3,lwd=lwd,col=par('col'))
  if(4 %in% side) thermo.axis(NULL,side=4,lwd=lwd,col=par('col'))
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
      logK <- subcrt(c('H2O','oxygen','hydrogen'),c(-1,0.5,1),T=T,P=P,convert=FALSE)$out$logK 
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
      logK <- subcrt(c('H2O','oxygen','hydrogen'),c(-1,0.5,1),T=T,P=P,convert=FALSE)$out$logK 
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
  # 20090826: just return the argument if a comma is already present
  if(length(grep(",",x)) > 0) return(x)
  lab <- x
  if(missing(opt)) do.opt <- TRUE else do.opt <- FALSE
  if(!is.null(opt)) if(is.na(opt)) do.opt <- TRUE
  if(lab %in% c('T','P','Eh','pH','pe','logK','IS')) {
    # the label is one of these properties
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
    # the label is a chemical activity or fugacity
    if(is.null(thermo$basis)) rn <- '' else rn <- rownames(basis())
    if(lab %in% rn) {
      # 20090215: the state this basis species is in
      if(do.opt) opt <- as.character(thermo$basis$state)[match(lab,rn)]
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
              if(lc==-1) lc <- "-"
              else if(lc==1) lc <- "+"
              else if(lc > 0) lc <- paste("+",as.character(lc),sep="")
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
      # e.g. axis.label('DG0r','k') for standard molal Gibbs energy
      # change of reaction in kcal/mol
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

examples <- function(do.png=FALSE,long=TRUE) {
  # run all the examples in CHNOSZ
  .ptime <- proc.time()
  # "examples" refers to "utilities" help file
  # "hkf" refers to "eos"
  topics <- c("CHNOSZ","thermo","examples","info","hkf","water","subcrt",
    "nuts","makeup","basis","species","affinity","diagram","buffer",
    "protein","ionize","get.protein","revisit","transfer")
  do.plot <- FALSE
  if(is.character(do.png))
    png(paste(do.png,"%d.png",sep=""),width=700,height=700,pointsize=18)
  else if(do.png) do.plot <- TRUE
  for(i in 1:length(topics)) {
    if(do.plot) png(paste(topics[i],"%d.png",sep=""),width=700,height=700,pointsize=18)
    myargs <- list(topic=topics[i],ask=FALSE)
    do.call(example,myargs)
    if(long) longex(topics[i])
    if(do.plot) dev.off()
  }
  if(is.character(do.png)) dev.off()
  cat("Time elapsed: ", proc.time() - .ptime, "\n")
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

grep.file <- function(file,x="",y=NULL,ignore.case=TRUE,startswith=">",lines=NULL,grep="grep") {
  # return the line numbers of the file that contain
  # the search term x and optionally don't contain y
  # use the system's grep if available and y is NULL
  # TODO: better test for system (include eg MacOSX etc)
  if(is.null(y) & Sys.info()[[1]]=="Linux" & is.null(lines)) {
    if(is.null(startswith)) startswith <- "" else startswith <- paste("^",startswith,".*",sep="")
    if(ignore.case) ic <- "-i" else ic <- ""
    # 20091021 changed grep to egrep
    sysexp <- paste(grep,'-n',ic,paste('"',startswith,x,'"',sep=""),file,' | cut -f 1 -d ":"')
    ix <- as.integer(system(sysexp,intern=TRUE))
  } else {
    if(is.null(lines)) lines <- readLines(file)
    ix <- grep(x,lines,ignore.case=ignore.case)
    if(!is.null(y)) {
      iy <- grep(y,lines,ignore.case=ignore.case)
      ix <- ix[!ix %in% iy] 
    }
    if(!is.null(startswith)) {
      ihead <- which(substr(lines,1,1)==startswith)
      ix <- ix[ix %in% ihead]
    }
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

spearman <- function(a,b) {
  # calculate Spearman's rho (rank correlation coefficient)
  # based on help(dSpearman) in package SuppDists
  if(length(a)!=length(b)) stop("a and b must be same length")
  ra <- rank(a)
  rb <- rank(b)
  dr <- rb - ra
  d <- sum(dr^2)
  r <- length(a)
  x <- 1-6*d/(r*(r^2-1))
  return(x)
}

unitize <- function(logact=NULL,length=NULL,logact.tot=0) {
  # scale the logarithms of activities given in loga
  # so that the logarithm of total activity of residues
  # is zero (i.e. total activity of residues is one),
  # or some other value set in loga.tot.
  # length indicates the number of residues in each species.
  # if loga is NULL, take the logarithms of activities from
  # the current species definition. if any of those species
  # are proteins, get their lengths using protein.length.
  if(is.null(logact)) {
    if(is.null(thermo$species)) stop("loga is NULL and no species are defined")
    ts <- species()
    logact <- ts$logact
    length <- rep(1,length(logact))
    ip <- grep("_",ts$name)
    if(length(ip) > 0) length[ip] <- protein.length(ts$name[ip])
  }
  # the lengths of the species
  if(is.null(length)) length <- 1
  length <- rep(length,length.out=length(logact)) 
  # remove the logarithms
  act <- 10^logact
  # the total activity
  act.tot <- sum(act*length)
  # the target activity
  act.to.get <- 10^logact.tot
  # the factor to apply
  act.fact <- act.to.get/act.tot
  # apply the factor
  act <- act * act.fact
  # take the logarithms
  logact <- log10(act)
  # done!
}

mtitle <- function(main,...) {
  # make a possibly multi-line plot title 
  # useful for including expressions on multiple lines 
  l <- length(main)
  for(i in 1:l) mtext(main[i],line=l-i,...)
}

longex <- function(which) {
  if(which=="transfer") {
    ## react potassium feldspar in a closed system
    ## after Steinmann et al., 1994 and Helgeson et al., 1969
    basis(c("Al+3","H4SiO4","K+","H2O","H+","O2"),c(0,-6,-6,0,0,0))
    species(c("k-feldspar","muscovite","pyrophyllite","kaolinite","gibbsite"))
    a <- affinity(H4SiO4=c(-6,-2),"K+"=c(-3,8))
    diagram(a)
    basis("pH",4)
    species(1:5,c(-4,rep(-999,4)))
    t <- transfer(550,dmode="coupled",plot=c(2,3),devmax=0.2)
    # plot the output from transfer
    draw.transfer(t)
    # reset the plot layout
    layout(matrix(1))
    # note, the same calculation can be performed using
    # feldspar("closed")

    ## react APC2 from Saccharomyces cerevisiae
    ## to other proteins in the anaphase-promoting complex
    basis(c("CO2","H2O","NH3","H2","H2S"),c(-10,0,-4,-10,-7))
    species(c("APC1","APC2","APC5","CDC16","APC10","SWM1"),"YEAST")
    a <- affinity(CO2=c(-10,0),H2=c(-10,0))
    diagram(a,as.residue=TRUE)
    species(1:nrow(species()),-999)
    species("APC2_YEAST",0)
    t <- transfer(510,ibalance="PBB",plot=c(1,4),dmode="coupled",devmax=0.15)
    # this calculation is also available with
    # apc("closed")
  } else if(which=="get.protein") {
    ## Oxygen fugacity - activity of H2O predominance 
    ## diagrams for proteologs for 23 YeastGFP localizations
    # arranged by decreasing metastability:
    # order of this list of locations is based on the 
    # (dis)appearance of species on the current set of diagrams
    names <- c("vacuole","early.Golgi","ER","lipid.particle",
      "cell.periphery","ambiguous","Golgi","mitochondrion",
      "bud","actin","cytoplasm","late.Golgi",
      "endosome","nucleus","vacuolar.membrane","punctate.composite",
      "peroxisome","ER.to.Golgi","nucleolus","spindle.pole",
      "nuclear.periphery","bud.neck","microtubule")
    nloc <- c(4,5,3,4,4,3)
    inames <- 1:length(names)
    # define the system
    basis("CHNOS+")
    # calculate amino acid compositions using "get.protein" function 
    for(i in 1:length(names)) {
      y <- yeastgfp(names[i])
      p <- get.protein(y$yORF,"SGD",y$abundance,names[i])
      add.protein(p)
    }
    species(names,"SGD")
    t <- affinity(H2O=c(-5,0,256),O2=c(-80,-66,256))
    # setup the plot
    layout(matrix(c(1,1,2:7),byrow=TRUE,nrow=4),heights=c(0.7,3,3,3))
    par(mar=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,paste("Subcellular proteologs of S. cerevisiae,",
     "after Dick, 2009\n",describe(thermo$basis[-c(2,5),])),cex=1.5)
    opar <- par(mar=c(3,4,1,1),xpd=TRUE)
    for(i in 1:length(nloc)) {
      diagram(t,balance="PBB",names=names[inames],
        ispecies=inames,cex.axis=0.75)
      label.plot(letters[i])
      title(main=paste(length(inames),"locations"))
      # take out the stable species
      inames <- inames[-(1:nloc[i])]
    }
    # return to plot defaults
    layout(matrix(1))
    par(opar)

    ## Compare calculated and experimenal relative abundances
    ## of proteins in a subcellular location, after Dick, 2009
    # get the amino acid composition of the proteins
    loc <- "vacuolar.membrane"
    y <- yeastgfp(loc)
    ina <- which(is.na(y$abundance))
    p <- get.protein(y$yORF[-ina],"SGD")
    add.protein(p)
    # set up the system
    basis("CHNOS+")
    # this is the logfO2 value that gives the best fit (see paper)
    basis("O2",-74)
    is <- species(p$protein,p$organism)
    np <- length(is)
    pl <- protein.length(species()$name)
    # we use unitize so total activity of residues is unity
    loga <- rep(0,np)
    species(1:np,unitize(loga,pl))
    a <- affinity()
    d <- diagram(a,do.plot=FALSE)
    calc.loga <- as.numeric(d$logact)
    expt.loga <- unitize(log10(y$abundance[-ina]),pl)
    # which ones are outliers
    rmsd <- sqrt(sum((expt.loga-calc.loga)^2)/np)
    residuals <- abs(expt.loga - calc.loga)
    iout <- which(residuals > rmsd)
    pch <- rep(16,length(is))
    pch[iout] <- 1
    # the colors reflect average oxidation number of carbon
    ZC <- ZC(thermo$obigt$formula[species()$ispecies])
    col <- rgb(0.15-ZC,0,0.35+ZC,max=0.5)
    # there is a color-plotting error on line 567 of the plot.R file 
    # of Dick, 2009 that can be reproduced with
    #col <- rep(col,length.out=9)
    xlim <- ylim <- extendrange(c(calc.loga,expt.loga))
    thermo.plot.new(xlim=xlim,ylim=ylim,xlab=expression(list("log"*italic(a),
      "calc")),ylab=expression(list("log"*italic(a),"expt")))
    points(calc.loga,expt.loga,pch=pch,col=col)
    lines(xlim,ylim+rmsd,lty=2)
    lines(xlim,ylim-rmsd,lty=2)
    title(main=paste("Calculated and experimental relative abundances of\n",
      "proteins in ",loc,", after Dick, 2009",sep=""),cex.main=0.95)
    

    ### examples for stress response experiments

    ## predominance fields for overall protein 
    ## compositions induced by
    ## carbon, sulfur and nitrogen limitation
    ## (experimental data from Boer et al., 2003)
    expt <- c("low.C","low.N","low.S")
    for(i in 1:length(expt)) {
      p <- get.protein(expt[i],"SGD",abundance=1)
      add.protein(p)
    }
    basis("CHNOS+") 
    basis("O2",-75.29)
    species(expt,"SGD")
    a <- affinity(CO2=c(-5,0),H2S=c(-10,0))
    diagram(a,balance="PBB",names=expt,color=NULL)
    title(main=paste("Metastabilities of proteins induced by",
      "carbon, sulfur and nitrogen limitation",sep="\n"))

    ## predominance fields for overall protein compositions 
    ## induced and repressed in an/aerobic carbon limitation
    ## (experiments of Tai et al., 2005)
    # the activities of glucose, ammonium and sulfate
    # are similar to the non-growth-limiting concentrations
    # used by Boer et al., 2003
    basis(c("glucose","H2O","NH4+","hydrogen","SO4-2","H+"),
      c(-1,0,-1.3,999,-1.4,-7))
    # the names of the experiments in thermo$stress
    expt <- c("Clim.aerobic.down","Clim.aerobic.up",
      "Clim.anaerobic.down","Clim.anaerobic.up")
    # here we use abundance to indicate that the protein
    # compositions should be summed together in equal amounts
    for(i in 1:length(expt)) {
      p <- get.protein(expt[i],"SGD",abundance=1)
      add.protein(p)
    }
    species(expt,"SGD")
    a <- affinity(C6H12O6=c(-35,-20),H2=c(-20,0))
    diagram(a,color=NULL,as.residue=TRUE)
    title(main=paste("Metastabilities of average protein residues in",
      "an/aerobic carbon limitation in yeast",sep="\n"))
  } else if(which=="diagram") {
    ## Degrees of formation of ATP species as a function of 
    ## temperature, after LaRowe and Helgeson, 2007
    # Mg+2 here is the immobile component
    # cf. LH07, where bulk composition of Mg+2 is specified
    basis(c("CO2","NH3","H2O","H3PO4","O2","H+","Mg+2"),
      c(999,999,999,999,999,-5,-4))
    species(c("HATP-3","H2ATP-2","MgATP-2","MgHATP-"))
    diagram(affinity(T=c(0,120,25)),alpha=TRUE)  
    title(main=paste("Degrees of formation of ATP species,\n",
      "pH=5, log(aMg+2)=-3. After LaRowe and Helgeson, 2007"),
      cex.main=0.9)

    ## speciation of phosphate as a function of ionic strength
    basis("CHNOPS+")
    T <- c(25,100)
    species(c("HPO4-2","H2PO4-"))
    diagram(affinity(IS=c(0,0.25),T=T[1]),ylim=c(-3.2,-2.8))
    title(main=paste("Aqueous phosphate speciation, pH=7",
      "25 and 100 deg C - black and red lines",sep="\n"))
    diagram(affinity(IS=c(0,0.25),T=T[2]),ylim=c(-3.2,-2.8),add=TRUE,color="red")  
    ## phosphate predominance f(IS,pH)
    diagram(affinity(pH=c(6.8,7.2),IS=c(0,0.25),T=T[1]))
    title(main=paste("Aqueous phosphate predominance, 25 deg C",
      "and 100 deg C (dotted overlay)",sep="\n"))
    diagram(affinity(pH=c(6.8,7.2),IS=c(0,0.25),T=T[2]),add=TRUE,dotted=2,
      names=NULL)

    ## using strip(): equilibrium abundances of 
    ## proteins from different mammals
    organisms <- c("BOVIN","CANFA","HUMAN","MOUSE","PIG")
    proteins <- c("AQP1","MYG","P53")
    basis("CHNOS+")
    species(rep(proteins,each=5),organisms)
    a <- affinity(O2=c(-85,-65,250))
    ispecies <- list(1:5,6:10,11:15)
    desc <- c("(Aquaporin-1)","(Myoglobin)",
      "(Cellular tumor antigen p53)")
    names(ispecies) <- paste(proteins,desc)
    col <- rainbow(5)
    strip(a,ispecies=ispecies,ymin=-0.7,col=col)
    legend("bottomright",legend=organisms,col=col,
      lty=1,lwd=4,bty="n")
    title(main=paste("Equilibrium abundances of proteins",
      "from different mammals",sep="\n"))
    
    
    ### predominance diagrams (equal activities of
    ### species as a function of two variables) 

    ## T-P stability diagram for a unary system (SiO2)
    ## after Ernst, 1976, Fig. 4.4
    # the system is SiO2 (one component) but
    # the basis species require all the elements:
    # note that basis(c("SiO2","O2")) would identify SiO2(aq) 
    # which is okay but brings in calculations of properties 
    # of  water which are relatively slow.
    basis(c("quartz","O2"))  # cr, gas
    # browse the species of interest
    info("SiO2 ")
    # we'll take the crystalline ones
    t <- info("SiO2","cr")
    species(t)
    # calculate chemical affinities
    # the do.phases argument is necessary for 
    # dealing with the phase transitions of minerals
    t <- affinity(T=c(0,2000),P=c(0,30000),do.phases=TRUE)
    diagram(t)
    title(main="Phases of SiO2\nafter Ernst, 1976")

    ## Distribution of copper on Eh-pH diagram
    ## after Fig. 15 of Pourbaix, 1949
    basis(c("Cu+2","Cl-","H2O","H+","e-"))
    basis("Cl-",-2)
    # two aqueous species, three solid ones
    # (our CuCl is aq but it was cr in P's Fig. 15)
    species(c("CuCl","Cu+2","copper","cuprite","tenorite"))
    # (which is equivalent to ...)
    # species(c("CuCl","Cu+2","Cu","Cu2O","CuO"),c("","","","","","cr"))
    t <- affinity(pH=c(0,7),Eh=c(-0.1,1))
    # this one corresponds to activity contours of 
    # aqueous species at 10^-3 (the default aq activity in CHNOSZ)
    diagram(t,color=NULL)
    # here we set activities to unity; aq-cr boundaries change
    species(c("CuCl","Cu+2"),c(0,0))
    t <- affinity(pH=c(0,7),Eh=c(-0.1,1))
    diagram(t,add=TRUE,names=NULL,col="blue",color=NULL)
    water.lines()
    title(main=paste("H2O-Cl-(Cu); activities of 10^-3 (black)\n",
      "or 0 (blue); after Pourbaix, 1949"))

    ## Eh-pH diagrams for copper-water-glycine
    ## After Fig. 2 of Aksu and Doyle, 2001
    # update rows of the database
    i <- info(c("Cu(Gly)+","glycinate","e-","H+"))
    n <- nrow(thermo$obigt <<- rbind(thermo$obigt,thermo$obigt[rep(i[1],2),]))
    thermo$obigt$name[n-1] <<- "Cu(Gly)2-"
    thermo$obigt$name[n] <<- "HCu(Gly)+2"
    thermo$obigt$formula[n-1] <<- makeup(makeup(c(i[1],i[2],i[3]),""),"")
    thermo$obigt$formula[n] <<- makeup(makeup(c(i[1],i[4]),""),"")
    # In Fig 2b, total log activities of Cu (Cu_T) 
    # and glycine (L_T) are -4 and -1
    basis(c("Cu+2","H2O","H+","e-","glycine","CO2"),c(999,0,999,999,-1,999))
    # solid species
    species(c("copper","cuprite","tenorite"))
    # aqueous species
    species(c("glycinium","glycine","glycinate","Cu+","Cu+2","CuO2-2","HCuO2-",
      "Cu(Gly)+","Cu(Gly)2","Cu(Gly)2-","HCu(Gly)+2"),-4)
    ispecies <- species()$ispecies
    # update the Gibbs energies using A&D's Table 1 and Table II
    logK <- c(convert(convert(c(0,-146,-129.7,-384.061,-370.647,-314.833,
      49.98,65.49,-183.6,-258.5,-298.2)*1000,"cal"),"logK"),15.64,10.1,2.92) 
    # do it in order so later species take account of prev. species' values
    for(i in 1:length(logK)) {
      G <- convert(logK[i],"G")
      if(i==12) G <- G + thermo$obigt$G[ispecies[8]] + 
        2*thermo$obigt$G[ispecies[6]]
      if(i==13) G <- G + thermo$obigt$G[ispecies[7]] + 
        2*thermo$obigt$G[ispecies[6]]
      if(i==14) G <- G + thermo$obigt$G[ispecies[11]]
      thermo$obigt$G[ispecies[i]] <- G
    }  # done with changing Gibbs free energies!
    # we have to get some leftovers out of there or diagram() gets confused
    species(c("glycinium","glycine","glycinate"),delete=TRUE)
    # make a plot to see if it's working
    ispecies <- ispecies[-(1:6)]
    afun <- function(cu,gly) {
      # from above: our fifth basis species is glycine(-ate,-ium)
      basis(rownames(basis())[5],gly)
      t <- match(ispecies,species()$ispecies)
      species(t,cu)
      affinity(pH=c(0,16),Eh=c(-0.6,1.0))
    }
    diagram(afun(-4,-1))
    title(main=paste("Aqueous Copper + Glycine, 25 deg C, 1 bar",
      "After Aksu and Doyle, 2001 Fig. 2b",sep="\n"))
    # What's missing? Try glycinate not glycine in reactions at ph > ~9.8
    basis(c("Cu+2","H2O","H+","e-","glycinate","CO2"),
      c(999,0,999,999,-2,999))
    species(c("copper","cuprite","tenorite","Cu+","Cu+2","CuO2-2","HCuO2-",
      "Cu(Gly)+","Cu(Gly)2","Cu(Gly)2-","HCu(Gly)+2"))
    loga_Cu <- -4
    loga_Gly <- -1
    diagram(afun(loga_Cu,loga_Gly),color=NULL,col="blue",
      names=species()$name,col.names="blue",add=TRUE)
    water.lines()
    # the glycine ionization constants could be calculated using
    # subcrt, here they are taken from A&D Table II
    abline(v=c(2.35,9.778),lty=3)
    # now put glycinium (low pH) in the basis
    basis(c("Cu+2","H2O","H+","e-","glycinium","CO2"),c(999,0,999,999,-2,999))
    species(c("copper","cuprite","tenorite","Cu+","Cu+2","CuO2-2","HCuO2-",
      "Cu(Gly)+","Cu(Gly)2","Cu(Gly)2-","HCu(Gly)+2"))
    diagram(afun(loga_Cu,loga_Gly),color=NULL,col="green",
      names=NULL,col.names="green",add=TRUE)
    # let's forget our changes to 'thermo' so the example 
    # below that uses glycine will work as expected
    data(thermo)

    ## a pe-pH diagram for hydrated iron sulfides,
    ## goethite and pyrite, After Majzlan et al., 2006
    basis(c("Fe+2","SO4-2","H2O","H+","e-"),c(0,log10(3),log10(0.75),999,999))
    species(c("rhomboclase","ferricopiapite","hydronium jarosite",
      "goethite","melanterite","pyrite"))
    a <- affinity(pH=c(-1,4),pe=c(-5,23))
    diagram(a)
    water.lines(yaxis="pe")
    title(main=paste("Fe-S-O-H After Majzlan et al., 2006",
      describe(thermo$basis[2:3,],digits=2),sep="\n"))

    ## cysteine-cysteinate-cystine Eh-pH
    ## at 25 and 100 deg C
    basis("CHNOSe")
    species(c("cysteine","cysteinate","cystine"))
    a <- affinity(pH=c(5,10),Eh=c(-0.5,0))
    diagram(a,color=NULL)
    water.lines()
    a <- affinity(pH=c(5,10),Eh=c(-0.5,0),T=100)
    diagram(a,col="red",add=TRUE,names=NULL)
    water.lines(T=convert(100,"K"),col="red")
    title(main=paste("Cysteine Cysteinate Cystine",
      "25 and 100 deg C",sep="\n"))
    
    ## Soil Organic Matter CO2-H2O, O2-H2O (after Tardy et al., 1997)
    # NH3 is conserved, and H2O is on an axis of the 
    # diagrams, so their activities are nonsensical here.
    # formulas for aqueous species, names for phases ...
    basis(c("NH3","water","CO2","O2"),c(999,999,-2.5,-28))
    # switch to gaseous CO2 (aq is the default)
    basis("CO2","gas")
    # load the species of interest
    species(c("microflore","plantes","acide fulvique",
      "acide humique","humine"))
    # proceed with the diagrams
    diagram(affinity(H2O=c(-0.6,0.1),CO2=c(-3,-1)))
    title(main=paste("Relative stabilities of soil organic matter\n",
      "after Tardy et al., 1997"))
    # this would be the O2-H2O diagram
    # diagram(affinity(H2O=c(-1,0.5),O2=c(-28.5,-27.5)))

    ## Aqueous Aluminum Species F-/OH- (afer Tagirov and Schott, 2001)
    # in order to reproduce this calculation, we have to 
    # consider the properties of species given by these authors,
    # which are not the default ones in the database
    # (see last example of 'thermo' help page)
    thermo$opt$level <- 2
    # The coefficients on H+ and O2 in all the formation reactions
    # are zero, so the number of basis species here is three. Al+3 
    # becomes the conservant, and F- and OH- are being plotted ...
    # so their initial activities don't have to be set.
    basis(c("Al+3","F-","OH-","H+","O2"),rep(999,5))
    species(c("Al+3","Al(OH)4-","AlOH+2","Al(OH)2+","Al(OH)3",
      "AlF+2","AlF2+","AlF3","AlF4-","Al(OH)2F2-","Al(OH)2F","AlOHF2"))
    # Increase the x- and y- resolution from the default and calculate
    # and plot predominance limits. Names of charged basis species,
    # such as "H+", "e-" and the ones shown here, should be quoted
    # when given as arguments to affinity(). The OH- values shown here
    # correspond to pH=c(0,14) (at unit activity of water).
    a <- affinity("OH-"=c(-14,0),"F-"=c(-1,-8),T=200)
    diagram(a)
    title(main=paste("Aqueous aluminium species, T=200 C, P=Psat\n",
      "after Tagirov and Schott, 2001"))
    # We could do this to overlay boundaries for a different pressure
    #a.P <- affinity("OH-"=c(-14,0),"F-"=c(-1,-8),T=200,P=5000)
    #diagram(a.P,names=NULL,color=NULL,col="blue",add=TRUE)
    # restore data selection level to default
    thermo$opt$level <- 1

    ## Mineral equilibrium activity diagram
    ## (After Bowers et al., 1984, p. 246)
    # Consider the system denoted by BJH84 as 
    # HCl-H2O-CaO-CO2-MgO-(SiO2) at 300 deg C and 1000 bar
    # HCl, CO2, O2 have zero stoichiometric coeffs in the species
    # Ca+2, Mg+2 are going to be on the diagram
    # SiO2 will be conserved (the immobile component)
    # (try changing any of the 999's, the diagram will be unaffected)
    basis(c("HCl","H2O","Ca+2","CO2","Mg+2","SiO2","O2","H+"),
          c(999,0,999,999,999,999,999,-7))
    # we have to tell CHNOSZ which species to include; there are
    # others that could be in this system
    species(c("quartz","talc","forsterite","tremolite","diopside",
      "wollastonite","monticellite","merwinite"))
    # calculate the chemical affinities of formation reactions
    #t <- affinity("Mg+2"=c(-14,-8),"Ca+2"=c(-12,-4),T=300,P=1000)
    t <- affinity("Mg+2"=c(-12,-4),"Ca+2"=c(-8,0),T=300,P=1000)
    diagram(t)
    title(main=paste("HCl-H2O-CaO-CO2-MgO-(SiO2) at 300 deg C,\n",
      "1000 bar and pH=7. After Bowers et al., 1984"))
    # note: BJH84 use a different method for representing 
    # the axes of the diagrams, similar to (a_Ca+2)/(a_H+)^2,
    # so this is so far an approximate reproduction of their diagram.

    ## Fe-S-O at 200 deg C, After Helgeson, 1970
    basis(c("Fe","O2","S2"))
    species(c("iron","ferrous-oxide","magnetite",
      "hematite","pyrite","pyrrhotite"))
    a <- affinity(S2=c(-50,0),O2=c(-90,-10),T=200)
    diagram(a,color="heat")
    title(main=paste("Fe-S-O, 200 degrees C, 1 bar",
      "After Helgeson, 1970",sep="\n"))

    ## Nucleobase - Amino Acid Interaction Eh-H2O
    # for this example we try a unique basis definition
    basis(c("CO2","H2O","glutamine","e-","H+"),c(-3,0,-3,0,-7))
    species(c("uracil","cytosine","adenine","guanine",
      "phenylalanine","proline","lysine","glycine"),"aq")
    # this loaded four nucleobases and four related amino acids
    # (coded for by the homocodon triplets)
    # check out the predominance diagrams
    a.1 <- affinity(H2O=c(-5,0),Eh=c(-0.5,0))
    diagram(a.1,color=NULL)
    # overlay a different temperature
    a.2 <- affinity(H2O=c(-5,0),Eh=c(-0.5,0),T=100)
    diagram(a.2,col="red",add=TRUE,names=NULL)
    # start make a title for the plot
    tb <- thermo$basis   # includes activities of basis species
    # exclude those that are on the axes
    tb <- tb[!((rownames(tb) %in% c("e-","H2O"))),]
    title(main=paste("Nucleobases and amino acids,",
      "T=25 and 100 C at Psat\n",
      describe(tb,T=NULL,P=NULL)),cex.main=0.9)

    ## Temperature-Pressure: kayanite-sillimanite-andalusite
    basis(c("kyanite","quartz","oxygen"))
    species(c("kyanite","sillimanite","andalusite"))
    diagram(affinity(T=c(0,1000,200),P=c(500,5000,200)),color=NULL)
    title(main="Al2SiO5")
  } else if(which=="ionize") {
    ## Eh-pH diagrams for intra/extracellular proteins
    organism <- c("PYRFU","ECOLI","YEAST")
    intracellular <- c("AMPM","AMPM","AMPM1")
    extracellular <- c("O08452","AMY1","PST1")
    basis("CHNOSe")  # for Eh we need electrons
    mycol <- c("red","green","blue")
    for(i in 1:3) {
      species(delete=TRUE)
      species(c(intracellular[i],extracellular[i]),organism[i])
      if(i == 1) add <- FALSE else add <- TRUE
      t <- affinity(pH=c(0,14),Eh=c(-1,0))
      diagram(t,add=add,color=NULL,names=species()$name,
        col=mycol[i],col.names=mycol[i])
    }
    title(main=paste("Intracellular (AMPM) and extracellular proteins\n",
      describe(thermo$basis[1:4,])))
    
    ## Buffer + ionization: Metastablilities of
    ## thiol peroxidases from model bactera
    ## (ECOLI, BACSU mesophile; AQUAE thermophile,
    ## THIDA acidophile, BACHD alkaliphile)
    basis("CHNOS+")
    organisms <- c("ECOLI","AQUAE","BACSU","BACHD","THIDA")
    species("TPX",organisms)
    # create a buffer with our proteins in it
    mod.buffer("TPX",paste("TPX",organisms,sep="_"))
    # set up the buffered activities
    basis(c("CO2","H2O","NH3","O2"),"TPX")
    t <- affinity(return.buffer=TRUE,T=50)
    basis(c("CO2","H2O","NH3","O2"),as.numeric(t[1:4]))
    t <- affinity(pH=c(0,14,200),T=c(25,80,200))
    diagram(t,color=NULL)
    title(main=paste("Thiol peroxidases from bacteria\n",
      describe(thermo$basis[-6,],T=NULL)),cex.main=0.9)

    ## Buffer + ionization: Metastable assemblage
    ## for E. coli sigma factors on a T-pH diagram
    # (sigma factors 24, 32, 38, 54, 70, i.e.
    # RpoE, RpoH, RpoS, RpoN, RpoD)
    proteins <- c("RPOE","RP32","RPOS","RP54","RPOD")
    basis("CHNOS+")
    basis("pH",7.4)
    # define and set the buffer
    change("_sigma",paste(proteins,"ECOLI",sep="_"))
    basis(c("CO2","NH3","H2S","O2"),"sigma")
    logact <- affinity(return.buffer=TRUE,T=25)
    # Set the activities of the basis species to constants 
    # corresponding to the buffer, and diagram the relative
    # stabilities as a function of T and pH
    basis(c("CO2","NH3","H2S","O2"),as.numeric(logact))
    species(paste(proteins,"ECOLI",sep="_"))
    a <- affinity(pH=c(5,10),T=c(10,40))
    diagram(a,balance="PBB",residue=FALSE)
    title(main=paste("Sigma factors in E. coli\n",
      describe(thermo$basis[-6,],T=NULL)),cex.main=0.95)
  } else if(which=="revisit") {
    ## carboxylases from different organisms
    # ribulose phosphate carboxylase/oxygenase
    rubisco <- c("Q6JAI0_9RHOD","A5CKC7_9CHRO","RBL_SYNJA","A1E8R4_9CHLO",
      "A8C9T6_9MYCO","A3EQE1_9BACT","A6YF84_9PROT", "RBL_BRAJA",
      "A1RZJ5_THEPD","RBL_METJA","A3DND9_STAMF","RBL_PYRHO") 
    # acetyl-coenzyme A carboxylase
    accoaco <- c("Q05KD0_HYDTH","Q9F7M8_PRB01","A6CDM2_9PLAN","A0GZU2_9CHLR",
      "ACCA_DEIRA","A1VC70_DESVV","A7WGI1_9AQUI","Q2JSS7_SYNJA",
      "A4AGS7_9ACTN","ACCA_AQUAE","ACCA_CAUCR","A6VIX9_METM7")  
    # calculate affinities as a function of T with
    # buffered logfH2 and CO2.. adjust the glutathione buffer 
    # (more oxidizing than default)
    mod.buffer("GSH-GSSG",c("GSH","GSSG"),logact=c(-3,-7))
    # add a CO2 gas saturation buffer
    mod.buffer("CO2","CO2","gas",-1)
    basis(c("CO2","H2O","NH4+","SO4-2","H2","H+"),
      c("CO2",0,-4,-4,"GSH-GSSG",-7))
    species(c(rubisco,accoaco))
    a <- affinity(T=c(0,160))
    # create speciation diagram with colors
    col <- c(rep("red",12),rep("blue",12))
    d <- diagram(a,residue=TRUE,color=col,ylim=c(-6,-1),legend.x=NULL)
    legend("topleft",col=c("red","blue"),lty=1,
      legend=c("ribulose phosphate carboxylase",
      "acetyl-coenzyme A carboxylase"))
    title(main=paste("Calculated relative abundances of",
      "carboxylases from different organisms",sep="\n"))
    # ... now on to a species richness diagram
    # all the proteins, then rubisco and accoaco
    draw.diversity(d,"richness",logactmin=-3.6)
    draw.diversity(d,"richness",logactmin=-3.6,
      ispecies=1:12,col="red",add=TRUE)
    draw.diversity(d,"richness",logactmin=-3.6,
      ispecies=13:24,col="blue",add=TRUE)
    legend("bottomleft",col=c("red","blue","black"),lty=1,
      legend=c("ribulose phosphate carboxylase",
      "acetyl-coenzyme A carboxylase","all"))
    title(main=paste("Carboxylases with activities",
      "greater than 10^(-3.6)",sep="\n"))

    ## continuing from above ... make a rank-abundance
    ## diagram and fit with a lognormal distribution
    #if(require(vegan)) {
    #  basis("H2",-4)
    #  a <- affinity()
    #  logact <- diagram(a,residue=TRUE,do.plot=FALSE)$logact
    #  act <- 10^as.numeric(logact)
    #  # we use family=Gamma because our species have activities
    #  # (i.e., proportional abundances) and not integer counts
    #  mod <- rad.lognormal(act,family=Gamma)
    #  plot(mod,main=paste("Relative abundances of carboxylases",
    #    "fit with lognormal distribution",sep="\n"))
    #  # calculate Shannon diversity index
    #  # using revisit (CHNOSZ)
    #  H1 <- revisit(act,"shannon",as.is=TRUE)
    #  legend("topright",legend=paste("H =",round(H1,2)),pch=1)
    #  # using diversity (vegan)
    #  H2 <- diversity(matrix(act,nrow=1))
    #  stopifnot(isTRUE(all.equal(H1,H2)))
    #}

    ### using grep.file, read.fasta, add.protein
    # calculations for Pelagibacter ubique
    f <- system.file("HTCC1062.faa",package="CHNOSZ")
    # line numbers of all entries in the file
    j <- grep.file(f)  # length = 1354
    # locate entries whose names contain DNA
    j <- grep.file(f,"hypothetical")
    # get the amino acid compositions of these protein
    p <- read.fasta(f,j)
    # add these proteins to CHNOSZ's inventory
    i <- add.protein(p)
    # set up a thermodynamic system
    basis("CHNOS+")
    # calculate affinities in logfO2-pH space
    a <- affinity(H2O=c(-10,-2),O2=c(-82,-76),iprotein=i)
    # calculate the logarithms of activities
    d <- diagram(a,do.plot=FALSE,mam=FALSE)
    # show the protein richness
    draw.diversity(d,"richness",logactmin="")
    mtitle(c("Richness of hypothetical proteins in",
      expression(italic("Pelagibacter ubique"))))
    # show the coefficient of variation of activities
    draw.diversity(d,"cv")
    mtitle(c("Coefficient of variation of hypothetical",
      expression("protein activities in"~italic("P. ubique"))))
  } else if(which=="protein") {
    ## subcellular homologs of yeast glutaredoxin
    ## as a function of logfO2 - logaH2O, after Dick, 2009
    basis("CHNOS+")
    protein <- c("GLRX1","GLRX2","GLRX3","GLRX4","GLRX5")
    loc <- c("(C)","(M)","(N)","(N)","(M)")
    species(protein,"YEAST")
    t <- affinity(H2O=c(-10,0),O2=c(-85,-60))
    diagram(t,names=paste(protein,loc))
    title(main=paste("Yeast glutaredoxins (black) and residues (blue)\n",
      describe(thermo$basis[-c(2,5),])))
    # note the difference when we set as.residue=TRUE to
    # plot stability fields for the residue equivalents of the
    # proteins instead of the proteins themselves ...
    # the residue equivalent for one of the larger proteins appears
    diagram(t,names=paste(protein,loc),as.residue=TRUE,
      add=TRUE,col="blue")

    ## surface-layer proteins from Methanococcus and others:
    ## a speciation diagram for surface layer proteins
    ## as a function of oxygen fugacity after Dick, 2008
    # make our protein list
    organisms <- c("METSC","METJA","METFE","HALJP","METVO",
      "METBU","ACEKI","BACST","BACLI","AERSA")
    proteins <- c(rep("CSG",6),rep("SLAP",4))
    proteins <- paste(proteins,organisms,sep="_")
    # set some graphical parameters
    lwd <- c(rep(3,6),rep(1,4))
    lty <- c(1:6,1:4)
    # load the basis species and proteins
    basis("CHNOS+")
    species(proteins)
    # calculate affinities
    a <- affinity(O2=c(-100,-65))
    # make diagram
    d <- diagram(a,ylim=c(-5,-1),legend.x=NULL,lwd=lwd,
      ylab=as.expression(quote(log~italic(a[j]))),yline=1.7)
    # label diagram
    text(-80,-1.9,"METJA")
    text(-74.5,-1.9,"METVO")
    text(-69,-1.9,"HALJP")
    text(-78,-2.85,"METBU",cex=0.8,srt=-22)
    text(-79,-3.15,"ACEKI",cex=0.8,srt=-25)
    text(-81,-3.3,"METSC",cex=0.8,srt=-25)
    text(-87,-3.1,"METFE",cex=0.8,srt=-17)
    text(-79,-4.3,"BACST",cex=0.8)
    text(-85.5,-4.7,"AERSA",cex=0.8,srt=38)
    text(-87,-4.25,"BACLI",cex=0.8,srt=30)
    # add water line
    abline(v=-83.1,lty=2)
    title(main=paste("Surface-layer proteins",
      "After Dick, 2008",sep="\n"))

    ## relative metastabilities of bovine proteins, 
    ## as a function of temperature along a glutathione redox buffer
    mod.buffer("GSH-GSSG",c("GSH","GSSG"),logact=c(-3,-7))   
    basis(c("CO2","H2O","NH4+","SO4-2","H2","H+"),
      c(-1,0,-4,-4,"GSH-GSSG",-7)) 
    basis("CO2","gas")
    species(c("CYC","RNAS1","BPT1","ALBU","INS","PRIO"),"BOVIN")
    a <- affinity(T=c(0,200))
    diagram(a,as.residue=TRUE,ylim=c(-2,0.5))
    title(main="Bovine proteins")
  } else if(which=="buffer") {
    ## log fH2 - temperature diagram for mineral buffers
    ## and for given activities of aqueous CH2O and HCN
    ## After Schulte and Shock, 1995, Fig. 6: 
    ## 300 bars, log fCO2=1, log fN2=0, log aH2O=0
    # the mineral buffers FeFeO, QFM, PPM and HM are already
    # included in the thermo$buffer table so let's plot them.
    basis.logacts <- c(999,1,0,0,999,999,999)
    basis(c("Fe","CO2","H2O","nitrogen","hydrogen",
      "H2S","SiO2"),basis.logacts)
    basis(c("CO2","N2"),"gas")
    # initialize the plot
    xlim <- c(0,350)
    thermo.plot.new(xlim=xlim,ylim=c(-4,4),
      xlab=axis.label("T"),ylab=axis.label("H2"))
    res <- 50
    Tseq <- seq(xlim[1],xlim[2],length.out=res) 
    # a function to plot the log fH2 of buffers and label the lines
    logfH2plot <- function(buffer,lty,where) {
      basis("H2",buffer)
      t <- as.numeric(affinity(T=c(xlim,res),P=300,return.buffer=TRUE)$H2)
      lines(Tseq,t,lty=lty)
      # "where" is the percent distance along the x-axis to plot the label
      wherethis <- seq(xlim[1],xlim[2],length.out=100)[where]
      if(length(grep("_",buffer)) > 0) tt <- 
        thermo$buffer$logact[thermo$buffer$name==buffer] else tt <- buffer
      text(wherethis,splinefun(Tseq,t)(wherethis),tt)
    }
    # plot log fH2 of each mineral buffer
    logfH2plot("FeFeO",1,16)
    logfH2plot("QFM",1,30)
    logfH2plot("PPM",1,80)
    logfH2plot("HM",1,40)
    anotherplotfunction <- function(mybuff,lty,logact,where) {
      for(i in 1:length(logact)) {
        # update the species activity
        mod.buffer(mybuff,logact=logact[i])
        logfH2plot(mybuff,lty,where[i])
      }
    }
    # add and plot new buffers (formaldehyde and HCN)
    mod.buffer("mybuffer_1","formaldehyde","aq")
    logact <- c(-6,-10,-15); where <- c(10,10,25)
    anotherplotfunction("mybuffer_1",3,logact,where)
    mod.buffer("mybuffer_2","HCN","aq")
    where <- c(20,73,50)
    anotherplotfunction("mybuffer_2",2,logact,where)
    # title
    title(main=paste("Mineral buffers (solid), HCN (dashed), formaldehyde",
      "(dotted)\n",describe(thermo$basis[c(2,4),],T=NULL,P=300),
      "After Schulte and Shock, 1995"),cex.main=0.9)
  } else if(which=="affinity") {
    ## Activity of glycine as a function of those of
    ## HCN and formaldehyde (200 C, 300 bar)
    ## After Schulte and Shock, 1995, Fig. 5
    # we can define the basis as this:
    basis(c("formaldehyde","H2O","HCN","O2"))
    species("glycine")
    a <- affinity(HCHO=c(-10,-2,9),HCN=c(-18,-2,9),T=200,P=300)
    # that gave us *affinities* (dimensionless) for logact(glycine)=-3
    # (the default). we can now find the *activities* that
    # correspond to affinity=0
    logact.glycine <- species()$logact + a$values[[1]]
    # note transposition of the z-value matrix in the following command
    contour(x=-10:-2,y=seq(-18,-2,by=2),z=t(logact.glycine),
      xlab=axis.label("HCHO"),ylab=axis.label("HCN"),
      labcex=1,xaxs="i",yaxs="i")
    title(main=paste("log activity glycine, after Schulte and Shock, 1995",
      "200 deg C, 300 bar, logaH2O = 1",sep="\n"))

    ## amino acid synthesis at low and high temperatures
    ## after Amend and Shock, 1998
    # select the basis species and species of interest
    # and set their activities (first for the 18 degree C case)
    basis(c("H2O","CO2","NH4+","H2","H+","H2S"),
      log10(c(1,1e-4,5e-8,2e-9,5e-9,1e-15)))
    species(c("alanine","argininate","asparagine","aspartate","cysteine",
      "glutamate","glutamine","glycine","histidine","isoleucine",
      "leucine","lysinium","methionine","phenylalanine","proline",
      "serine","threonine","tryptophan","tyrosine","valine"),
      log10(c(3.9,0.7,1.1,3.3,0.5,3.8,1.0,5.8,1.2,0.7,
      0.8,1.0,2.8,0.5,0.5,4.6,5.8,0.6,0.9,2.8)/1e9))
    Tc <- 18
    T <- convert(Tc,"K")
    # converting A (dimensionless) to G of reaction (cal/mol) 
    # is like converting log K to standard G of reaction 
    AS98.18 <- 
      convert(convert(as.numeric(affinity(T=Tc)$values),"G",T=T),"J")/1000
    # the 100 degree C case
    Tc <- 100
    T <- convert(Tc,"K")
    basis(c("H2O","CO2","NH4+","H2","H+","H2S"),
      log10(c(1,2.2e-3,2.9e-6,3.4e-4,1.9e-6,1.6e-3)))
    species(1:20,log10(c(2.8e-9,5.0e-10,7.9e-10,2.4e-9,3.6e-10,
                         2.7e-9,7.2e-10,4.2e-9,8.6e-10,5.0e-10,
                         5.7e-10,7.2e-10,2.0e-9,3.6e-10,3.6e-10,
                         3.3e-9,4.2e-9,4.3e-10,6.5e-10,2.0e-9)))
    AS98.100 <- 
      convert(convert(as.numeric(affinity(T=Tc)$values),"G",T=T),"J")/1000
    # the nominal carbon oxidation state
    Z.C <- ZC(as.character(thermo$obigt$formula[thermo$species$ispecies]))
    # put them together
    print(data.frame(T100=AS98.100,T18=AS98.18,Z.C=Z.C))
    # values not exactly reproducing AS98 - different amino acid parameters
    # forget species to run next example
    species(delete=TRUE)

    ## affinities of metabolic reactions
    ## after Amend and Shock, 2001, Fig. 7
    basis(c("CO2","H2","NH3","O2","H2S","H+"))
    basis(c("O2","H2"),"aq")   # O2 and H2 were gas
    species("H2O")
    doplot <- function(T) {
      res <- 20
      a <- affinity(H2=c(-10,0,res),O2=c(-10,0,res),T=T)
      T.K <- convert(T,"K")   # temperature in Kelvin
      a <- convert(a$values[[1]],"G",T.K)  # affinities (cal/mol)
      a <- convert(a,"J")  # affinities (Joule)
      contour(x=seq(-10,0,length.out=res),
        y=seq(-10,0,length.out=res),z=t(a/1000),
        labcex=1,xlab=axis.label("H2"),ylab=axis.label("O2"))
    }
    layout(matrix(c(1,1,2,3,4,5),ncol=2,byrow=TRUE),heights=c(1,4,4))
    T <- c(25,55,100,150)
    par(mar=c(0,0,0,0))
    plot.new()
    text(0.5,0.1,paste(c("H2(aq) + 0.5O2(aq) = H2O(liq)\n\n",
      "after Amend and Shock, 2001")),cex=2)
    par(mar=c(3,3,0.5,0.5),cex=1.3,mgp=c(2,1,0))
    for(i in 1:length(T)) doplot(T[i])
    # this is so the plots in the next examples show up OK
    layout(matrix(1))

    ## continuation of last example, affinity calculations 
    ## in three dimensions
    print(affinity(H2=c(-10,0,3),O2=c(-10,0,3),T=c(25,150,4))$values)

    ## calculations on a transect
    # suppose that temperature and oxygen fugacity
    # both change in space (say from 1 to 6 meters),
    # that we have six values for each but want to
    # interpolate them to make a plot with smooth curves
    T <- splinefun(1:6,c(0,25,30,40,55,75))(seq(1,5,length.out=100))
    O2 <- splinefun(1:6,c(-90,-83,-78,-73,-68,-63))(seq(1,5,length.out=100))
    # what system could this be? 
    basis("CHNOS+")
    species(paste("CSG",c("METBU","METVO","METTL","METJA"),sep="_"))
    # now pass T and O2 to affinity, which because their lengths
    # are greater than three, treats them like coordinates for a
    # transect through chemical potential space rather than
    # the definition of a 2-dimensional grid
    a <- affinity(T=T,O2=O2)
    diagram(a,ylim=c(-4,-2))
    title(main=paste("Computed abundances of surface-layer proteins",
      "as a function of T and logfO2",sep="\n"))
  } else if(which=="subcrt") {
    ## heat capacity of Fe(cr)
    # compare calculated values of heat capacity with
    # values from Robie and Hemingway, 1995, (from which 
    # the parameters in the database were derived)
    nuts(c("J","K"))  # set the units
    # we set pressure here otherwise subcrt goes for 
    # Psat (saturation vapor pressure of H2O above 
    # 100 degrees C) which can not be calculated above 
    # the critical point of H2O (~647 K)
    t <- subcrt("Fe",T=seq(300,1800,20),P=1)
    plot(t$out[[1]]$T,t$out[[1]]$Cp,type="l",
      xlab=axis.label("T"),ylab=axis.label("Cp"))
    # add points from RH95
    RH95 <- thermo$expt$RH95
    points(RH95[,1],RH95[,2])
    title(main=paste("Heat capacity of Fe(cr)\n",
      "(points - Robie and Hemingway, 1995)"))
    # reset the units
    nuts(c("C","cal"))
    
    ## Skarn example after Johnson et al., 1992
    P <- seq(500,5000,500)
    # this is like the temperature specification used 
    # in the example by Johnson et al., 1992
    # T <- seq(0,1000,100)
    # we use this one to avoid warnings at 0 deg C, 5000 bar
    T <- c(2,seq(100,1000,100))
    subcrt(c("andradite","carbon dioxide","H2S","Cu+","quartz","calcite",
      "chalcopyrite","H+","H2O"),coeff=c(-1,-3,-4,-2,3,3,2,2,3),
      state=c("cr","g","aq","aq","cr","cr","cr","aq","liq"),
      P=P,T=T,grid="P")
    # the volumes are significantly different from SUPCRT92

    ## Standard Gibbs energy of reactions with HCN and 
    ## formaldehyde, after Schulte and Shock, 1995 Fig. 1
    rxn1 <- subcrt(c("formaldehyde","HCN","H2O","glycolic acid","NH3"),
      c(-1,-1,-2,1,1),P=300)
    rxn2 <- subcrt(c("formaldehyde","HCN","H2O","glycine"),
      c(-1,-1,-1,1),P=300)
    plot(x=rxn1$out$T,rxn1$out$G/1000,type="l",ylim=c(-40,-10),
      xlab=axis.label("T"),ylab=axis.label("DG0r","k"))
    lines(rxn1$out$T,rxn2$out$G/1000)
    # write the reactions on the plot
    text(150,-14,describe(rxn1$reaction,
      use.name=c(TRUE,FALSE,FALSE,TRUE,FALSE)))
    text(200,-35,describe(rxn2$reaction,
      use.name=c(TRUE,FALSE,FALSE,TRUE)))
    title(main=paste("Standard Gibbs energy of reactions",
      "after Schulte and Shock, 1995",sep="\n"))

    ## Calculation of chemical affinities
    # after LaRowe and Helgeson, 2007
    # Fig. 3 (a): reduction of nicotinamide adenine 
    # dinucleotide (NAD) coupled to oxidation of glucose
    # list the available NAD species
    info("NAD ")
    T <- seq(0,120,10)
    # oxidation of glucose (C6H12O6)
    basis(c("glucose","H2O","NH3","CO2","H+"),c(-3,0,999,-3,-7))
    t <- subcrt(c("NAD(ox)-","NAD(red)-2"),c(-12,12),logact=c(0,0),T=T)
    # LH07's diagrams are shown per mole of electron (24 e- per 12 NAD)
    A <- t$out$A/24/1000
    plot(x=T,y=A,xlim=range(T),ylim=c(1.4,5.4),
      xlab=axis.label("T"),ylab=axis.label("A",opt="k"),type="l")
    text("NAD(ox)-/NAD(red)-2 = 1",x=median(T),y=median(A))
    # different activity ratio
    t <- subcrt(c("NAD(ox)-","NAD(red)-2"),c(-12,12),logact=c(1,0),T=T)
    A <- t$out$A/24/1000
    lines(x=T,y=A)
    text("NAD(ox)-/NAD(red)-2 = 10",x=median(T),y=median(A))
    # different activity ratio
    t <- subcrt(c("NAD(ox)-","NAD(red)-2"),c(-12,12),logact=c(0,1),T=T)
    A <- t$out$A/24/1000
    lines(x=T,y=t$out$A/24/1000)
    text("NAD(ox)-/NAD(red)-2 = 0.1",x=median(T),y=median(A))
    # this command prints the reaction on the plot
    text(40,4.5,c2s(s2c(describe(t$reaction,
      use.name=c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)),
      sep="=",move.sep=TRUE),sep="\n"))
    # label the plot
    title(main=paste("Chemical affinity of NAD reduction",
     "after LaRowe and Helgeson, 2007",sep="\n"),
     sub=describe(thermo$basis,T=NULL))

    ### non-ideality calculations -- activity coefficients of 
    ### aqueous species as a function of charge, temperature,
    ### and ionic strength -- after Alberty, 2003 
    ## p. 16 Table 1.3  apparent pKa of acetic acid with
    ## changing ionic strength
    subcrt(c("acetic acid","acetate","H+"),c(-1,1,1),
      IS=c(0,0.1,0.25),T=25,property="logK")
    # note that these *apparent* values of G and logK approach
    # their *standard* counterparts as IS goes to zero.
    ## p. 95: basis and elemental stoichiometries of species 
    ## (a digression here from the nonideality calculations) 
    # note coefficient of O2 and NH3 will be zero for these species
    basis(c("ATP-4","H+","H2O","HPO4-2","O2","NH3"))
    # cf Eq. 5.1-33: (basis composition) 
    species(c("ATP-4","H+","H2O","HPO4-2","ADP-3","HATP-3","HADP-2","H2PO4-"))
    lb <- length(basis())
    # cf Eq. 5.1-32: (elemental composition)
    as.matrix(species()[,1:lb]) %*% as.matrix(basis()[,1:lb]) 
    ## p. 273-275: activity coefficient (gamma)
    ## as a function of ionic strength and temperature
    ## (as of 20080304, these do look quantitatively different 
    ## from the plots in Alberty's book.)
    iplotfun <- function(T,col,add=TRUE) {
      IS <- seq(0,0.25,0.0025)
      s <- subcrt(c("H2PO4-","HADP-2","HATP-3","ATP-4"),IS=IS,grid="IS",T=T)
      if(!add) thermo.plot.new(xlim=range(IS),ylim=c(0,1),
        xlab=axis.label("IS"),ylab="gamma")
      for(i in 1:4) lines(IS,10^s$out[[i]]$loggam,col=col)
    }
    iplotfun(0,"blue",add=FALSE)
    iplotfun(25,"black")
    iplotfun(40,"red")
    title(main=paste("activity coefficient (gamma) of -1,-2,-3,-4",
      "charged species at 0, 25, 40 deg C, after Alberty, 2003",
      sep="\n"),cex.main=0.95)
  }

}

.First.lib <- function(lib,pkg) {
  # version figuring adapted from package mgcv
  pkghelp <- library(help=CHNOSZ)$info[[1]]
  # things are different for older versions of R
  if(length(pkghelp)==1) pkghelp <- library(help=CHNOSZ)$info[[2]]
  version <- pkghelp[pmatch("Version:",pkghelp)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  date <- pkghelp[pmatch("Date:",pkghelp)]
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
    else add.protein("protein.csv")
  }
  if(file.exists('obigt.csv')) {
    tr <- try(rbind(read.csv('obigt.csv',as.is=TRUE),thermo$obigt),silent=TRUE)
    if(identical(class(tr),'try-error')) cat("thermo: obigt.csv in current directory is not compatible with thermo$obigt data table.\n")
    else add.obigt("obigt.csv")
  }
}


