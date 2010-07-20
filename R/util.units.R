# CHNOSZ/util.units.R
# Copyright (C) 2006-2009 Jeffrey M. Dick
# convert units and set preferred units 

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

