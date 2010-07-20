# CHNOSZ/thermo.R
# Copyright (C) 2006-2009 Jeffrey M. Dick
# utility functions for the CHNOSZ package
# speciate/thermo.R 20051021 jmd

protein.length <- function(protein) {
  # character, name of protein
  # positive numeric, index of protein in thermo$species
  # negative numeric, index of protein in thermo$protein
  # dataframe, amino acid composition in format of thermo$protein
  # 20090331 added negative numeric option and multicore support
  # 20100329 added support for dataframe

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
  if(is.data.frame(protein)) return(as.numeric(rowSums(protein[,6:25])))
  else return(as.numeric(mylapply(protein,x2length)))
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
  ## load data files in user's directory
  #if(file.exists('protein.csv')) {
  #  tr <- try(rbind(read.csv('protein.csv',as.is=TRUE),thermo$protein),silent=TRUE)
  #  if(identical(class(tr),'try-error')) cat("thermo: protein.csv in current directory is not compatible with thermo$protein data table.\n")
  #  else add.protein("protein.csv")
  #}
  #if(file.exists('obigt.csv')) {
  #  tr <- try(rbind(read.csv('obigt.csv',as.is=TRUE),thermo$obigt),silent=TRUE)
  #  if(identical(class(tr),'try-error')) cat("thermo: obigt.csv in current directory is not compatible with thermo$obigt data table.\n")
  #  else add.obigt("obigt.csv")
  #}
}

