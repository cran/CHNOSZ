# CHNOSZ/util.element.R
# Copyright (C) 2006-2009 Jeffrey M. Dick
# functions to deal with chemical formulas

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

expand.formula <- function(elements,makeup) {
  # a function to expand formulas (makeups) to any
  # number of elements, i.e. 
  # expand.formula(c('C','H','O'),makeup('H2O')) returns c(0,2,1)
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

ZC <- function(x) {
  # calculate nominal carbon oxidation state of chemical formulas
  # FIXME: add more elements, warn about missing ones
  if(!is.null(names(x))) {
    # to deal with matrix or data frame arguments
    if(!is.data.frame(x)) x <- as.data.frame(x)
    # are carbon and other elements on rows or columns?
    iC <- match("C",colnames(x))
    if(is.na(iC)) {
      iC <- match("C",rownames(x))
      if(is.na(iC)) stop("carbon not found in the formula matrix!")
      else x <- t(x)
    }
    # carbon and other elements are now on the columns
    elements <- c("H","N","O","S","Z")
    charges <- c(-1,3,2,2,1)
    ielement <- match(elements,colnames(x))
    ina <- is.na(ielement)
    ielement <- ielement[!ina]
    elements <- elements[!ina]
    charges <- charges[!ina]
    xc <- x[,ielement,drop=FALSE]
    # sum of the charges
    xc <- rowSums(t(t(xc)*charges))
    # numbers of carbons
    nC <- x[,iC]
    # average oxidation number
    ZC <- as.numeric(xc/nC)
  } else {
    # to deal with character arguments
    # defactorize
    if(is.factor(x[1])) x <- as.character(x)
    # make sure all the elements we want are listed (even zeros)
    ZC <- numeric()
    for(i in 1:length(x)) {
      m <- as.data.frame(t(makeup(c(x[i],'C0H0N0O0S0Z0'))))
      ZC <- c(ZC,(-1*m$H+3*m$N+2*m$O+2*m$S+m$Z)/m$C)
    }
  }
  return(ZC)
}

## return a matrix with chemical formulas of residues
residue.formula <- function() {
  groups <- paste("[",colnames(thermo$protein)[6:25],"]",sep="")
  # formula of the sidechain groups
  f.groups <- as.character(thermo$obigt$formula[match(groups,thermo$obigt$name)])
  f.H2O <- "H2O"
  formulas <- c(f.H2O,f.groups)
  # an empty formula so that we have all the elements we need
  m.0 <- makeup("C0H0N0O0S0")
  # formula of the protein backbone group
  f.bb <- "C2H2NO"
  m.bb <- makeup(m.0,makeup(f.bb))
  for(i in 1:length(formulas)) {
    # sum the empty formula with that of the residue
    if(i==1) f <- makeup(m.0,makeup(formulas[i]))
    else f <- makeup(m.bb,makeup(formulas[i]))
    # create the output
    if(i==1) out <- f else out <- cbind(out,f)
  }
  out <- t(out)
  row.names(out) <- c(f.H2O,groups)
  out <- as.matrix(out)
  return(out)
}

protein.formula <- function(proteins,as.residue=FALSE) {
  # return a data frame with chemical formulas of proteins
  # proteins is a data frame in the format of thermo$protein
  rf <- residue.formula()
  out <- as.matrix(proteins[,5:25]) %*% as.matrix(rf)
  if(as.residue) out <- out/rowSums(proteins[,6:25])
  row.names(out) <- paste(proteins$protein,proteins$organism,sep="_")
  out <- as.data.frame(out)
  return(out)
}

