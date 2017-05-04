# CHNOSZ/util.args.R
# functions to create argument lists and get name of calling function

### unexported functions ###

# These functions are used to normalize user-input arguments, which are case-insensitive.

# return a list with elements:
#   `prop` for all the properties available for the specified equations-of-state
#   `prop` for the lower-case version of property
#   `Prop` for the upper-case (of first letter) version of property
# produces an error if any of `property` is not in the list of available properties.
# (See water() and subcrt() for the available properties for different species.)
eos.args <- function(eos='',property=NULL,T=NULL,P=NULL) {
  # the available properties for supcrt, probably
  props <- c('G','H','S','Cp','V','kT','E')
  if(eos=='water') {
    # things we also get with water
    props <- c(props,'A','U','Cv','Psat','rho','Q','X','Y','epsilon','w')
    # they keep on coming: things we also get with SUPCRT92
    if(get("thermo")$opt$water == "SUPCRT92")
      props <- c(props,'Z','visc','tcond','tdiff','Prndtl','visck','albe','daldT','alpha','beta')
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

# force T and P to equal length
# also looks for the keyword Psat in the value of P and substitutes calculated values of the saturation vapor pressure
TP.args <- function(T=NULL, P=NULL) {
  # keep the [1] here because some functions (e.g. subcrt) will repeat "Psat"
  if(identical(P[1], "Psat")) {
    P <- water("Psat", T, P="Psat")[, 1]
    # water.SUPCRT92 issues its own warnings about 
    # exceeding Psat's temperature limit
    if(get("thermo")$opt$water == "IAPWS95")
      if(length(which(is.na(P)))>0) 
        warning('TP.args: NAs in Psat (likely T > Tc where Tc = 647.096 K)',call.=FALSE)
  }
  if(length(P) < length(T) & !is.null(P)) P <- rep(P, length.out=length(T))
  else if(length(T) < length(P) & !is.null(T)) T <- rep(T, length.out=length(P))
  # something we do here so the SUPCRT water calculations work
  T[T==273.15] <- 273.16
  return(list(T=T, P=P))
}

# make the argument lowercase, then transform a, c, g, and l to aq, gas, cr, and liq
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

caller.name <- function(n=2) {
  # returns the name of the calling function n frames up
  # (n=2: the caller of the function that calls this one)
  # or character() if called interactively
  if(sys.nframe() < n) name <- character()
  else {
    sc <- sys.call(-n)[[1]]
    name <- try(as.character(sc),silent=TRUE)
    # also return character() if the value from sys.call is
    # the function itself (why does this sometimes happen,
    # e.g. when called from affinity()?)
    if(class(name)=="try-error") name <- character()
  }
  return(name)
}
