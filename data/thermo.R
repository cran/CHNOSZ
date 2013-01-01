# CHNOSZ/data/thermo.R
# load the "thermo" data object into an environment named CHNOSZ:thermo

# data(thermo) can be called multiple times during an interactive session,
# first call is at package loading (by .onAttach); subsequent calls are from the user
if("CHNOSZ:thermo" %in% search()) {
  # the environment already exists; restore the "thermo" object to default values
  source("xxx.R") 
} else {
  # the environment doesn't exist; create it and the "thermo" object
  sys.source("yyy.R", attach(NULL, name="CHNOSZ:thermo"))
}

# give a summary of some of the data
packageStartupMessage(paste("thermo$obigt:",
  nrow(thermo$obigt[thermo$obigt$state=="aq",]),
  "aqueous,", nrow(thermo$obigt), "total species"))

# note if there are duplicated species
idup <- duplicated(paste(thermo$obigt$name, thermo$obigt$state))
if(any(idup)) warning("thermo$obigt: duplicated species: ", 
  paste(thermo$obigt$name[idup], "(", thermo$obigt$state[idup], ")", sep="", collapse=" "))
rm(idup)
