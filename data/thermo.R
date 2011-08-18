# CHNOSZ/data/thermo.R
# create the "thermo" data object from data files


thermo <- list(
  # as.is: keep character values as character and not factor
  opt = as.list(read.csv("opt.csv",as.is=TRUE)),
  element = read.csv("element.csv",as.is=1:3),
  obigt = read.csv("OBIGT.csv",as.is=1:7),
  refs = read.csv("refs.csv",as.is=TRUE),
  buffers = read.csv("buffer.csv",as.is=1:3),
  protein = read.csv("protein.csv",as.is=1:4),
  stress = read.csv("stress.csv",check.names=FALSE,as.is=TRUE),
  groups = read.csv("groups.csv",row.names=1,check.names=FALSE),
  basis = NULL,
  species = NULL,
  water = NULL,
  water2 = NULL
)

packageStartupMessage(paste("thermo$obigt has",
  nrow(thermo$obigt[thermo$obigt$state=="aq",]),
  "aqueous,",nrow(thermo$obigt),"total species"))

#cat(paste("thermo: loaded",nrow(thermo$ECO),"proteins to thermo$ECO\n"))
#cat(paste("thermo: loaded",nrow(thermo$SGD),"proteins to thermo$SGD\n"))
#cat(paste("thermo: loaded",nrow(thermo$yeastgfp),"localizations and",
#  length(thermo$yeastgfp$abundance[!is.na(thermo$yeastgfp$abundance)]),"abundances to thermo$yeastgfp\n"))

