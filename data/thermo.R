# vim:syntax=r
# CHNOSZ/data/thermo.R
# create the "thermo" data object from text files

# assemble the data object from text files.

# yeastgfp preprocessing
ygfp <- read.csv('yeastgfp.csv')
# convert factors to numeric w/o NA coercion warnings
warn <- options(warn=-1)
ygfp$abundance <- as.numeric(as.character(ygfp$abundance))
options(warn)

thermo <- list(
  # as.is: suppressing conversion of character columns 
  # to factors makes it easier to add new values later 
  # and just makes everything easier
  opt = as.list(read.csv("opt.csv",as.is=TRUE)),
  element = read.csv("element.csv",as.is=1:3),
  obigt = read.csv("OBIGT.csv",as.is=1:7),
  source = read.csv("source.csv",as.is=TRUE),
  buffer = read.csv("buffer.csv",as.is=1:3),
  protein = read.csv("protein.csv",as.is=1:4),
  stress = read.csv("stress.csv",check.names=FALSE,as.is=TRUE),
  groups = read.csv("groups.csv",row.names=1,check.names=FALSE),
  basis = NULL,
  species = NULL,
  water = NULL,
  water2 = NULL,
  expt = list(
    PM90 = read.csv("PM90.csv"),
    RH95 = read.csv("RH95.csv"),
    RT71 = read.csv("RT71.csv")
  ),
  SGD = read.csv("SGD.csv",as.is=TRUE),
  ECO = read.csv("ECO.csv",as.is=TRUE),
  yeastgfp = ygfp
)

# yeastgfp cleanup
rm(ygfp)
rm(warn)

cat(paste("thermo: loaded",nrow(thermo$obigt[thermo$obigt$state=="aq",]),"aqueous,",nrow(thermo$obigt),
  "total species to thermo$obigt\n"))
cat(paste("thermo: loaded",nrow(thermo$ECO),"proteins to thermo$ECO\n"))
cat(paste("thermo: loaded",nrow(thermo$SGD),"proteins to thermo$SGD\n"))
cat(paste("thermo: loaded",nrow(thermo$yeastgfp),"localizations and",
  length(thermo$yeastgfp$abundance[!is.na(thermo$yeastgfp$abundance)]),"abundances to thermo$yeastgfp\n"))


