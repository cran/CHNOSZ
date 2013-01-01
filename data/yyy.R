# CHNOSZ/data/yyy.R
# create the "thermo" data object from data files
# to be called by thermo.R (using sys.source)

thermo <- list(
  # as.is: keep character values as character and not factor
  opt = as.list(read.csv("opt.csv", as.is=TRUE)),
  element = read.csv("element.csv", as.is=1:3),
  obigt = read.csv("OBIGT.csv", as.is=1:7),
  refs = read.csv("refs.csv", as.is=TRUE),
  buffers = read.csv("buffer.csv", as.is=1:3),
  protein = read.csv("protein.csv", as.is=1:4),
  groups = read.csv("groups.csv", row.names=1, check.names=FALSE),
  basis = NULL,
  species = NULL,
  water = NULL,
  water2 = NULL,
  Psat = NULL,
  opar = NULL
)
