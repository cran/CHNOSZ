## ----CHNOSZ_reset, include=FALSE----------------------------------------------
library(CHNOSZ)
reset()

## ----setfile, include=FALSE---------------------------------------------------
# Assign the file name to a variable and print the file name and number of species
setfile <- function(csvfile, dat=NULL) {
  # Assign csvfile outside this function
  assign("csvfile", csvfile, parent.frame())
  file <- system.file(paste0("extdata/OBIGT/", csvfile), package="CHNOSZ")
  dat <- read.csv(file, as.is=TRUE)
  ## Exclude entries for high-T polymorphs
  #dat <- dat[!dat$state %in% c("cr2", "cr3", "cr4", "cr5", "cr6", "cr7", "cr8", "cr9"), ]
  # The state and class of substance (used as section header), followed by number of species
  basename <- gsub(".csv", "", csvfile)
  class <- strsplit(basename, "_")[[1]][1]
  substr(class, 1, 1) <- toupper(substr(class, 1, 1))
  state <- strsplit(basename, "_")[[1]][2]
  if(identical(state, "aq")) state <- "Aqueous "
  else if(identical(state, "cr")) state <- "Solid "
  else if(identical(state, "gas")) state <- "Gas "
  else if(identical(state, "liq")) state <- "Liquid "
  else state <- "Optional "
  paste0(state, class, " (", nrow(dat), " species)")
}

## ----filerefs, include=FALSE--------------------------------------------------
filerefs <- function(csvfile, dat=NULL, message=FALSE) {
  # With dat, look for ref2 in dat
  whichref <- "ref2"
  # Without dat, look for ref1 in csvfile
  if(is.null(dat)) {
    file <- system.file(paste0("extdata/OBIGT/", csvfile), package="CHNOSZ")
    dat <- read.csv(file, as.is=TRUE)
    whichref <- "ref1"
  }
  ## Exclude entries for high-T polymorphs
  #dat <- dat[!dat$state %in% c("cr2", "cr3", "cr4", "cr5", "cr6", "cr7", "cr8", "cr9"), ]
  # Count number of times each reference is used
  tab <- table(dat[, whichref])
  # In case there are not references (previously for H2O_aq.csv) we return the species here
  if(length(tab)==0) return(paste(dat$name, dat$state))
  # The reference keys
  keys <- names(tab)
  # Warn if any keys aren't in thermo()$ref$key
  ikey <- match(keys, thermo()$ref$key)
  ina <- is.na(ikey)
  if(any(ina)) cat(paste("**WARNING: key(s)", paste(names(tab)[ina], collapse=" "), "not found in `thermo()$ref$key`**\n\n"))
  # Put the table in chronological order, according to thermo()$ref
  ikey <- order(match(keys, thermo()$ref$key))
  tab <- tab[ikey]
  keys <- keys[ikey]
  xxx <- lapply(seq_along(tab), function(i){
    thiskey <- keys[i]
    # Read thermo()$ref$note
    iref <- match(thiskey, thermo()$ref$key)
    note <- thermo()$ref$note[iref]
    # Show the note in italics
    if(!identical(note, "")) note <- paste0(" *", note, "* ")
    # Use bullets for ref2
    if(whichref=="ref2") bullet <- "- " else bullet <- ""
    # Convert key (e.g. LD12.2) to ref in OBIGT.bib (e.g. LD12)
    thisref <- gsub("\\..*$", "", thiskey)
    # Replace SLOP98 with slop98.dat, etc.
    # (we don't actually cite them here to keep the year from showing -- it's annoying to see e.g. "slop98.dat (1998)")
    citemark <- "@"
    if(thisref=="SLOP16") { thisref <- "slop16.dat"; citemark <- "" }
    if(thisref=="SLOP07") { thisref <- "slop07.dat"; citemark <- "" }
    if(thisref=="SLOP98") { thisref <- "slop98.dat"; citemark <- "" }
    if(thisref=="SPRONS92") { thisref <- "sprons92.dat"; citemark <- "" }
    if(thisref=="OBIGT") { thisref <- paste0("OBIGT (", thermo()$ref$year[iref], ")"); citemark <- "" }
    cat(bullet, citemark, thisref, " -- ", tab[i], note, "\n\n", sep="")
    # Get ref2 if we're in the outer list
    if(whichref!="ref2") filerefs(dat=dat[dat$ref1==names(tab)[i], ])
  })
  # Return all the species listed
  paste(dat$name, dat$state)
}

## ----used, include=FALSE------------------------------------------------------
# Initialize the list of used species
used <- character()
# Initialize the list of used optional species
optused <- character()

## ----reflist, results="asis", echo=FALSE--------------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE--------------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE--------------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE--------------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE--------------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE--------------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE--------------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE--------------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE--------------------------------------
used <- c(used, filerefs(csvfile))

## ----optreflist, results="asis", echo=FALSE-----------------------------------
optused <- c(optused, filerefs(csvfile))

## ----optreflist, results="asis", echo=FALSE-----------------------------------
optused <- c(optused, filerefs(csvfile))

## ----optreflist, results="asis", echo=FALSE-----------------------------------
optused <- c(optused, filerefs(csvfile))

## ----optreflist, results="asis", echo=FALSE-----------------------------------
optused <- c(optused, filerefs(csvfile))

## ----optreflist, results="asis", echo=FALSE-----------------------------------
optused <- c(optused, filerefs(csvfile))

## ----optreflist, results="asis", echo=FALSE-----------------------------------
optused <- c(optused, filerefs(csvfile))

