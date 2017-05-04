## ----library_CHNOSZ, include=FALSE---------------------------------------
library(CHNOSZ)

## ----data_thermo, include=FALSE------------------------------------------
data(thermo)

## ----nspecies, include=FALSE---------------------------------------------
# assign the file name to a variable and print the file name and number of species
setfile <- function(csvfile, dat=NULL) {
  # assign csvfile outside this function
  assign("csvfile", csvfile, parent.frame())
  file <- system.file(paste0("extdata/OBIGT/", csvfile, ".xz"), package="CHNOSZ")
  dat <- read.csv(file, as.is=TRUE)
  # the class of substance, followed by number of species
  class <- gsub("_.*", "", csvfile)
  substr(class, 1, 1) <- toupper(substr(class, 1, 1))
  paste0("`", class, "` (", nrow(dat), ")")
}

## ----filerefs, include=FALSE---------------------------------------------
filerefs <- function(csvfile, dat=NULL, message=FALSE) {
  # with dat, look for ref2 in dat
  whichref <- "ref2"
  # without dat, look for ref1 in csvfile
  if(is.null(dat)) {
    file <- system.file(paste0("extdata/OBIGT/", csvfile, ".xz"), package="CHNOSZ")
    dat <- read.csv(file, as.is=TRUE)
    whichref <- "ref1"
  }
  # count number of times each reference is used
  tab <- table(dat[, whichref])
  # there are no references in H2O_aq.csv so we return the species here
  if(length(tab)==0) return(paste(dat$name, dat$state))
  # the keys only (no [S92] etc.)
  keys <- sapply(strsplit(names(tab), " "), "[", 1)
  # warn if any keys aren't in thermo$ref$key
  ikey <- match(keys, thermo$ref$key)
  ina <- is.na(ikey)
  if(any(ina)) cat(paste("**WARNING: key(s)", paste(names(tab)[ina], collapse=" "), "not found in `thermo$ref$key`**\n\n"))
  # put the table in chronological order, according to thermo$ref
  #ikey <- na.omit(match(thermo$ref$key, keys))
  ikey <- order(match(keys, thermo$ref$key)) # works for duplicated keys (e.g. "Sho92" and "Sho92 [S98]")
  tab <- tab[ikey]
  keys <- keys[ikey]
  xxx <- lapply(seq_along(tab), function(i){
    thiskey <- keys[i]
    # read thermo$ref$note
    iref <- match(thiskey, thermo$ref$key)
    note <- thermo$ref$note[iref]
    if(!identical(note, "")) note <- paste0(" *", note, "* ")
    # append symbol for [S92], [S98], or [S15]
    if(grepl("[S92]", names(tab)[i], fixed=TRUE) | grepl("SPRONS92", names(tab)[i], fixed=TRUE)) note <- paste0(note, "(ø)")
    if(grepl("[S98]", names(tab)[i], fixed=TRUE) | grepl("SLOP98", names(tab)[i], fixed=TRUE)) note <- paste0(note, "(\\*)")
    if(grepl("[S07]", names(tab)[i], fixed=TRUE) | grepl("SLOP07", names(tab)[i], fixed=TRUE)) note <- paste0(note, "(†)")
    if(grepl("[S15]", names(tab)[i], fixed=TRUE) | grepl("SLOP15", names(tab)[i], fixed=TRUE)) note <- paste0(note, "(‡)")
    # use bullets for ref2
    if(whichref=="ref2") bullet <- "- " else bullet <- ""
    # convert key (e.g. LD12.2) to ref in OBIGT.bib (e.g. LD12)
    thisref <- gsub("\\..*$", "", thiskey)
    # replace SLOP98 with slop98.dat, etc.
    # (we don't actually cite them here to keep the year from showing -- it's annoying to see e.g. "slop98.dat (1998)")
    citemark <- "@"
    if(thisref=="SLOP15") { thisref <- "slop15.dat"; citemark <- "" }
    if(thisref=="SLOP07") { thisref <- "slop07.dat"; citemark <- "" }
    if(thisref=="SLOP98") { thisref <- "slop98.dat"; citemark <- "" }
    if(thisref=="SPRONS92") { thisref <- "sprons92.dat"; citemark <- "" }
    if(thisref=="CHNOSZ") { citemark <- "" }
    cat(bullet, citemark, thisref, " -- ", tab[i], note, "\n\n", sep="")
    # get ref2 if we're in the outer list
    if(whichref!="ref2") filerefs(dat=dat[dat$ref1==names(tab)[i], ])
  })
  # return all the species listed
  paste(dat$name, dat$state)
}

## ----H2O_aq, results="asis", echo=FALSE----------------------------------
cat('This file contains H<sub>2</sub>O, *e*<sup>-</sup>, and H<sup>+</sup>. The properties of H<sub>2</sub>O are listed as NA; CHNOSZ calculates its properties using a Fortran subroutine taken from SUPRCT92 [@JOH92]. The properties of the proton (H<sup>+</sup>) are 0. The properties of the electron (*e*<sup>-</sup>) are 0, except for *S*°, which is the opposite of *S*° for the "element" of charge, Z (see `?thermo`).\n')

## ----used, include=FALSE-------------------------------------------------
# initialize the list of used species
used <- character()

## ----reflist, results="asis", echo=FALSE---------------------------------
used <- c(used, filerefs(csvfile))

## ----inorganic_aq, results="asis", echo=FALSE----------------------------
cat("Note: ZnCl4-2 was present in sprons92.dat but not in slop98.dat or later files, and is not included in CHNOSZ.")

## ----reflist, results="asis", echo=FALSE---------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE---------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE---------------------------------
used <- c(used, filerefs(csvfile))

## ----CHNOSZ_aq, results="asis", echo=FALSE-------------------------------
cat("The primary aqueous silica species in CHNSOZ is SiO<sub>2</sub> [@SHS89]. The pseudospecies H<sub>4</sub>SiO<sub>4</sub> is used to make activity diagrams with *a*<sub>H<sub>4</sub>SiO<sub>4</sub></sub> as a variable. The GHS and HKF parameters for this pseudospecies were calculated using CHNOSZ; see the vignette [*Regressing thermodynamic data*](eos-regress.html) for more information.")

## ----reflist, results="asis", echo=FALSE---------------------------------
used <- c(used, filerefs(csvfile))

## ----inorganic_cr, results="asis", echo=FALSE----------------------------
cat("Note: chamosite,7A and witherite were present in sprons92.dat but not in slop98.dat or later files, and are not included in CHNOSZ.\n\n")
cat("Note: parameters used here for goethite differ slightly from those listed in the slop files [@Sho09].")

## ----reflist, results="asis", echo=FALSE---------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE---------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE---------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE---------------------------------
used <- c(used, filerefs(csvfile))

## ----reflist, results="asis", echo=FALSE---------------------------------
used <- c(used, filerefs(csvfile))

## ----H2O_aq, results="asis", echo=FALSE----------------------------------
cat('This file contains H<sub>2</sub>O, *e*<sup>-</sup>, and H<sup>+</sup>. The properties of H<sub>2</sub>O are listed as NA; CHNOSZ calculates its properties using a Fortran subroutine taken from SUPRCT92 [@JOH92]. The properties of the proton (H<sup>+</sup>) are 0. The properties of the electron (*e*<sup>-</sup>) are 0, except for *S*°, which is the opposite of *S*° for the "element" of charge, Z (see `?thermo`).\n')

## ----used2, include=FALSE------------------------------------------------
# initialize the list of used species
used2 <- character()

## ----reflist2, results="asis", echo=FALSE--------------------------------
used2 <- c(used2, filerefs(csvfile))

## ----inorganic_aq, results="asis", echo=FALSE----------------------------
cat("Note: ZnCl4-2 was present in sprons92.dat but not in slop98.dat or later files, and is not included in CHNOSZ.")

## ----reflist2, results="asis", echo=FALSE--------------------------------
used2 <- c(used2, filerefs(csvfile))

## ----reflist2, results="asis", echo=FALSE--------------------------------
used2 <- c(used2, filerefs(csvfile))

## ----reflist2, results="asis", echo=FALSE--------------------------------
used2 <- c(used2, filerefs(csvfile))

## ----CHNOSZ_aq, results="asis", echo=FALSE-------------------------------
cat("The primary aqueous silica species in CHNSOZ is SiO<sub>2</sub> [@SHS89]. The pseudospecies H<sub>4</sub>SiO<sub>4</sub> is used to make activity diagrams with *a*<sub>H<sub>4</sub>SiO<sub>4</sub></sub> as a variable. The GHS and HKF parameters for this pseudospecies were calculated using CHNOSZ; see the vignette [*Regressing thermodynamic data*](eos-regress.html) for more information.")

## ----reflist2, results="asis", echo=FALSE--------------------------------
used2 <- c(used2, filerefs(csvfile))

## ----inorganic_cr, results="asis", echo=FALSE----------------------------
cat("Note: chamosite,7A and witherite were present in sprons92.dat but not in slop98.dat or later files, and are not included in CHNOSZ.\n\n")
cat("Note: parameters used here for goethite differ slightly from those listed in the slop files [@Sho09].")

## ----reflist2, results="asis", echo=FALSE--------------------------------
used2 <- c(used2, filerefs(csvfile))

## ----reflist2, results="asis", echo=FALSE--------------------------------
used2 <- c(used2, filerefs(csvfile))

## ----reflist2, results="asis", echo=FALSE--------------------------------
used2 <- c(used2, filerefs(csvfile))

## ----reflist2, results="asis", echo=FALSE--------------------------------
used2 <- c(used2, filerefs(csvfile))

## ----reflist2, results="asis", echo=FALSE--------------------------------
used2 <- c(used2, filerefs(csvfile))

## ----check_used_used2, results="asis", echo=FALSE------------------------
if(length(used) != length(used2)) cat(paste0("**WARNING: Tabbed list has ", length(used), " species but 'All at once list' has ", length(used2),".**\n\n"))

