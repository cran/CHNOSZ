\name{info}
\alias{info}
\alias{info.character}
\alias{info.numeric}
\alias{info.approx}
\alias{info.text}
\title{Search the Thermodynamic Database}
\description{
  Search for species by name or formula, retrieve their thermodynamic properties and parameters, and add proteins to the thermodynamic database.
}

\usage{
  info(species = NULL, state = NULL, check.it=TRUE)
  info.character(species, state = NULL, check.protein=TRUE)
  info.numeric(ispecies, check.it=TRUE)
  info.approx(species, state = NULL)
  info.text(ispecies)
}

\arguments{
  \item{species}{character, names or formulas of species, or (for \code{info} only) numeric with same meaning as \code{ispecies}}
  \item{state}{character, physical states of the species}
  \item{check.it}{logical, check GHS and EOS parameters for self-consistency?}
  \item{check.protein}{logical, check if a matching protein can be found?}
  \item{ispecies}{numeric, index of species in the thermodynamic database}
}

\details{
  \code{info} is the primary function used for querying the thermodynamic database (\code{\link{thermo}$obigt}). In common usage, it is called recursively; first with a character value (or values) for \code{species} indicating the name(s) or formula(s) of the species of interest. The result of this call is a numeric value, which can be provided as an argument in a second call to \code{info} in order to retrieve a data frame of the thermodynamic properties of the species. For its work, \code{info} calls on the other functions that are described below, which unlike \code{info} all expect arguments with length=1.

  \code{info.character} searches for matches of the indicated \code{species} to names, chemical formulas, and abbreviations (in the \samp{abbrv} column) in the thermodynamic database. If the text of the \code{species} is matched the index of that species is returned. If there are multiple matches for the \code{species} and \code{state} is NULL, the index of first match is returned. The order of entries in \code{thermo$obigt} is grouped by states in the order \samp{aq}, \samp{cr}, \samp{gas}, \samp{liq}, so for species in both aqueous and gaseous states the index of the aqueous species is returned, unless \code{state} is set to \samp{gas}. The two exceptions are species identified by \samp{O2} or \samp{oxygen} (which without any indicated \code{state} matches the gaseous species) and \samp{H2O} (which matches the liquid species even if the indicated state is \samp{aq}). Normally, if a species match can not be located, the function then looks for proteins with the name of \code{species} (using \code{\link{iprotein}}), computes its properties if found (\code{\link{ip2aa}}) and adds this to the thermodynamic database (\code{\link{mod.obigt}}). \code{check.protein} prevents the processing of proteins and is provided to avoid an infinite loop in the interaction with \code{mod.obigt}.

  \code{info.character} has additional logic for dealing with proteins and with multiple matches for the \samp{cr} state. If the \code{state} is \samp{cr}, matches will be counted for states entered as \samp{cr1}, \samp{cr2} etc in the database, and all of the species indices will be returned. Note, however, that \code{info} only ever returns a single species index, which becomes NA in the case of multiple matches to \samp{cr}. This functionality of \code{info.character} is used in \code{\link{subcrt}} to handle minerals with phase transitions. 

  Names of species including an underscore character are indicative of proteins, e.g. \samp{LYSC_CHICK}. If the name of a protein is provided to \code{info.character} and the composition of the protein can be found using \code{\link{protein}}, the thermodyamic properties and parameters of the nonionized protein (calculated using amino acid group additivity) are added to the thermodynamic database. Included in the return value, as for other species, is the index of the protein in the thermodynamic database or \code{NA} if the protein is not found. Names of proteins and other species can be mixed.

  \code{info.approx} searches the database for similar names or formulas using \code{\link{agrep}}. If one or more of these is found, the results are summarized on the screen, and the indices of the approximately matching species are returned. Species that have no approximate matches are indicated by \code{NA} in the return value. When invoked by \code{info}, the latter function accepts the species index only for a single approximate match; multiple matches are translated to NA.

  \code{info.numeric} returns the rows of \code{thermo$obigt} indicated by \code{ispeices}, after removing any order-of-magnitude scaling factors. If these species are all aqueous or are all not aqueous, the compounded column names used in \code{thermo$obigt} are replaced with names appropriate for the corresponding equations of state. A missing value of one of the standard molal Gibbs energy (\code{G}) or enthalpy (\code{H}) of formation from the elements or entropy (\code{S}) is calculated from the other two, if available. If \code{check.it} is TRUE, several checks of self consistency among the thermodynamic properties and parameters are performed using \code{\link{checkGHS}} and \code{\link{checkEOS}} (this depends on the completeness of the data entry).


}


\seealso{
  \code{\link{thermo}} for the thermodynamic database (specifically, \code{thermo$obigt}). \code{\link{check.obigt}} for checking self-consistency of individual entries in the database. \code{\link{protein}} for gathering compositions and thermodynamic properties of proteins. 
}


\examples{
\dontshow{data(thermo)}
## summary of available data
info()
\dontrun{
## run a consistency check on each species in the database
# (marked dontrun because it takes a while)
info(check=TRUE) }

## species information
# search for something named (or whose formula is) "Fe"
si <- info("Fe")
# use the number to get the full record
info(si)
# it is possible to get a range of records
info(si:(si+3))

## dealing with states
# default order of preference for names: aq > gas > cr,liq
info(c("methane","ethanol","glycinate"))  # aq, aq, aq
info(c("adenosine","alanine","hydroxyapatite"))  # aq, aq, cr
# state argument overrides the default
info(c("ethanol","adenosine"),state=c("gas","cr"))
# formulas default to aqueous species, if available
info(c("CH4","CO2","CS2","MgO"))  # aq, aq, gas, cr
# state argument overrides the default
info(c("CH4","CO2","MgO"),"gas")  # gas, gas, NA
# exceptions to the aqueous default is O2
info("O2")  # gas

## partial name or formula searches
info("ATP")
info("thiol")
info("MgC")
# add an extra character to refine a search
# or to search using terms that have exact matches
info("MgC ")
info("acetate ")
info(" H2O")
}

\keyword{primary}
