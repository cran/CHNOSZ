\name{info}
\alias{info}
\title{Search the Thermodynamic Database}
\description{
  Search for species by name or formula, and retrieve their thermodynamic properties and parameters.
}

\usage{
  info(species = NULL, states = NULL, quiet = FALSE, return.approx = TRUE)
}

\arguments{
  \item{species}{character, names or formulas of species, or numeric, indices of species in the thermodynamic database.}
  \item{states}{numeric, physical states of the species.}
  \item{quiet}{logical, produce fewer messages and make fewer tests?}
  \item{return.approx}{logical, return the indices of approximately matching species?}
}

\details{

  This function searches the names or chemical formulas of \code{species} in the thermodynamic database (\code{\link{thermo}$obigt}). Searches can optionally be limited to certain physical \code{state}s (among \samp{aq}, \samp{cr}, \samp{gas}, \samp{liq} as well as \samp{cr1}, \samp{cr2} etc.) For each of the \code{species} that is matched, the index of that species is appended to the return value. If an exact match of a species is not located, \code{info} searches the database for similar names or formulas; if any of these are found, the results are summarized on the screen, and the indices of the approximately matching species are included in the return value (if \code{return.approx} is TRUE). Species that have no exact or approximate matches are indicated by \code{NA} in the return value. 

  If \code{states} is NULL (the default) and there are multiple matches for the name of a species, the species whose state is that in \code{thermo$opt$state} (\samp{aq} by default) is selected, or failing that, the first matching species is taken. The exceptions to this rule is \samp{O2}, for which the state defaults to \samp{gas}.

  If \code{species} is numeric instead of character, the corresponding rows of the thermodynamic database are returned, after removing any order-of-magnitude scaling factors. If these species are all aqueous or are all not aqueous, the compounded column names used in \code{thermo$obigt} are replaced with names appropriate for the corresponding equations of state.

  The names of proteins are distinguished from those of other species insofar as they contain an underscore character, as in \samp{LYSC_CHICK}. If the name of a protein is provided to \code{info} and the composition of the protein can be found using \code{\link{protein}}, the thermodyamic properties and parameters of the nonionized protein (calculated using \code{protein}) are added to the thermodynamic database. Included in the return value, as for other species, is the index of the protein in the thermodynamic database or \code{NA} if the protein is not found. Names of proteins and other species can be mixed.

  If \code{quiet} is \code{FALSE}, several checks of self consistency among the thermodynamic properties and parameters are performed. A missing value of one of the standard molal Gibbs energy (\code{G}) or enthalpy (\code{H}) of formation from the elements or entropy (\code{S}) is calculated from the other two, if available. In addition, the value of \code{G} calculated from those of \code{H} and \code{S} is compared with the tabulated value and a message printed if the difference exceeds 1000 cal mol\eqn{^{-1}}{^-1}. Likewise, if the values of standard molal heat capacity or volume calculated from the equations-of-state parameters differ (by greater than 1 cal mol K\eqn{^{-1}}{^-1} or 1 cm\eqn{^3} mol\eqn{^{-1}}{^-1}) from the tabulated values, messages to this effect are printed. 

}


\value{

  If \code{species} is character of length \eqn{n} and no approximate matches are found, the \code{\link{invisible}} return is a numeric vector of length \eqn{n} containing for each species its index in the thermodynamic database or \code{NA}. If two or more approximate matches are located for any of the species, the return value is instead a list of length \eqn{n}, each element of which contains the matching index or indices of the species or \code{NA}. If \code{species} is numeric, the corresponding rows of the thermodynamic database are visibly returned, after removing order-of-magnitude multipliers. 
  
}

\seealso{
  \code{\link{thermo}} for the thermodynamic database (specifically, \code{thermo$obigt}). \code{\link{check.obigt}} for checking self-consistency of individual entries in the database. \code{\link{protein}} for gathering compositions and thermodynamic properties of proteins. 
}


\examples{
  \dontshow{data(thermo)}
  ## basic operation
  \dontrun{
  # run a consistency check on each species in the database
  # (marked dontrun because it takes a while)
  info(check=TRUE)
  }

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
