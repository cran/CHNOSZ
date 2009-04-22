\name{species}
\alias{species}
\title{Species of Interest}
\description{
  Define the species of interest in a system; modify their physical states and logarithms of activities.
}

\usage{
  species(species = NULL, state = NULL, delete = FALSE, quiet = TRUE)
}

\arguments{
  \item{species}{character, names or formulas of species to add to the species definition; numeric, rownumbers of species to modify or delete.}
  \item{state}{character, physical states; numeric, logarithms of activities or fugacities.}
  \item{delete}{logical, delete the species identified by numeric values of \code{species} (all species if that argument is missing)?}  
  \item{quiet}{logical, make less noise?}
}

\details{

  Only after defining the \code{\link{basis}} species of your system can you use this function to identify for the program the species of interest. The primary effect of this function is to add to a globally accessible dataframe (\code{\link{thermo}$species}) identifying information about the species you have selected and the coefficients in the formation reactions of the species from the basis species.

  The command to add species operates on one or more character values of \code{species}. For each match of \code{species} (optionally restricted to a \code{state} among \samp{aq}, \samp{cr}, \samp{gas}, \samp{liq}) to the name of a species or a formula or abbreviation in the thermodynamic database, a row is added to the dataframe in \code{thermo$species}.  This dataframe holds the identifying characteristics of the species (including the logarithms of their activities or fugacities) as well as the stoichiometric reaction coefficients for the formation of each of the species from the basis species. The default values for logarithms of activities are \eqn{-3} for aqueous (\samp{aq}) species and \eqn{0} for others. If the basis species do not contain all the elements in any of \code{species}, a message to this effect is generated but the species is nevertheless added with its incomplete stoichiometric definition, even if all zeros.

  If \code{state} is \code{NULL} (the default), species in any state can be matched in the thermodynamic database. If there are multiple matches for a species, the one that is in the state given by \code{thermo$opt$state} is chosen, otherwise the first matching species is used. Note that the \code{state}s of species representing phases of minerals that undergo phase transitions are coded as \samp{cr1}, \samp{cr2}, \samp{cr3}, \code{...} (phases with increasing temperature). If \code{state} is \samp{cr} when one of these minerals is matched, all the species/phases are added, otherwise only the one corresponding to \samp{cr1} (or a higher-temperature phase if that was specified).

  The command to modify the logarithms of activities of species (logarithms of fugacities for gases) operates on one or more numeric values of \code{species} referring to the rownumbers of the species dataframe. This regards the numeric values in \code{state} as the logarithms of activities; if instead the values in \code{state} are character the function attempts to appropriately alter the state(s) of the species identified in the first argument. If \code{species} is numeric and \code{delete} is \code{TRUE}, the rows representing these species are deleted from the dataframe; if the only argument is \code{delete} and it is \code{TRUE}, all the species are removed.

  The \code{quiet} option turns off most, but not all, printed messages. In particular, the display of the species dataframe after an update remains in effect.
}

\value{

  With no arguments, the value of \code{thermo$species}. Successful addition of species returns the rownumbers of the species added to \code{thermo$species}; for every species that is not located in the thermodynamic database a warning is generated, and if no species are found the value of \code{thermo$species} (which is \code{NULL} at program start-up) is returned. A successful deletion returns the number of species remaining after deletion, or \code{NULL} if no species remain.

}

\seealso{

  \code{\link{basis}} for what must happen before \code{\link{species}}, \code{\link{info}} for searching the thermodynamic database, \code{\link{basis.comp}} (and \code{\link{makeup}}) for generating the coefficients of formation reactions. \code{\link{affinity}} depends on the definitions made using \code{\link{species}}, but \code{\link{subcrt}}, used to calculate standard molal properties, does not.

}

\examples{
  \dontshow{data(thermo)}

  ## add, modify, delete species
  basis('CHNOS')
  species(c('CO2','NH3'))  # aqueous species
  species(c('CO2','NH3'),'gas')  # gases
  # delete the first couple of species
  species(1:2,delete=TRUE)
  # modify the logarithms of activities (actually
  # fugacities) of the remaining species
  species(1:2,c(-2,-5))
  # set the species to aqueous
  species(1:2,'aq')
  # delete all the species
  species(delete=TRUE)

  ## add and delete species
  basis(c('CaO','CO2','H2O','SiO2','MgO','O2'))
  # a species with Fe (faylite) can be loaded,
  # but its initial logarithm of activity is NA to
  # indicate that it falls outside the basis definition
  species(c('dolomite','quartz','calcite','forsterite','fayalite'))
  # changing the elements in the basis definition
  # causes species to be deleted
  basis(c('CO2','H2O','O2'))
  species()  # NULL
  
}

\keyword{misc}
