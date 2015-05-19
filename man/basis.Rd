\name{basis}
\alias{basis}
\alias{put.basis}
\alias{mod.basis}
\alias{preset.basis}
\alias{preset.logact}
\title{Define Basis Species}
\description{
  Define the basis species of a chemical system. Change their physical states or chemical activities or fugacities. Get the stoichiometries of the basis species in selected species of interest.
}

\usage{
  basis(species = NULL, state = NULL, logact = NULL, delete = FALSE)
  put.basis(ispecies, logact)
  mod.basis(species, state, logact)
  preset.basis(key = NULL)
  preset.logact(species)
}

\arguments{
  \item{species}{character, names or formulas of species, or numeric, indices of species}
  \item{state}{character, physical states or names of buffers}
  \item{logact}{numeric, logarithms of activities or fugacities}
  \item{delete}{logical, delete the current basis species definition?}
  \item{ispecies}{numeric, indices of species (rownumbers of \code{thermo$obigt})}
  \item{key}{character, a keyword identifying a preset basis definition}
}

\details{
  \code{basis} is used to define the basis species in a system of interest, and in many workflows is followed by calls to \code{\link{species}}, \code{\link{affinity}} and \code{\link{diagram}} for making equilibrium chemical activity diagrams. The other functions documented here are supporting functions for \code{basis} and generally are not intended to be called by the user.

  The basis species represent the possible range of chemical compositions for all the species of interest. Any valid set of basis species used here must meet two conditions: 1) the number of basis species is the same as the number of chemical elements (including charge) in those species and 2) the square matrix representing the elemental stoichiometries of the basis species has a real inverse. Basis species might, but do not always (and not if a charged basis species is present), correspond to the thermodynamic components of a system. 

  To create a basis definition, call \code{basis} with the names or formulas of the basis species in the first argument. Alternatively, the first argument may consist of numeric values indicating the species indices (rownumbers in \code{\link{thermo}$obigt}), but a mixture of character and numeric values will generate an error. The special names \samp{pH}, \samp{pe} and \samp{Eh} can be included in the \code{species} argument; they get translated into the names of the proton (\samp{H+}) and electron (\samp{e-}) as appropriate.

   If the new basis definition meets all requirements, it is stored in \code{\link{thermo}$basis}, replacing any previous basis definition; \code{put.basis} does the actual storing of the basis definition.

  The physical states or logarithms of activities of species in the basis definition can be changed directly using \code{mod.basis} but usually more conveniently by calling \code{basis} with the formulas of species that are in the basis set, or their species indices. If either of the second or third arguments to \code{basis} is of type character, it refers to the name of a state (if present in \code{thermo$obigt$state}) or to the name of a chemical activity \code{\link{buffer}} (if present in \code{thermo$buffers$name}). If either of these arguments is numeric it specifies the logarithms of activities (or fugacities for gases) of the basis species. In case \samp{pH}, \samp{pe} or \samp{Eh} are named, the logarithm of activity of the basis species is converted from these values. (For example, a value of 7 for pH is stored as a logarithm of activity of -7).

  Whenever \code{basis} is called with NULL values of both \code{state} and \code{logact}, the new set of species, if they are a valid basis set, completely replaces any existing basis definition. If this occurs, any existing species definition (created by the \code{species} function) is deleted. However, \code{\link{swap.basis}} can be used to change the species (the compositions and/or physical states thereof) in the basis set while maintaining the list of species of interest, with the added benefit of equivalence of the chemical potentials of the elements before and after the swap.

  Call \code{basis} with \code{delete} set to TRUE to clear the basis definition. Any current basis definition (before being deleted) is returned by this call or by calling \code{basis} with all default arguments.

  If the value of \code{basis} is one of the keywords in the following table, the corresponding set of basis species is loaded (defined in \code{preset.basis}), and their activities set to reference values (defined in \code{preset.logact}). This approach is used by many of the examples in the package. The basis species identified by these keywords are aqueous except for \eqn{\mathrm{H_2O}}{H2O} (liq), \eqn{\mathrm{O_2}}{O2} (gas) and \eqn{\mathrm{Fe_2O_3}}{Fe2O3} (hematite, cr1).

  \tabular{ll}{
    \code{CHNOS} \tab \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, \eqn{\mathrm{NH_3}}{NH3}, 
                      \eqn{\mathrm{H_2S}}{H2S}, \eqn{\mathrm{O_2}}{O2} \cr
    \code{CHNOS+} \tab \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, \eqn{\mathrm{NH_3}}{NH3}, 
                       \eqn{\mathrm{H_2S}}{H2S}, \eqn{\mathrm{O_2}}{O2}, \eqn{\mathrm{H^+}}{H+} \cr
    \code{CHNOSe} \tab \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, \eqn{\mathrm{NH_3}}{NH3}, 
                       \eqn{\mathrm{H_2S}}{H2S}, \eqn{e^-}{e-}, \eqn{\mathrm{H^+}}{H+} \cr
    \code{CHNOPS+} \tab \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, \eqn{\mathrm{NH_3}}{NH3}, 
                        \eqn{\mathrm{H_3PO_4}}{H3PO4}, \eqn{\mathrm{H_2S}}{H2S}, \eqn{e^-}{e-}, \eqn{\mathrm{H^+}}{H+} \cr
    \code{MgCHNOPS+} \tab \eqn{\mathrm{Mg^{+2}}}{Mg+2}, \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, 
                          \eqn{\mathrm{NH_3}}{NH3}, \eqn{\mathrm{H_3PO_4}}{H3PO4}, \eqn{\mathrm{H_2S}}{H2S}, 
                          \eqn{e^-}{e-}, \eqn{\mathrm{H^+}}{H+} \cr
    \code{FeCHNOS} \tab \eqn{\mathrm{Fe_2O_3}}{Fe2O3}, \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, 
                        \eqn{\mathrm{NH_3}}{NH3}, \eqn{\mathrm{H_2S}}{H2S}, \eqn{\mathrm{O_2}}{O2} \cr
    \code{FeCHNOS+} \tab \eqn{\mathrm{Fe_2O_3}}{Fe2O3}, \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, 
                         \eqn{\mathrm{NH_3}}{NH3}, \eqn{\mathrm{H_2S}}{H2S}, \eqn{\mathrm{O_2}}{O2}, \eqn{\mathrm{H^+}}{H+} \cr
  }

}

\value{
  \code{basis} returns the value of \code{thermo$basis} after any modifications; or, if \code{delete} is TRUE, its value before deletion.


}

\seealso{
  \code{\link{info}} to query the thermodynamic database in order to find what species are available. 
  \code{\link{makeup}} is used by \code{basis} to generate the stoichiometric matrix from chemical formulas.
  \code{\link{species}} is the usual next step after \code{basis}. 
  \code{\link{swap.basis}} is used to change the chemical compounds (species formulas) used in the basis definition while keeping the chemical potentials of the elements unaltered.
}

\examples{
\dontshow{data(thermo)}
## define basis species
# one, two and three element examples
basis("O2")
basis(c("H2O", "O2"))
basis(c("H2O", "O2", "H+"))
## clear the basis species
basis(delete=TRUE)

\dontrun{
## marked dontrun because they produce errors
# fewer species than elements
basis(c("H2O", "H+"))
# more species than elements
basis(c("H2O", "O2", "H2", "H+"))
# non-independent species
basis(c("CO2", "H2O", "HCl", "Cl-", "H+"))}

## specify activities and states
basis(c("H2O", "O2", "CO2"), c(-2, -78, -3), c("liq", "aq", "aq"))
# change logarithms of activities/fugacities	
basis(c("H2O", "O2"), c(0, -80))	
# change state of CO2
basis("CO2", "gas")
}

\concept{Primary thermodynamic modeling}
