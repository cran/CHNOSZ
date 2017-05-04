\encoding{UTF-8}
\name{basis}
\alias{basis}
\title{Define Basis Species}
\description{
Define the basis species of a chemical system.
Change their physical states or chemical activities or fugacities.
Get the reaction coefficients of basis species in formation reactions of species.
}

\usage{
  basis(species = NULL, state = NULL, logact = NULL, delete = FALSE)
}

\arguments{
  \item{species}{character, names or formulas of species, or numeric, indices of species}
  \item{state}{character, physical states or names of buffers}
  \item{logact}{numeric, logarithms of activities or fugacities}
  \item{delete}{logical, delete the current basis species definition?}
}

\details{
\code{basis} is used to define the basis species in a thermodynamic system.

The basis species represent the possible range of chemical compositions for all the species of interest.
Any valid set of basis species used here must meet two conditions: 1) the number of basis species is the same as the number of chemical elements (including charge) in those species and 2) the square matrix representing the elemental stoichiometries of the basis species has a real inverse.

To create a basis definition, call \code{basis} with the names or formulas of the basis species in the first argument.
Alternatively, the first argument may consist of numeric values indicating the species indices (rownumbers in \code{\link{thermo}$obigt}), but a mixture of character and numeric values will generate an error.
The special names \samp{pH}, \samp{pe} and \samp{Eh} can be included in the \code{species} argument; they get translated into the names of the proton (\samp{H+}) and electron (\samp{e-}) as appropriate.

The physical states or logarithms of activities of species in the basis definition can be changed by calling \code{basis} with the formulas of species that are in the basis set, or their species indices.
If either of the second or third arguments to \code{basis} is of type character, it refers to the name of a state (if present in \code{thermo$obigt$state}) or to the name of a chemical activity \code{\link{buffer}} (if present in \code{thermo$buffers$name}).
If either of these arguments is numeric it specifies the logarithms of activities (or fugacities for gases) of the basis species.
In case \samp{pH}, \samp{pe} or \samp{Eh} is named, the logarithm of activity of the basis species is converted from these values.
For example, a value of 7 for pH is stored as a logarithm of activity of -7.

Whenever \code{basis} is called with NULL values of both \code{state} and \code{logact}, the new set of species, if they are a valid basis set, completely replaces any existing basis definition.
If this occurs, any existing species definition (created by the \code{species} function) is deleted.
However, \code{\link{swap.basis}} can be used to change the species (the compositions and/or physical states thereof) in the basis set while maintaining the list of species of interest, with the added benefit of equivalence of the chemical potentials of the elements before and after the swap.

Call \code{basis} with \code{delete} set to TRUE to clear the basis definition.
Any current basis definition (before being deleted) is returned by this call or by calling \code{basis} with all default arguments.

If the value of \code{basis} is one of the keywords in the following table, the corresponding set of basis species is loaded, and their activities are given preset values.
This approach is used by many of the examples in the package.
The basis species identified by these keywords are aqueous except for \H2O (liq), \O2 (gas) and \Fe2O3 (hematite, cr1).

  \tabular{ll}{
    \code{CHNOS} \tab \CO2, \H2O, \NH3, \H2S, \O2 \cr
    \code{CHNOS+} \tab \CO2, \H2O, \NH3, \H2S, \O2, \Hplus \cr
    \code{CHNOSe} \tab \CO2, \H2O, \NH3, \H2S, \eminus, \Hplus \cr
    \code{CHNOPS+} \tab \CO2, \H2O, \NH3, \H3PO4, \H2S, \eminus, \Hplus \cr
    \code{MgCHNOPS+} \tab \Mgplus2, \CO2, \H2O, \NH3, \H3PO4, \H2S, \eminus, \Hplus \cr
    \code{FeCHNOS} \tab \Fe2O3, \CO2, \H2O, \NH3, \H2S, \O2 \cr
    \code{FeCHNOS+} \tab \Fe2O3, \CO2, \H2O, \NH3, \H2S, \O2, \Hplus \cr
    \code{QEC4} \tab cysteine, glutamic acid, glutamine, \H2O, \O2 \cr
    \code{QEC} \tab cysteine, glutamic acid, glutamine, \H2O, \O2 \cr
    \code{QEC+} \tab cysteine, glutamic acid, glutamine, \H2O, \O2, \Hplus \cr
  }

The logarithms of activities of amino acids in the \samp{QEC4} basis are -4 (i.e., basis II in Dick, 2016); those in \samp{QEC} and \samp{QEC+} are set to approximate concentrations in human plasma (see Dick, 2017).
}

\value{
Returns the value of \code{thermo$basis} after any modifications; or, if \code{delete} is TRUE, its value before deletion.
}

\section{Side Effects}{
If the new basis definition meets all requirements, it is stored in \code{\link{thermo}$basis}, replacing any previous basis definition.
}

\seealso{
  \code{\link{info}} to query the thermodynamic database in order to find what species are available. 
  \code{\link{makeup}} is used by \code{basis} to generate the stoichiometric matrix from chemical formulas.
  \code{\link{swap.basis}} is used to change the chemical compounds (species formulas) used in the basis definition while keeping the chemical potentials of the elements unaltered.
  \code{\link{species}} for setting up the formation reactions from basis species. 
}

\references{
Dick, J. M. (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}:e2238. \url{https://doi.org/10.7717/peerj.2238}

Dick, J. M. (2017) Chemical composition and the potential for proteomic transformation in cancer, hypoxia, and hyperosmotic stress. \emph{bioRxiv} \url{https://doi.org/10.1101/097667}
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
