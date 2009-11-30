\name{basis}
\alias{basis}
\title{Basis Species}
\description{
  Define the basis species of a system. Change their physical states or chemical activities or fugacities. Get the stoichiometries of the basis species in selected species of interest.
}

\usage{
  basis(basis = NULL, values = NULL, values2=NULL, delete = FALSE)
}

\arguments{
  \item{basis}{character, names or formulas of basis species, or numeric, indices of species (rownumbers of \code{thermo$obigt}).}
  \item{values}{character or numeric, physical states or names of buffers, or logarithms of activities or fugacities.}
  \item{values2}{character or numeric, physical states, or logarithms of activities or fugacities.}
  \item{delete}{logical, delete the current basis species definition?}
}

\details{

  This is the first function of two (\code{\link{species}} is the second) used to define for the program the system of interest. The basis species represent the possible range of chemical compositions for all the species of interest. Any valid set of basis species used here must meet two conditions: 1) the number of basis species is the same as the number of chemical elements (including charge) in those species and 2) the square matrix representing the elemental stoichiometries of the basis species has a real inverse. Basis species can, but do not always (and not if a charged basis species is present), correspond to the thermodynamic components of a system.

  To create a basis definition, call this function with the names or formulas of the basis species in the first argument. If successful, \code{\link{thermo}$basis} is replaced with the new basis definition. If the second and/or third arguments are present: if either one is of type numeric it refers to values for logarithms of activities (or fugacities) of basis species identified in the first argument. If either of these arguments is character, it refers to the name of a state (if present in \code{thermo$obigt$state}) or to the name of a \code{\link{buffer}} (if present in \code{thermo$buffer$name}).

  To update the logarithms of activities of basis species, provide their formulas or names in the first argument and the values in the second. To change the basis definition, specify the names or formulas of the new basis species in the first argument. When the basis definition is changed, any species of interest that were present are deleted, unless the new basis definition has exactly the same elements as before. In that case, the species are kept and the activities of the new basis species are set so that the chemical potentials of the elements at 25 \eqn{^{\circ}}{degrees }C and 1 bar are unchanged.

  If the value of \code{basis} is one of the keywords in the following table, the corresponding set of basis species is loaded, and their activities set to reference values useful for carrying out the examples. The basis species identified by these keywords are aqueous except for \eqn{\mathrm{H_2O}}{H2O} (liq), \eqn{\mathrm{O_2}}{O2} (gas) and \eqn{\mathrm{Fe_2O_3}}{Fe2O3} (hematite, cr1).

  \tabular{ll}{
    \code{CHNOS} \tab \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, \eqn{\mathrm{NH_3}}{NH3}, \eqn{\mathrm{H_2S}}{H2S}, \eqn{\mathrm{O_2}}{O2} \cr
    \code{CHNOS+} \tab \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, \eqn{\mathrm{NH_3}}{NH3}, \eqn{\mathrm{H_2S}}{H2S}, \eqn{\mathrm{O_2}}{O2}, \eqn{\mathrm{H^+}}{H+} \cr
    \code{CHNOSe} \tab \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, \eqn{\mathrm{NH_3}}{NH3}, \eqn{\mathrm{H_2S}}{H2S}, \eqn{e^-}{e-}, \eqn{\mathrm{H^+}}{H+} \cr
    \code{CHNOPS+} \tab \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, \eqn{\mathrm{NH_3}}{NH3}, \eqn{\mathrm{H_3PO_4}}{H3PO4}, \eqn{\mathrm{H_2S}}{H2S}, \eqn{e^-}{e-}, \eqn{\mathrm{H^+}}{H+} \cr
    \code{MgCHNOPS+} \tab \eqn{\mathrm{Mg^{+2}}}{Mg+2}, \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, \eqn{\mathrm{NH_3}}{NH3}, \eqn{\mathrm{H_3PO_4}}{H3PO4}, \eqn{\mathrm{H_2S}}{H2S}, \eqn{e^-}{e-}, \eqn{\mathrm{H^+}}{H+} \cr
    \code{FeCHNOS} \tab \eqn{\mathrm{Fe_2O_3}}{Fe2O3}, \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, \eqn{\mathrm{NH_3}}{NH3}, \eqn{\mathrm{H_2S}}{H2S}, \eqn{\mathrm{O_2}}{O2} \cr
    \code{FeCHNOS+} \tab \eqn{\mathrm{Fe_2O_3}}{Fe2O3}, \eqn{\mathrm{CO_2}}{CO2}, \eqn{\mathrm{H_2O}}{H2O}, \eqn{\mathrm{NH_3}}{NH3}, \eqn{\mathrm{H_2S}}{H2S}, \eqn{\mathrm{O_2}}{O2}, \eqn{\mathrm{H^+}}{H+} \cr
  }

  Set \code{delete} to TRUE to clear the current basis definition.

}

\value{
  Upon defining or changing the basis definition, the function returns the value of \code{thermo$basis}. 
  With no arguments present, the function invisibly returns the square dataframe of stoichiometries of elements in the basis species (or NULL if \code{thermo$basis} does not exist). This return value is only the stoichiometric matrix contained in \code{thermo$basis}.
  If \code{basis} is numeric, it refers to indices of species (i.e., rownumbers of \code{thermo$obigt}). As long as the basis species were previously defined, the function in this case calculates and returns the basis stoichiometry(s) of the corresponding specie(s).
}

\seealso{
  \code{\link{info}} to query the thermodynamic database in order to find what species are available. \code{\link{species}} is the usual next step after \code{basis}. \code{\link{makeup}} is used by \code{basis} to generate the stoichiometric matrix from chemical formulas.
}

\examples{
  \dontshow{data(thermo)}
  ## define basis species
  basis("O2")  # one element
  basis(c("H2O","O2"))  # two elements
  basis(c("H2O","O2","H+"))

  \dontrun{
    # marked dontrun because they produce intentional errors
    # fewer basis species than elements
    basis(c("H2O","H+"))
    # more basis species than elements
    basis(c("H2O","O2","H2","H+"))
    # non-independent basis species
    basis(c("CO2","H2O","HCl","Cl-","H+"))
  }

  ## specify activities and states
  basis(c("H2O","O2","CO2"),c(-2,-78,-3),c("liq","gas","aq"))
  # change logarithms of activities/fugacities	
  basis(c("H2O","O2"),c(0,-80))	
  # change state of CO2
  basis("CO2","gas")

  ## get the basis stoichiometry of species
  basis("CHNOS")			
  ispecies <- info(c("glutamic acid","phenylalanine"))
  basis(ispecies)
			
  ## changing the basis species
  basis("CHNOPS+")
  # load species whose names are like "ATP"
  species(info("ATP ")[[1]])
  # notice missing elements; try a new basis
  basis("MgCHNOPS+")  # species are deleted
  species(info("ATP ")[[1]])
  # different basis species, but same elements
  basis(c("MgSO4","HCHO","CH4","NO","CO","H3PO4","CS2","OH-"))
  species()  # formation reactions were updated
}

\keyword{misc}
