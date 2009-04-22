\name{basis}
\alias{basis}
\title{Basis Species}
\description{
  Define the basis species of a system under study. Change their physical states or chemical activities or fugacities. Get the stoichiometries of the basis species in selected species of interest.
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

  This is the first function of two (\code{\link{species}} is the second) used to identify for the program the system of interest. The basis species represent the possible range of compositions for all the species of interest. Any valid set of basis species used here must meet two conditions: 1) the number of basis species is the same as the number of chemical elements (including charge) in those species and 2) the square matrix representing the elemental stoichiometries of the basis species has a real inverse. Basis species can, but do not always (and not if a charged basis species is present), correspond to the thermodynamic components of a system.

  To create a basis definition, call this function with the names or formulas of the basis species in the first argument. If successful, \code{\link{thermo}$basis} is replaced with the new basis definition. If the second or third arguments are present: if either one is of type numeric it refers to values for logarithms of activities (or fugacities) of basis species identified in the first argument. If either of these arguments is character, it refers to the name of a state (if present in \code{thermo$obigt$state}) or to the name of a \code{\link{buffer}} (if present in \code{thermo$buffer$name}).

  To update the logarithms of activities of basis species, provide their formulas or names in the first argument with the values in the second. To change the basis definition, specify the names or formulas of the new basis species in the first argument. When the basis definition is changed, any species of interest that were present are deleted, unless the new basis definition has exactly the same elements as before. In that case, the species are kept and the activities of the new basis species are set so that the chemical potentials of the elements at 25 \eqn{^{\circ}}{degrees }C and 1 bar are unchanged.

  If the value of \code{basis} is one of the keywords in the following table, the corresponding set of basis species is loaded, and their activities set to reference values useful for carrying out the examples. The basis species identified by these keywords are aqueous except for H2O (liq), O2 (gas) and Fe2O3 (hematite, cr1).

  \tabular{ll}{
    \code{CHNOS} \tab CO2, H2O, NH3, H2S, O2 \cr
    \code{CHNOS+} \tab CO2, H2O, NH3, H2S, O2, H+ \cr
    \code{CHNOSe} \tab CO2, H2O, NH3, H2S, e-, H+ \cr
    \code{CHNOPS+} \tab CO2, H2O, NH3, H3PO4, H2S, e-, H+ \cr
    \code{MgCHNOPS+} \tab Mg+2, CO2, H2O, NH3, H3PO4, H2S, e-, H+ \cr
    \code{FeCHNOS} \tab Fe2O3, CO2, H2O, NH3, H2S, O2 \cr
    \code{FeCHNOS+} \tab Fe2O3, CO2, H2O, NH3, H2S, O2, H+ \cr
  }

  If \code{basis} is numeric, calculate and return the basis stoichiometry(s) of the corresponding specie(s) identified by their rownumbers in \code{thermo$obigt} (the set of basis species must have been previously defined).

  Set \code{delete} to TRUE to clear the current basis definition.

}

\value{
  With no arguments present, the dataframe of stoichiometries of elements in the basis species (or NULL if \code{thermo$basis} does not exist). Note that the \code{\link{invisible}} return value is different from (is only the stoichiometric part of) the screen output of this function, which prints the entire dataframe in \code{thermo$basis}.
}

\seealso{
  \code{\link{species}}, \code{\link{info}}
}

\examples{
  \dontshow{data(thermo)}

  ## define basis species
  basis('O2')  # one element
  basis(c('H2O','O2'))  # two elements
  basis(c('H2O','O2','H+'))

  \dontrun{
    # marked dontrun because they produce errors
    # fewer basis species than elements
    basis(c('H2O','H+'))
    # more basis species than elements
    basis(c('H2O','O2','H2','H+'))
    # non-independent basis species
    basis(c('CO2','H2O','HCl','Cl-','H+'))
  }

  ## specify activities and states
  basis(c('H2O','O2','CO2'),c(-2,-78,-3),c('liq','gas','aq'))
  # change logarithms of activities/fugacities	
  basis(c('H2O','O2'),c(0,-80))	
  # change state of CO2
  basis('CO2','gas')

  ## get the basis stoichiometry of a species
  basis('CHNOS')			
  ispecies <- info(c('glutamic acid','phenylalanine'))
  basis(ispecies)
			
  ## changing the basis species
  basis('CHNOPS+')
  # load species whose names are like 'ATP'
  species(info('ATP ')[[1]])
  # notice missing elements; try a new basis
  basis('MgCHNOPS+')  # species are deleted
  species(info('ATP ')[[1]])
  # different basis species, but same elements
  basis(c('MgSO4','HCHO','CH4','NO','CO','H3PO4','CS2','OH-'))
  species()  # formation reactions were updated
  
  
}

\keyword{misc}
