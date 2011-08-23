\name{util.misc}
\alias{util.misc}
\alias{protein.length}
\alias{MP90.cp}
\alias{dPdTtr}
\alias{Ttr}
\alias{describe}
\alias{basis.comp}
\alias{nonideal}
\alias{which.balance}
\alias{unitize}
\alias{caller.name}
\title{Functions for Miscellaneous Tasks}
\description{
  Calculate \eqn{dP/dT}{dP/dT} and temperature of phase transitions, calculate heat capacities of unfolded proteins using an equation from the literature, calculate non-ideal contributions to apparent standard molal properties, identify a conserved basis species, execute multicore calculations, scale logarithms of activity to a desired total activity, get name of calling function.
}

\usage{
  protein.length(protein)
  MP90.cp(T, protein)
  dPdTtr(x)
  Ttr(x, P = 1, dPdT = NULL)
  describe(x = NULL, T = NULL, P = NULL, use.name = FALSE, 
    as.reaction = NULL, digits = 1)
  basis.comp(basis)
  nonideal(species, proptable, IS, T)
  which.balance(species)
  unitize(logact = NULL, length = NULL, logact.tot = 0)
  caller.name(n)
}

\arguments{
  \item{protein}{character, name of protein species; numeric, species index of protein.}
  \item{T}{numeric, temperature (K) (\code{lines.water}, \code{describe}, \code{MP90.cp}, \code{nonideal}).}
  \item{x}{numeric index of a mineral phase (\code{dPdTtr}, \code{Ttr}), or dataframe, definition of basis species or reaction (\code{describe}).}
  \item{P}{numeric, pressure (bar).}
  \item{dPdT}{numeric, values of (\eqn{dP/dT}{dP/dT}) of phase transitions (\code{Ttr}).}
  \item{use.name}{logical, write names instead of formulas? (\code{describe}).}
  \item{as.reaction}{logical, interpret input as a reaction?}
  \item{digits}{numeric, how many digits to round logarithms of activities.}
  \item{basis}{numeric or character, species number or formula.}
  \item{species}{Names or indices of species for which to calculate nonideal properties (\code{nonideal}), or dataframe, species definition such as that in \code{\link{thermo}$species} (\code{which.balance}).}
  \item{proptable}{list of dataframes of species properties.}
  \item{IS}{numeric, ionic strength(s) used in nonideal calculations, mol kg\eqn{^{-1}}{^-1}.}
  \item{logact}{numeric, logarithms of activity.}
  \item{length}{numeric, numbers of residues.}
  \item{logact.tot}{numeric, logarithm of total activity.}
  \item{n}{numeric, number of frame to go up.}
}

\details{

  The argument of \code{protein.length}, if it is character, refers to the name of protein(s) (e.g., \samp{LYSC_CHICK}) for which to calculate the length (number of amino acid residues). If the argument is numeric, it refers to the species index of a protein (value in \code{thermo$species$ispecies}). For a positive numeric argument to work, the protein information must have been previously loaded into the species list (using \code{\link{info}}). If the numeric value is negative, it refers to the rownumber of the protein in \code{thermo$protein}.
  
  \code{MP90.cp} takes \code{T} (one or more temperatures in \eqn{^{\circ}}{degrees }C) and \code{protein} (name of protein) and returns the heat capacity of the unfolded protein using values of heat capacities of the residues taken from Makhatadze and Privalov, 1990. Those authors provided values of heat capacity at six points between 5 and 125 \eqn{^{\circ}}{degrees }C; this function interpolates (using \code{\link{splinefun}}) values at other temperatures.

  \code{dPdTtr} returns values of \eqn{(dP/dT)_{Ttr}}{(dP/dT)Ttr}, where \eqn{Ttr}{Ttr} represents the transition temperature, of the phase transition at the high-\eqn{T}{T} stability limit of the \code{x}th species in \code{thermo$obigt} (no checking is done to verify that the species represents in fact one phase of a mineral with phase transitions). \code{dPdTtr} takes account of the Clapeyron equation, \eqn{(dP/dT)_{Ttr}}{(dP/dT)Ttr}=\eqn{{\Delta}S/{\Delta}V}{DS/DV}, where \eqn{{\Delta}S}{DS} and \eqn{{\Delta}V}{DV} represent the changes in entropy and volume of phase transition, and are calculated using \code{subcrt} at Ttr from the standard molal entropies and volumes of the two phases involved. Using values of \code{dPdT} calculated using \code{dPdTtr} or supplied in the arguments, \code{Ttr} returns as a function of \code{P} values of the upper transition temperature of the mineral phase represented by the \code{x}th species.

  \code{describe} generates a textual representation of the temperature, pressure, and logarithms of activities of the basis species, given in the arguments by \code{x} (i.e. the dataframe in \code{thermo$basis}) and \code{T} and \code{P} (given in Kelvin and bar and converted by the function to those specified by \code{\link{nuts}}). The \code{digits} argument tells to what decimal place the logarithms of activities should be rounded. If any of the supplied arguments is NULL its specification is not printed in the output; T and P, if present, are prepended to the basis summary. If \code{x} instead is a dataframe representing a chemical reaction (as output by \code{\link{subcrt}} and identified by having a column named \code{coeff}), the function returns a textual summary of that reaction (i.e., showing reactants on the left, an equal sign, and products on the right; reactants and products are preceded by their reaction coefficient unless it is \eqn{1}). However, if only reactants or products can be found, or \code{as.reaction} is set to FALSE, the names or formulas of the species are printed with their coefficients and interceding plus or minus signs, as appropriate. Whether the names or formulas are printed is controlled by \code{use.name} (\code{FALSE} by default), a logical vector the length of which corresponds to the number of rows in \code{x} (but is expanded to the right length if needed). 

  \code{basis.comp} calculates the composition of the given 

  \code{nonideal} takes a list of dataframes (in \code{proptable}) containing the standard molal properties of the identified \code{species}. For those species whose charge (determined by the number of Z in their \code{\link{makeup}}) is not equal to zero, the values of \code{IS} are combined with Alberty's (2003) equation 3.6-1 (Debye-Huckel equation) and its derivatives, to calculate apparent molal properties at the specified ionic strength(s) and temperature(s). The lengths of \code{IS} and \code{T} supplied in the arguments should be equal to the number of rows of each dataframe in \code{proptable}, or one to use single values throughout. The apparent molal properties that can be calculated include \code{G}, \code{H}, \code{S} and \code{Cp}; any columns in the dataframes of \code{proptable} with other names are left untouched. A column named \code{loggam} (logarithm of gamma, the activity coefficient) is appended to the output dataframe of species properties.

  \code{which.balance} returns, in order, which column(s) of \code{species} all have non-zero values. It is used by \code{\link{diagram}} and \code{\link{transfer}} to determine a conservant (i.e. basis species that are conserved in transformation reactions) if none is supplied by the user.

  \code{spearman} calculates Spearman's rank correlation coefficient for \code{a} and \code{b}.

  \code{unitize} scales the logarithms of activities given in \code{logact} so that the logarithm of total activity of residues is equal to zero (i.e. total activity of residues is one), or to some other value set in \code{logact.tot}.  \code{length} indicates the number of residues in each species. If \code{logact} is NULL, the function takes the logarithms of activities from the current species definition. If any of those species are proteins, the function gets their lengths using \code{protein.length}.

  \code{caller.name} returns the name of the calling function \code{n} frames up (i.e., for \code{n} equal to 2, the caller of the function that calls this one). If called interactively, returns character().
}

\value{

  Numeric returns are made by, \code{protein.length}, \code{dPdTtr}, \code{Ttr}, \code{MP90.cp}, \code{spearman}, \code{mod.obigt} Functions with no (or unspecified) returns are \code{thermo.plot.new}, \code{thermo.postscript}, \code{label.plot} and \code{water.lines}.

}


\seealso{
  For some of the functions on which these utilities depend or were modeled, see  \code{\link{paste}}, \code{\link{substr}}, \code{\link{tolower}}, \code{\link{par}} and \code{\link{text}}.
}

\examples{
  \dontshow{data(thermo)}
  ## calculate protein length
  protein.length("LYSC_CHICK")
  # another way to do it
  basis("CHNOS")
  species("LYSC_CHICK")
  protein.length(species()$ispecies)
  # another way to do it
  ip <- protein("LYSC","CHICK")
  protein.length(-ip)

  ## heat capacity as a function of temperature
  ## (Makhatadze & Privalov, 1990) units: J mol-1
  MP90.cp(c(5,25,50,75,100,125),"LYSC_CHICK")

  ## properties of phase transitions
  si <- info("enstatite")
  # (dP/dT) of transitions
  dPdTtr(si)  # first transition
  dPdTtr(si+1) # second transition
  # temperature of transitions (Ttr) as a function of P
  Ttr(si,P=c(1,10,100,1000))
  Ttr(si,P=c(1,10,100,1000))
  
  ## the basis stoichiometry of a made-up species
  # warns because P isn't in our basis
  basis("CHNOS")
  basis.comp("SPONCH")

  ## describing the basis species
  basis("CHNOSe")
  describe(thermo$basis)
  describe(thermo$basis,T=NULL,P=NULL)

  ## scale logarithms of activity
  # suppose we have two proteins whose lengths are 100 and 
  # 200; what are the logarithms of activity of the proteins 
  # that are equal to each other and that give a total 
  # activity of residues equal to unity?
  logact <- c(-3,-3)  # could be any two equal numbers
  length <- c(100,200)
  logact.tot <- 0
  loga <- unitize(logact,length,logact.tot)
  # the proteins have equal activity
  stopifnot(identical(loga[1],loga[2]))
  # the sum of activity of the residues is unity
  stopifnot(isTRUE(all.equal(sum(10^loga * length),1)))
  ## now, what if the activity of protein 2 is ten
  ## times that of protein 1?
  logact <- c(-3,-2)
  loga <- unitize(logact,length,logact.tot)
  # the proteins have unequal activity
  stopifnot(isTRUE(all.equal(loga[2]-loga[1],1)))
  # but the activities of residues still add up to one
  stopifnot(isTRUE(all.equal(sum(10^loga * length),1)))
}

\references{

  Alberty, R. A. (2003) \emph{Thermodynamics of Biochemical Reactions}, John Wiley & Sons, Hoboken, New Jersey, 397 p. \url{http://www.worldcat.org/oclc/51242181}

  Makhatadze, G. I. and Privalov, P. L. (1990) Heat capacity of proteins. 1. Partial molar heat capacity of individual amino acid residues in aqueous solution: Hydration effect \emph{J. Mol. Biol.} \bold{213}, 375--384. \url{http://dx.doi.org/10.1016/S0022-2836(05)80197-4}

 
}

\keyword{util}
