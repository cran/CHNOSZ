\name{util.misc}
\alias{util.misc}
\alias{dPdTtr}
\alias{Ttr}
\alias{nonideal}
\alias{which.balance}
\alias{unitize}
\title{Functions for Miscellaneous Tasks}
\description{
  Calculate \eqn{dP/dT}{dP/dT} and temperature of phase transitions, calculate heat capacities of unfolded proteins using an equation from the literature, calculate non-ideal contributions to apparent standard molal properties, identify a conserved basis species, scale logarithms of activity to a desired total activity, calculate Gibbs energy of transformation of a system.
}

\usage{
  dPdTtr(x)
  Ttr(x, P = 1, dPdT = NULL)
  nonideal(species, proptable, IS, T)
  which.balance(species)
  unitize(logact = NULL, length = NULL, logact.tot = 0)
}

\arguments{
  \item{x}{numeric index of a mineral phase (\code{dPdTtr}, \code{Ttr})}
  \item{P}{numeric, pressure (bar)}
  \item{dPdT}{numeric, values of (\eqn{dP/dT}{dP/dT}) of phase transitions (\code{Ttr})}
  \item{species}{Names or indices of species for which to calculate nonideal properties (\code{nonideal}), or dataframe, species definition such as that in \code{\link{thermo}$species} (\code{which.balance})}
  \item{proptable}{list of dataframes of species properties}
  \item{T}{numeric, temperature (K) (\code{lines.water}, \code{nonideal})}
  \item{IS}{numeric, ionic strength(s) used in nonideal calculations, mol kg\eqn{^{-1}}{^-1}}
  \item{logact}{numeric, logarithms of activity}
  \item{length}{numeric, numbers of residues}
  \item{logact.tot}{numeric, logarithm of total activity}
}

\details{

  

  \code{dPdTtr} returns values of \eqn{(dP/dT)_{Ttr}}{(dP/dT)Ttr}, where \eqn{Ttr}{Ttr} represents the transition temperature, of the phase transition at the high-\eqn{T}{T} stability limit of the \code{x}th species in \code{thermo$obigt} (no checking is done to verify that the species represents in fact one phase of a mineral with phase transitions). \code{dPdTtr} takes account of the Clapeyron equation, \eqn{(dP/dT)_{Ttr}}{(dP/dT)Ttr}=\eqn{{\Delta}S/{\Delta}V}{DS/DV}, where \eqn{{\Delta}S}{DS} and \eqn{{\Delta}V}{DV} represent the changes in entropy and volume of phase transition, and are calculated using \code{subcrt} at Ttr from the standard molal entropies and volumes of the two phases involved. Using values of \code{dPdT} calculated using \code{dPdTtr} or supplied in the arguments, \code{Ttr} returns as a function of \code{P} values of the upper transition temperature of the mineral phase represented by the \code{x}th species.

  \code{nonideal} takes a list of dataframes (in \code{proptable}) containing the standard molal properties of the identified \code{species}. For those species whose charge (determined by the number of Z in their \code{\link{makeup}}) is not equal to zero, the values of \code{IS} are combined with Alberty's (2003) equation 3.6-1 (Debye-Huckel equation) and its derivatives, to calculate apparent molal properties at the specified ionic strength(s) and temperature(s). The lengths of \code{IS} and \code{T} supplied in the arguments should be equal to the number of rows of each dataframe in \code{proptable}, or one to use single values throughout. The apparent molal properties that can be calculated include \code{G}, \code{H}, \code{S} and \code{Cp}; any columns in the dataframes of \code{proptable} with other names are left untouched. A column named \code{loggam} (logarithm of gamma, the activity coefficient) is appended to the output dataframe of species properties.

  \code{which.balance} returns, in order, which column(s) of \code{species} all have non-zero values. It is used by \code{\link{diagram}} and \code{\link{transfer}} to determine a conservant (i.e. basis species that are conserved in transformation reactions) if none is supplied by the user.

  \code{spearman} calculates Spearman's rank correlation coefficient for \code{a} and \code{b}.

  \code{unitize} scales the logarithms of activities given in \code{logact} so that the logarithm of total activity of residues is equal to zero (i.e. total activity of residues is one), or to some other value set in \code{logact.tot}.  \code{length} indicates the number of residues in each species. If \code{logact} is NULL, the function takes the logarithms of activities from the current species definition. If any of those species are proteins, the function gets their lengths using \code{protein.length}.
}

\value{
  Numeric returns are made by \code{dPdTtr}, \code{Ttr}, \code{spearman}, \code{mod.obigt} Functions with no (or unspecified) returns are \code{thermo.plot.new}, \code{thermo.postscript}, \code{label.plot} and \code{water.lines}.
}


\seealso{
  For some of the functions on which these utilities depend or were modeled, see \code{\link{paste}}, \code{\link{substr}}, \code{\link{tolower}}, \code{\link{par}} and \code{\link{text}}.
}

\examples{\dontshow{data(thermo)}
## properties of phase transitions
si <- info("enstatite")
# (dP/dT) of transitions
dPdTtr(si)  # first transition
dPdTtr(si+1) # second transition
# temperature of transitions (Ttr) as a function of P
Ttr(si,P=c(1,10,100,1000))
Ttr(si,P=c(1,10,100,1000))

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
}

\keyword{util}
