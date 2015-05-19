\encoding{UTF-8}
\name{util.affinity}
\alias{util.affinity}
\alias{energy}
\alias{energy.args}
\alias{A.ionization}
\title{Functions to Work With Chemical Affinities}
\description{
  The \dQuote{guts} of the affinity calculations: compute affinities or other thermodynamic properties of formation reactions of species on multidimensional arrays with dimensions corresponding to the ranges of intesive variables requested by the user; calculate affinities of ionization reactions of proteins.
}

\usage{
  energy(what, vars, vals, lims, T=get("thermo")$opt$Tr, P="Psat", IS=0, 
    sout=NULL, exceed.Ttr=FALSE, transect = FALSE)
  energy.args(args)
  A.ionization(iprotein, vars, vals, T=get("thermo")$opt$Tr,
    P="Psat", pH=7, transect=FALSE)
}

\arguments{
  \item{what}{character, name of property to calculate}
  \item{vars}{character, names of variables over which to calculate a property}
  \item{vals}{list of numeric, values for each variable}
  \item{lims}{list of numeric, limits of the values for each variable}
  \item{T}{numeric, temperature. Default is to take the temperature from \code{\link{thermo}$opt$Tr}, which corresponds to 25 \eqn{^{\circ}}{°}C}
  \item{P}{numeric, pressure, or character "Psat" (default), which denotes 1 bar or the saturation vapor pressure of \eqn{\mathrm{H_2O}}{H2O} above 100 \eqn{^{\circ}}{°}C (see \code{\link{water}})}
  \item{IS}{numeric, ionic strength; default is 0 mol kg\eqn{^{-1}}{^-1}}
  \item{sout}{list, output from subcrt function}
  \item{exceed.Ttr}{logical, allow subcrt to compute properties for phases beyond their transition temperature?}
  \item{transect}{logical, perform calculations on a transect instead of a grid?}
  \item{args}{list, defines the variables over which to calculate properties}
  \item{iprotein}{numeric, rownumber in \code{thermo$protein}}
  \item{pH}{numeric, pH}
}

\details{
  \code{energy} is the engine for the calculations of chemical affinity. Given \eqn{n} (which can be zero, one, or more) names of basis species and/or \samp{T}, \samp{P}, or \samp{IS} as the \code{var}s, it calculates the property given in \code{what} on an \eqn{n}-dimensional grid  or transect for each of the values (\code{vals}) of the corresponding variable. The limits for each variable given in \code{lims} indicate the minimum and maximum value and, if a third value is supplied, the resolution, or number of points in the given dimension. If \samp{T}, \samp{P}, and/or \samp{IS} are not among the \code{var}s, their constant values can be supplied in \code{T} (in Kelvin), \code{P} (in bar, or \samp{Psat}), and \code{IS} (in mol kg\eqn{^{-1}}{^-1}). \code{sout}, if provided, replaces the call to \code{\link{subcrt}} which can greatly speed up the calculations if this intermediate step is stored by other functions (e.g., \code{\link{transfer}}). \code{exceed.Ttr} is passed to \code{\link{subcrt}} so that the properties of mineral phases beyond their transition temperatures can optionally be calculated.

  The \code{what} argument of \code{energy} is analogous to the \code{property} argument of \code{affinity}.

  \code{energy.args} is used by \code{affinity} to generate the argument list for \code{energy}. \code{energy.args} also has the job of converting \samp{Eh} to \samp{pe} as a function of temperature (see \code{\link{convert}}), and converting \samp{pe} and \samp{pH} to logarithms of activities of the electron and protein, respectively (i.e., negating the values). 


  In CHNOSZ version 0.9, \code{energy} gained a new argument \samp{transect} which is set to TRUE by \code{energy.args} when the length(s) of the variables is(are) greater than three. In this mode of operation, instead of performing the calculations on an \eqn{n}{n}-dimensional grid, the affinities are calculated on an \eqn{n}{n}-dimensional transect through chemical potential (possibly including T and/or P) space. 

  \code{A.ionization} builds a list of values of A/2.303RT of the ionization reactions of proteins that are a function of \code{T}, \code{P} and \code{pH} but are expanded to as many dimensions as defined in \code{vars} in order to be included by the calculations by \code{energy}. These calculations are invoked if proteins are in the species definition, and the basis species contain \samp{H+}.
}

\value{
  For \code{energy}, a list the first element of which is \code{sout} (the results from \code{\link{subcrt}}) and the second element of which is \code{a}, which contains the calculated properties. The latter itself is a list, one element for each species of interest, which have dimensions that are the number of \code{var}iables passed to the function.

  For \code{energy.args}, a list with elements \code{what}, \code{vars}, \code{vals}, \code{lims}, \code{T}, \code{P}, \code{IS} that are appropriate for the corresponding arguments in \code{energy}.
}

\seealso{
  In most cases, \code{\link{affinity}} is used interactively instead of these functions.
}


\examples{\dontshow{data(thermo)}
basis("CHNOS")
species("acetic acid")
eargs <- energy.args(list(O2=c(-90, -60, 5), T=c(0, 100, 5)))
ea <- do.call(energy, eargs)
}

\keyword{utilities}
