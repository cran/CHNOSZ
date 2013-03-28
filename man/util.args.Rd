\name{util.args}
\alias{util.args}
\alias{eos.args}
\alias{TP.args}
\alias{state.args}
\title{Functions for Processing Argument Lists}

\description{
Handle arguments referring to temperature, pressure, states, and equations of state. 
}

\usage{
  eos.args(eos, property = NULL, T = NULL, P = NULL)
  TP.args(T = NULL, P = NULL)
  state.args(state = NULL)
}

\arguments{
  \item{eos}{character, name of equation of state (one of \samp{hkf}, \samp{cgl}, \samp{water}).}
  \item{property}{character, name(s) of thermodynamic properties.}
  \item{T}{numeric, temperature (K).}
  \item{P}{numeric, pressure (bar) (can also be character, \samp{Psat} in \code{TP.args}).}
  \item{state}{character, name(s) of states (e.g., \samp{cr}, \samp{aq}).}
}

\details{
  The \code{*.args} functions are used to normalize user-input arguments, which are case-insensitive. \code{eos.args} returns a list with elements named \code{props}, for all the properties available for the specified equations-of-state, \code{prop} for the lower-case version of \code{property}, and \code{Prop}, for the upper-case (of first letter) version of \code{property}. \code{eos.args} produces an error if one of the \code{property}s is not in the list of available properties. (See \code{\link{water}} and \code{\link{subcrt}} for the available properties for different species.) \code{TP.args} forces \code{T} and \code{P} to equal length. This function also looks for the keyword \samp{Psat} in the value of \code{P} and substitutes calculated values of the saturation vapor pressure (see \code{\link{water}}). \code{state.args} makes its argument lowercase, then transforms \samp{a}, \samp{c}, \samp{g}, and \samp{l} to \samp{aq}, \samp{gas}, \samp{cr}, and \samp{liq}, respectively.
}

\value{
A list is return by \code{eos.args} and \code{TP.args}, and character is returned by \code{state.args}.
}

\examples{\dontshow{data(thermo)}
## argument processing
eos.args("hkf",c("g","H","S","cP","V","kT","e"))
## produces an error because "Q" is not allowed in water.SUPCRT92
\dontrun{
  eos.args("hkf",c("G","H","S","Cp","V","kT","E","Q")) }
  thermo$opt$water <<- "IAPWS"  # needed for p and n in next line
  eos.args("water",c("p","u","cv","psat","rho","n","q","x","y","epsilon"))
  TP.args(c(273.15,373.15))
  TP.args(c(273.15,373.15),"Psat")
  TP.args(c(273.15,373.15),c(100,100,200,200))
  state.args(c("AQ","GAS"))
  state.args(c("a","l","liq"))
}

\keyword{util}
