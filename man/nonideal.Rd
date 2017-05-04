\encoding{UTF-8}
\name{nonideal}
\alias{nonideal}
\title{Activity coefficients of aqueous species}
\description{
Calculate activity coefficients and non-ideal contributions to apparent standard molal properties of aqueous species.
}

\usage{
  nonideal(species, proptable, IS, T)
}

\arguments{
  \item{species}{names or indices of species for which to calculate nonideal properties}
  \item{proptable}{list of dataframes of species properties}
  \item{IS}{numeric, ionic strength(s) used in nonideal calculations, mol kg\eqn{^{-1}}{^-1}}
  \item{T}{numeric, temperature (K) (\code{lines.water}, \code{nonideal})}
}

\details{
\code{nonideal} takes a list of dataframes (in \code{proptable}) containing the standard molal properties of the identified \code{species}.
The function bypasses (leaves unchanged) properties of all species whose charge (determined by the number of Z in their \code{\link{makeup}}) is equal to zero.
The proton (H+) and electron (e-) are also bypassed by default; to apply the calculations to H+ and/or e-, change \code{\link{thermo}$opt$ideal.H} or \code{ideal.e} to FALSE.
The values of \code{IS} are combined with Alberty's (2003) equation 3.6-1 (extended Debye-Hückel equation) and its derivatives, to calculate apparent molal properties at the specified ionic strength(s) and temperature(s).
The lengths of \code{IS} and \code{T} supplied in the arguments should be equal to the number of rows of each dataframe in \code{proptable}, or one to use single values throughout.
The apparent molal properties that can be calculated include \code{G}, \code{H}, \code{S} and \code{Cp}; any columns in the dataframes of \code{proptable} with other names are left untouched.
A column named \code{loggam} (logarithm of gamma, the activity coefficient) is appended to the output dataframe of species properties.
}

\section{Warning}{
The logarithms of activity coefficients (\code{loggam}) returned by \code{nonideal} use the natural logarithm (cf. Alberty, 2003 Eq. 3.6-1).
To maintain consistency with the conventions used elsewhere in the package (i.e. for logarithms of equilibrium constants and of chemical activities), the values of \code{loggam} returned by \code{\link{subcrt}} are expressed using the common (base 10) logarithm.
Note that the first example below uses \code{loggam} returned by \code{subcrt}, therefore requiring a base of 10 for calculating gamma.
}

\examples{\dontshow{data(thermo)}
### Examples following Alberty, 2003 

## p. 273-276: activity coefficient (gamma)
## as a function of ionic strength and temperature
T <- c(0, 25, 40)
col <- c("blue", "black", "red")
IS <- seq(0, 0.25, 0.0025)
thermo.plot.new(xlim=range(IS), ylim=c(0, 1), xlab=axis.label("IS"), ylab="gamma")
for(j in 1:3) {
  s <- subcrt(c("H2PO4-", "HADP-2", "HATP-3", "ATP-4"), IS=IS, grid="IS", T=T[j])
  for(i in 1:4) lines(IS, 10^s$out[[i]]$loggam, col=col[j])
}
text(0.125, 0.8, "Z = -1")
text(0.1, 0.42, "Z = -2")
text(0.075, 0.18, "Z = -3")
text(0.05, 0.08, "Z = -4")
title(main=paste("activity coefficient (gamma) of -1,-2,-3,-4",
  "charged species at 0, 25, 40 deg C, after Alberty, 2003",
  sep="\n"), cex.main=0.95)
legend("topright", lty=c(NA, 1, 1, 1), col=c(NA, "blue", "black", "red"),
  legend=c(as.expression(axis.label("T")), 0, 25, 40))

## p. 16 Table 1.3: apparent pKa of acetic acid with
## changing ionic strength
# we set this option to FALSE so that nonideal() will calculate activity
# coefficients for the proton (makes for better replication of the values
# in Alberty's book)
thermo$opt$ideal.H <<- FALSE
subcrt(c("acetic acid", "acetate", "H+"), c(-1, 1, 1),
  IS=c(0, 0.1, 0.25), T=25, property="logK")
# note that *apparent* values equal *standard* values at IS=0
# reset option to default
thermo$opt$ideal.H <<- TRUE

## p. 95: basis and elemental stoichiometries of species 
# (this example doesn't use activity coefficients)
basis(c("ATP-4", "H+", "H2O", "HPO4-2", "O2", "NH3"))
# cf Eq. 5.1-33: basis composition
species(c("ATP-4", "H+", "H2O", "HPO4-2", "ADP-3", "HATP-3", "HADP-2",
  "H2PO4-"))

### A different example

# speciation of phosphate as a function of ionic strength
opar <- par(mfrow=c(2, 1))
basis("CHNOPS+")
Ts <- c(25, 100)
species(c("PO4-3", "HPO4-2", "H2PO4-"))
for(T in Ts) {
  a <- affinity(IS=c(0, 0.14), T=T)
  e <- equilibrate(a)
  if(T==25) diagram(e, ylim=c(-3.0, -2.6), legend.x=NULL)
  else d <- diagram(e, ylim=c(-3.0, -2.6), add=TRUE, col="red")
}
title(main="Non-ideality model for phosphate species")
dp <- describe.property(c("pH", "T", "T"), c(7, Ts))
legend("topright", lty=c(NA, 1, 1), col=c(NA, "black", "red"), legend=dp)
text(0.07, -2.76, expr.species("HPO4-2"))
text(0.07, -2.90, expr.species("H2PO4-"))
#
# phosphate predominance f(IS,pH)
a <- affinity(IS=c(0, 0.14), pH=c(6, 13), T=Ts[1])
d <- diagram(a, fill=NULL)
a <- affinity(IS=c(0, 0.14), pH=c(6, 13), T=Ts[2])
d <- diagram(a, add=TRUE, names=NULL, col="red")
par(opar)
}

\references{
Alberty, R. A. (2003) \emph{Thermodynamics of Biochemical Reactions}, John Wiley & Sons, Hoboken, New Jersey, 397 p. \url{http://www.worldcat.org/oclc/51242181}
}

\concept{Secondary thermodynamic modeling}
