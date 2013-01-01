\encoding{UTF-8}
\name{EOSregress}
\alias{EOSregress}
\alias{EOSvar}
\alias{EOScalc}
\alias{EOSplot}
\alias{EOSlab}
\alias{EOScoeffs}
\title{Regress Equations-of-State Parameters for Aqueous Species}
\description{
  Functions for fitting experimental volume and heat capacities using regression equations. Possible models include the Helgeson-Kirkham-Flowers (HKF) equations of state and other equations defined using any combination of terms derived from the temperature, pressure and thermodynamic and electrostatic properties of water. 
}

\usage{
  EOSregress(exptdata, var = "", T.max = 9999)
  EOSvar(var, T, P)
  EOScalc(coefficients, T, P)
  EOSplot(exptdata, var = NULL, T.max = 9999, T.plot = NULL, 
    P = NULL, fun.legend = "topleft", coefficients = NULL)
  EOSlab(var, coeff = "")
  EOScoeffs(species, property)
}

\arguments{
  \item{exptdata}{dataframe, experimental data}
  \item{var}{character, name(s) of variables in the regression equations}
  \item{T.max}{numeric, maximum temperature for regression, in degrees Kelvin}
  \item{T}{numeric, temperature in degrees Kelvin}
  \item{P}{numeric, pressure in bars}
  \item{T.plot}{numeric, upper limit of temperature range to plot}
  \item{fun.legend}{character, where to place legend on plot}
  \item{coefficients}{dataframe, coefficients to use to make line on plot}
  \item{coeff}{numeric, value of equation of state parameter for plot legend}
  \item{species}{character, name of aqueous species}
  \item{property}{character, \samp{Cp} or \samp{V}}
}

\details{
  \code{EOSregress} uses \code{\link{lm}} to regress the experimental heat capacity or volume data in \code{exptdata}, which is a data.frame with columns \samp{T} (temperature in degrees Kelvin), \samp{P} (pressure in bars), and \samp{Cp} or \samp{V} (heat capacity in cal/mol.K or volume in cm3/mol). Only data below the temperature of \code{T.max} are included in the regression. The regression formula is specified by a vector of names in \code{var}; these names correspond to variables identified below:

  \tabular{ll}{
    \code{T}           \tab  \eqn{T}{T} (temperature) \cr
    \code{P}           \tab  \eqn{P}{P} (pressure) \cr
    \code{TTheta}      \tab  \eqn{(T-\Theta)}{(T-Theta)} (\eqn{\Theta}{Theta} = 228 K) \cr
    \code{invTTheta}   \tab  \eqn{1/(T-\Theta)}{1/(T-Theta)} \cr
    \code{TTheta2}     \tab  \eqn{(T-\Theta)^2}{(T-Theta)^2} \cr
    \code{invTTheta2}  \tab  \eqn{1/(T-\Theta)^2}{1/(T-Theta)^2} \cr
    \code{invPPsi}     \tab  \eqn{1/(P+\Psi)}{1/(P+Psi)} \cr
    \code{invPPsiTTheta}   \tab  \eqn{1/((P+\Psi)(T-\Theta))}{1/((P+Psi)(T-Theta))} \cr
    \code{V}           \tab  \eqn{V}{V} (volume of water) \cr
    \code{E}           \tab  \eqn{E}{E} (isobaric expansivity of water) \cr
    \code{kT}          \tab  \eqn{\kappa_T}{kT} (isothermal compressibility of water) \cr
    \code{alpha}       \tab  \eqn{\alpha}{alpha} (coefficient of isobaric expansivity of water) \cr
    \code{beta}        \tab  \eqn{\beta}{beta} (coefficients of isothermal compressibility of water) \cr
    \code{X}           \tab  \eqn{X}{X} (Born function \eqn{X}{X}) \cr
    \code{Q}           \tab  \eqn{Q}{Q} (Born function \eqn{Q}{Q}) \cr
    \code{TX}          \tab  \eqn{TX}{TX} (temperature times \eqn{X}{X}) \cr
    \code{drho.dT}     \tab  \eqn{d\rho/dT}{drho/dT} (temperature derivative of density of water) \cr
    \code{V.kT}        \tab  \eqn{V\kappa_T}{V.kT} (volume times isothermal compressibility of water) 
  }

  \code{EOSvar} takes as input \code{var} (one of the names of the variables listed above), and \code{T} (temperature in degrees Kelvin), \code{P} (pressure in bars). It returns the value of the variable at the specified temperature-pressure condition(s). This function is used by \code{EOSregress} to get the values of the variables used in the regression.

  \code{EOScalc} calculates the predicted heat capacities or volumes using coefficients provided by the result of \code{EOSregress}, at the temperatures and pressures specified by \code{T} and \code{P}.

  \code{EOSplot} takes a table of data in \code{exptdata}, runs \code{EOSregress} and \code{EOSpred} and plots the results. The experimental data are plotted as points, and the calculated values as a smooth line. The point symbols are filled circles where the calculated value is within 10\% of the experimental value; open circles otherwise.

  \code{EOSlab} produces labels for the variables listed above that can be used \code{\link{as.expression}}s in plots. The value of \code{coeff} is prefixed (using \code{\link{substitute}}) to the name of the variable.

  \code{EOScoeffs} retrieves coefficients in the Helgeson-Kirkham-Flowers equations from the thermodynamic database (\code{\link{thermo}$obigt}) for the given aqueous \code{species}. If the \code{property} is \samp{Cp}, the resulting dataframe has column names of \samp{(Intercept)}, \samp{invTTheta2} and \samp{TX}, respectively holding the coefficients \eqn{c_1}{c1}, \eqn{c_2}{c2} and \eqn{\omega}{omega} in equation \eqn{Cp^\circ = c_1 + c_2/(T-\Theta)^2 + {\omega}TX}{Cp = c1 + c2/(T-Theta)^2 + omega*TX}. If the \code{property} is \samp{V}, the data frame has column names of \samp{(Intercept)}, \samp{invTTheta} and \samp{Q}, respectively holding the coefficients \eqn{\sigma}{sigma}, \eqn{\xi}{xi} and \eqn{-\omega}{-omega} in \eqn{V^\circ = \sigma + \xi/(T-\Theta) - {\omega}Q}{V = sigma + xi/(T-Theta) - omega*Q}.

  The motivation for writing these functions is to explore alternatives or possible modifications to the revised Helgeson-Kirkham-Flowers equations applied to aqueous nonelectrolytes. As pointed out by Schulte et al., 2001, the functional forms of the equations do not permit retrieving values of the solvation parameter (\eqn{\omega}{omega}) that closely represent the observed trends in both heat capacity and volume at high temperatures (above ca. 200 \eqn{^{\circ}}{°}C).

}

\value{
  For \code{EOSregress}, an object of class \dQuote{lm}. \code{EOSvar} and \code{EOScalc} both return numeric values. \code{EOScoeffs} returns a data frame.
}

\seealso{
  See \code{\link{lm}} for the details of the regression calculations.
}

\references{
  Schulte, M. D., Shock, E. L. and Wood, R. H. (1995) The temperature dependence of the standard-state thermodynamic properties of aqueous nonelectrolytes. \emph{Geochim. Cosmochim. Acta} \bold{65}, 3919--3930. \url{http://dx.doi.org/10.1016/S0016-7037(01)00717-7}
}


\examples{
\dontshow{data(thermo)}
## regress calculated heat capacities and volumes of CH4(aq)
## to test that regressions return known HKF parameters
# calculate the properties of CH4 on a T-P grid
T <- convert(seq(0, 350, 25), "K")
P <- seq(200, 1000, 100)
# convert=FALSE means that temperature has units of K
CH4.prop <- subcrt("CH4", T=T, P=P, grid="T", convert=FALSE)$out[[1]]
# terms in the HKF equations for Cp
Cp.var <- c("invTTheta2", "TX")
# get coefficients in Cp regression
Cp.lm <- EOSregress(CH4.prop[, c("T", "P", "Cp")], Cp.var)
Cp.coeff <- Cp.lm$coefficients
# terms in the HKF equations for V
V.var <- c("invPPsi", "invTTheta", "invPPsiTTheta", "Q")
# get coefficients in V regression
V.lm <- EOSregress(CH4.prop[, c("T", "P", "V")], V.var)
# use same units as HKF: convert from cm3.bar to calories (divide by 41.84)
V.coeff <- convert(V.lm$coefficients, "calories")
## the tests: did we get the HKF parameters that are in the database?
CH4.par <- info(info("CH4"))
# c1 and c2
stopifnot(all.equal(Cp.coeff[1], CH4.par$c1, check.attr=FALSE))
stopifnot(all.equal(Cp.coeff[2], CH4.par$c2, check.attr=FALSE))
# omega (from Cp)
stopifnot(all.equal(Cp.coeff[3], CH4.par$omega, check.attr=FALSE))
# a1, a2, a3 and a4
stopifnot(all.equal(V.coeff[1], CH4.par$a1, check.attr=FALSE))
stopifnot(all.equal(V.coeff[2], CH4.par$a2, check.attr=FALSE))
stopifnot(all.equal(V.coeff[3], CH4.par$a3, check.attr=FALSE))
stopifnot(all.equal(V.coeff[4], CH4.par$a4, check.attr=FALSE))
# omega (from V) - note negative sign
stopifnot(all.equal(-V.coeff[5], CH4.par$omega, check.attr=FALSE))

## regress experimental heat capacities of CH4
## using revised Helgeson-Kirkham-Flowers equations
# read the data from Hnedkovsky and Wood, 1997
f <- system.file("extdata/cpetc/Cp.CH4.HW97.csv", package="CHNOSZ")
d <- read.csv(f)
# have to convert J to cal and MPa to bar
d$Cp <- convert(d$Cp, "cal")
d$P <- convert(d$P, "bar")
# specify the terms in the HKF equations
var <- c("invTTheta2", "TX")
# perform regression, with a temperature limit
EOSlm <- EOSregress(d, var, T.max=600)
# the result is within 10% of the accepted
# values of c1, c2 and omega for CH4(aq)
CH4coeffs <- EOScoeffs("CH4", "Cp")
dcoeffs <- EOSlm$coefficients - CH4coeffs
stopifnot(all(abs(dcoeffs/CH4coeffs) < 0.1))
## make plots comparing the regressions
## here with the accepted EOS parameters of CH4
par(mfrow=c(2,2))
EOSplot(d, T.max=600)
title("Cp of CH4(aq), fit to 600 K")
legend("bottomleft", pch=1, legend="Hnedkovsky and Wood, 1997")
EOSplot(d, coefficients=CH4coeffs)
title("Cp from EOS parameters in database")
EOSplot(d, T.max=600, T.plot=600)
title("Cp fit to 600 K, plot to 600 K")
EOSplot(d, coefficients=CH4coeffs, T.plot=600)
title("Cp from EOS parameters in database")

## model experimental volumes of CH4
## using HKF equation and an exploratory one
f <- system.file("extdata/cpetc/V.CH4.HWM96.csv", package="CHNOSZ")
d <- read.csv(f)
d$P <- convert(d$P, "bar")
# the HKF equation
varHKF <- c("invTTheta", "Q")
# alpha is the expansivity coefficient of water
varal <- c("invTTheta", "alpha")
par(mfrow=c(2,2))
# for both HKF and the expansivity equation
# we'll fit up to a temperature limit
EOSplot(d, varHKF, T.max=663, T.plot=625)
legend("bottomright", pch=1, legend="Hnedkovsky et al., 1996")
title("V of CH4(aq), HKF equation")
EOSplot(d, varal, T.max=663, T.plot=625)
title("V of CH4(aq), expansivity equation")
EOSplot(d, varHKF, T.max=663)
title("V of CH4(aq), HKF equation")
EOSplot(d, varal, T.max=663)
title("V of CH4(aq), expansivity equation")
# note that the volume regression using the HKF gives
# a result for omega (coefficient on Q) that is
# not consistent with the high-T heat capacities
}

\keyword{extra}
