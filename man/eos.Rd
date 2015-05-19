\encoding{UTF-8}
\name{eos}
\alias{hkf}
\alias{cgl}
\alias{gfun}
\title{Equations of State}
\description{
  Calculate thermodynamic properties using the revised Helgeson-Kirkham-Flowers (HKF) equations of state for aqueous species, or using a generic heat capacity equation for crystalline, gas, and liquid species.
}

\usage{
  cgl(property = NULL, T = 298.15, P = 1, ghs = NULL, eos = NULL)
  hkf(property = NULL, T = 298.15, P = 1, ghs = NULL, eos = NULL,
    contrib = c("n","s","o"), H2O.PT = NULL, H2O.PrTr = NULL, 
    domega = TRUE)
  gfun(rhohat, Tc, P, alpha, daldT, beta)
}

\arguments{
  \item{property}{character, name(s) of properties to calculate}
  \item{T}{numeric, temperature(s) at which to calculate properties (K)}
  \item{P}{numeric, pressure(s) at which to calculate properties (bar)}
  \item{ghs}{dataframe, values of the standard molal Gibbs energy and enthalpy of formation from the elements and entropy at 25 \eqn{^{\circ}}{°}C and 1 bar}
  \item{eos}{dataframe, values of the equations-of-state parameters}
  \item{contrib}{character, which contributions to consider in the revised HKF equations equations of state: (\code{n})onsolvation, (\code{s})olvation (the \eqn{\omega}{omega} terms), or (o)rigination contributions (i.e., the property itself at 25 \eqn{^{\circ}}{°}C and 1 bar). Default is \code{c("n","s","o")}, for all contributions}
  \item{H2O.PT}{dataframe, values of the electrostatic properties of water at the temperature(s) and pressure(s) of interest}
  \item{H2O.PrTr}{dataframe, values of the electrostatic properties of water at the reference temperature and pressure}
  \item{domega}{logical, calculate the T and P derivatives of omega?}
  \item{rhohat}{numeric, density of water (g/cm3)}
  \item{Tc}{numeric, temperature (\eqn{^{\circ}}{°}C)}
  \item{alpha}{numeric, coefficient of isobaric expansivity (K^-1)}
  \item{daldT}{numeric, temperature derivative of coefficient of isobaric expansivity (K^-2)}
  \item{beta}{numeric, coefficient of isothermal compressibility (bar^-1)}
}

\details{

  The equations of state permit the calculation of the standard molal properties of species as a function of temperature and pressure. For interactive use, \code{\link{subcrt}} is usually more convenient than calling these functions directly. 

  The \code{property} argument is required and refers to one or more of \samp{G}, \samp{H}, \samp{S}, \samp{Cp} and \samp{V}, and for aqueous species only, \samp{kT} and \samp{E}. The units of these properties are the first ones shown in the description for \code{\link{subcrt}}. The names of the properties are matched without regard to case. 

  The revised HKF equations of state (Helgeson et al., 1981; Tanger and Helgeson, 1988; Shock and Helgeson, 1988) are incorporated in \code{hkf}. The equations-of-state parameters are \code{a1}, \code{a2}, \code{a3}, \code{a4}, \code{c1}, \code{c2}, \code{omega} and \code{Z}; the units of these parameters are as indicated for \code{\link{thermo}$obigt}, \emph{sans} the order of magnitude multipliers. Note that the equation-of-state parameter \code{Z} (appearing in the \eqn{g}{g}-function for the temperature derivatives of the omega parameter; Shock et al., 1992) is taken from \code{thermo$obigt} and \emph{not} from the \code{\link{makeup}} of the species, although in most cases the two values are coincident. \code{H2O.PT} and \code{H2O.PrTr} are optional arguments that contain the electrostatic properties of \eqn{\mathrm{H_2O}}{H2O} required for the calculations. If either of these is \code{NULL} (the default), the required values are retrieved using \code{\link{water}}. 

  Unless \code{domega}, the value of which is recycled to the number of species (rows of ghs and eos), is FALSE for any of the species, the temperature and pressure derivatives of the \code{omega} parameter for charged species (where \code{Z != 0}) are calculated using the \eqn{g}{g}- and \eqn{f}{f}-functions described by Shock et al., 1992 and Johnson et al., 1992, and coded here in \code{gfun}. This option is currently blocked when the IAPWS-95 equations are activated (see \code{\link{water}}).

  The parameters in the \code{cgl} equations of state for crystalline, gas and liquid species (except liquid water) include \code{V}, \code{a}, \code{b}, \code{c}, \code{d}, \code{e}, \code{f} and \code{lambda}. The terms denoted by \code{a}, \code{b} and \code{c} correspond to the Maier-Kelley equation for heat capacity (Maier and Kelley, 1932); the additional terms are useful for representing heat capacities of minerals (Robie and Hemingway, 1995) and gaseous or liquid organic species (Helgeson et al., 1998). The standard molal volumes (\samp{V}) of species in these calculations are taken to be independent of temperature and pressure.

For both \code{hkf} and \code{cgl}, if at least one equations-of-state parameter for a species is provided, any NA values of the other parameters are reset to zero.
If all equations-of-state parameters are NA, but values of \samp{Cp} and/or \samp{V} are available, those values are used in the integration of \samp{G}, \samp{H} and \samp{S} as a function of temperature. 

The \code{T} and \code{P} arguments (and \code{rhohat}, \code{Tc}, \code{alpha}, \code{daldT}, \code{beta} for \code{gfun}) should all be the same length; the functions perform no argument validating.
}

\section{Warning}{
The temperature and pressure range of validity of the revised HKF equations of state for aqueous species corresponds to the stability region of liquid water or the supercritical fluid at conditions between 0 to 1000 \eqn{^{\circ}}{°}C and 1 to 5000 bar (Tanger and Helgeson, 1988; Shock and Helgeson, 1988).
The \code{hkf} function does not check these limits and will compute properties as long as the requisite electrostatic properties of water are available. There are conceptually no temperature limits (other than 0 Kelvin) for the validity of the \code{cgl} equations of state.
However, the actual working upper temperature limits correspond to the temperatures of phase transitions of minerals or to those temperatures beyond which extrapolations from experimental data become highly uncertain.
These temperature limits are stored in the thermodynamic database for some minerals, but \code{cgl} ignores them; however, \code{\link{subcrt}} warns if they are exceeded.
}

\value{
  A list of length equal to the number of species (i.e., number rows of supplied \code{ghs} and \code{eos} values). Each element of the list contains a dataframe, each column of which corresponds to one of the specified properties; the number of rows is equal to the number of pressure-temperature points.
}

\seealso{
  \code{\link{info}} for retrieving equations of state parameters from the thermodynamic database, \code{\link{water}} for equations of state of water, \code{\link{subcrt}} for calculations that use these equations.
}


\examples{
\dontshow{data(thermo)}
## aqueous species
a <- info(info("methane","aq"))	
hkf(property="Cp",ghs=a,eos=a)
# the non-solvation heat capacity
hkf(property="Cp",ghs=a,eos=a,contrib="n")	
# at different temperature and pressure
hkf(property="Cp",ghs=a,eos=a,T=c(373.15,473.15),P=1000)

## crystalline, gas, liquid species
a <- info(info("methane","gas"))	
cgl(property="Cp",ghs=a,eos=a)
# melting and vaporization of n-octane
a <- info(info(rep("n-octane",3),c("cr","liq","gas")))
b <- cgl(property="G",ghs=a,eos=a,T=seq(200,420,10),P=1)
which.pmax(b,pmin=TRUE)  # 1 = cr, 2 = liq, 3 = gas
# compare that result with the tabulated transition temperatures
print(a)
}

\references{

  Helgeson, H. C., Kirkham, D. H. and Flowers, G. C. (1981) Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures. IV. Calculation of activity coefficients, osmotic coefficients, and apparent molal and standard and relative partial molal properties to 600\eqn{^{\circ}}{°}C and 5 Kb. \emph{Am. J. Sci.} \bold{281}, 1249--1516. \url{http://www.ajsonline.org/cgi/content/abstract/281/10/1249}

  Helgeson, H. C., Owens, C. E., Knox, A. M. and Richard, L. (1998) Calculation of the standard molal thermodynamic properties of crystalline, liquid, and gas organic molecules at high temperatures and pressures. \emph{Geochim. Cosmochim. Acta} \bold{62}, 985--1081. \url{http://dx.doi.org/10.1016/S0016-7037(97)00219-6}

  Maier, C. G. and Kelley, K. K. (1932) An equation for the representation of high-temperature heat content data. \emph{J. Am. Chem. Soc.} \bold{54}, 3243--3246. \url{http://dx.doi.org/10.1021/ja01347a029}

  Robie, R. A. and Hemingway, B. S. (1995) \emph{Thermodynamic Properties of Minerals and Related Substances at 298.15 K and 1 Bar (\eqn{10^5} Pascals) Pressure and at Higher Temperatures}. U. S. Geol. Surv., Bull. 2131, 461 p. \url{http://www.worldcat.org/oclc/32590140}

  Shock, E. L. and Helgeson, H. C. (1988) Calculation of the thermodynamic and transport properties of aqueous species at high pressures and temperatures: Correlation algorithms for ionic species and equation of state predictions to 5 kb and 1000\eqn{^{\circ}}{°}C. \emph{Geochim. Cosmochim. Acta} \bold{52}, 2009--2036. \url{http://dx.doi.org/10.1016/0016-7037(88)90181-0}
  
  Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \eqn{^{\circ}}{°}C and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \url{http://dx.doi.org/10.1039/FT9928800803}

  Tanger, J. C. IV and Helgeson, H. C. (1988) Calculation of the thermodynamic and transport properties of aqueous species at high pressures and temperatures: Revised equations of state for the standard partial molal properties of ions and electrolytes. \emph{Am. J. Sci.} \bold{288}, 19--98. \url{http://www.ajsonline.org/cgi/content/abstract/288/1/19}

}

\concept{Secondary thermodynamic modeling}
