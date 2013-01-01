\encoding{UTF-8}
\name{water}
\alias{water}
\alias{water.AW90}
\alias{water.IAPWS95}
\alias{idealgas.IAPWS95}
\alias{residual.IAPWS95}
\alias{water.WP02}
\alias{water.SUPCRT92}
\title{Properties of Water}
\description{
  Calculate thermodynamic and electrostatic properties of water.
}

\usage{
  water(property = NULL, T = thermo$opt$Tr, P = "Psat")
  water.SUPCRT92(property, T = 298.15, P = 1, isat = 0)
  water.IAPWS95(property, T = 298.15, rho = 1000)
  idealgas.IAPWS95(p, delta, tau)
  residual.IAPWS95(p, delta, tau)
  water.WP02(property, T = 298.15)
  water.AW90(T = 298.15, rho = 1000, P = 0.1)
}

\arguments{
  \item{property}{character, name(s) of property(s) to calculate}
  \item{T}{numeric, temperature (K)}
  \item{P}{character or numeric, \samp{Psat}; or pressure (bar for \code{water}, \code{water.SUPCRT92}; MPa for \code{water.AW90}, \code{water.IAPWS95}, \code{water.WP02})}
  \item{isat}{numeric, if \eqn{1}, calculate values of \samp{Psat}}
  \item{rho}{numeric, density (kg m\eqn{^{-3}}{^-3})}
  \item{p}{character, name of property (Helmholtz free energy or its derivatives)}
  \item{delta}{numeric, density divided by critical density}
  \item{tau}{numeric, critical temperature divided by temperature}
}

\details{

  These functions compute the thermodynamic (Gibbs energy and it derivatives) and electrostatic (dielectric constant and its derivatives) properties of liquid or supercritical \eqn{\mathrm{H_2O}}{H2O} using equations of state taken from the literature. The wrapper function \code{water} accepts two major computational alternatives. The default option (for \code{\link{thermo}$opt$water} equal to \samp{SUPCRT}) indicates to retrieve thermodynamic and electrostatic properties as a function of temperature and pressure using a FORTRAN subroutine taken from the \acronym{SUPCRT92} software package (Johnson et al., 1992). If \code{thermo$opt$water} is set to \samp{IAPWS}, the thermodynamic properties are calculated using a (slower) implementation in \R code of the IAPWS-95 formulation (Wagner and Pruss, 2002) and electrostatic properties are calculated using the equations of Archer and Wang, 1990.

  The allowed \code{property}s for \code{water} are one or more of those given below, depending on the computational option (availability is shown by an asterisk). The names of properties in the arguments are not case sensitive, and some of the properties that can actually be calculated using the equations of state are not implemented here.

  \tabular{lllll}{
     Property \tab Description \tab Units \tab IAPWS \tab SUPCRT \cr
     \code{A} \tab Helmholtz energy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{G} \tab Gibbs energy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{S} \tab Entropy \tab cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{U} \tab Internal energy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{H} \tab Enthalpy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{Cv} \tab Isochoric heat capacity \tab cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{Cp} \tab Isobaric heat capacity \tab cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{w} (\code{Speed}) \tab Speed of sound \tab cm s\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{E} \tab Isobaric expansivity \tab cm\eqn{^3} K\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{kT} \tab Isothermal compressibility \tab cm\eqn{^3} bar\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{alpha} \tab Coefficient of isobaric expansivity \tab K\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{beta} \tab Coefficient of isothermal compressibility \tab bar\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{epsilon} (\code{diel}) \tab Dielectric constant \tab dimensionless \tab NA \tab * \cr
     \code{visc} \tab Dynamic viscosity \tab g cm\eqn{^{-1}}{^-1} s\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{tcond} \tab Thermal conductivity \tab cal cm\eqn{^{-1}}{^-1} s\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{tdiff} \tab Thermal diffusivity \tab cm\eqn{^2} s\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{Prndtl} \tab Prandtl number \tab dimensionless \tab NA \tab * \cr
     \code{visck} \tab Kinematic viscosity \tab cm\eqn{^2} s\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{albe} \tab Isochoric expansivity \tab bar K\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \tab -compressibility \tab \tab \tab \cr
     \code{Z} (\code{ZBorn}) \tab Z Born function \tab dimensionless \tab NA \tab * \cr
     \code{Y} (\code{YBorn}) \tab Y Born function \tab K\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{Q} (\code{QBorn}) \tab Q Born function \tab bar\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{daldT} \tab Isobaric temperature derivative \tab K\eqn{^{-2}}{^-2} \tab NA \tab * \cr
     \tab of expansibility \tab \tab \tab \cr
     \code{X} (\code{XBorn}) \tab X Born function \tab K\eqn{^{-2}}{^-2} \tab * \tab * \cr
     \code{N} \tab N Born function \tab bar\eqn{^{-2}}{^-2} \tab * \tab NA \cr
     \code{UBorn} \tab U Born function \tab bar\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1} \tab * \tab NA \cr
     \code{V} \tab Volume \tab cm\eqn{^3} mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{rho} \tab Density \tab kg cm\eqn{^3} \tab * \tab * \cr     
     \code{Psat} \tab Saturation vapor pressure \tab bar \tab * \tab * \cr
     \code{P} \tab Pressure \tab bar \tab * \tab NA \cr
     \code{de.dT} \tab Temperature derivative \tab K\eqn{^{-1}}{^-1} \tab * \tab NA \cr
     \tab of dielectric constant \tab \tab \tab \cr
     \code{de.dP} \tab Pressure derivative \tab bar\eqn{^{-1}}{^-1} \tab * \tab NA \cr
     \tab of dielectric constant \tab \tab \tab \cr
  }

  \code{UBorn} refers to the \eqn{U}{U} Born function, \samp{U} to internal energy. The coefficients of isobaric expansivity and of isothermal compressibility are defined as \eqn{(1/V)(dV/dT)_P}{(1/V)(dV/dT)P} and \eqn{-(1/V)(dV/dP)_T}{-(1/V)(dV/dP)T} respectively. All of the properties are calculated as a function of temperature and pressure except \samp{Psat} f(T) (values supplied in the argument \code{P} are ignored) and \samp{P} f(T,rho). Except for those of \code{rho}, the units used are as in Johnson and Norton, 1991. Names of properties that are used in \code{water.SUPCRT92} (but not in \code{water}) are shown in parentheses. 

  \code{water.SUPCRT92} interfaces to the FORTRAN subroutine taken from the \acronym{SUPCRT92} package (H2O92D.F) for calculating properties of water. These calculations are based on data and equations of Levelt-Sengers et al., 1983, Haar et al., 1984, and Johnson and Norton, 1991, among others (see Johnson et al., 1992). If \code{isat} is equal to \eqn{1} (or \code{TRUE}), the values of \code{P} are ignored and values of \samp{Psat} are returned. \samp{Psat} refers to one bar below 100 \eqn{^{\circ}}{°}C, otherwise to the vapor-liquid saturation pressure at temperatures below the critical point (\samp{Psat} is not available at temperatures above the critical point). \code{water.SUPCRT92} function provides a limited interface to the FORTRAN subroutine; some functions provided there are not made available here (e.g., using variable density instead of pressure, or calculating the properties of steam). The properties of steam in CHNOSZ, as in \acronym{SUPCRT92}, are calculated using general equations for crystalline, gaseous and liquid species (\code{\link{cgl}}). The IAPWS-95 formulation also has provisions for computing the properties of steam, but these are currently not used by CHNOSZ.

  \code{water.IAPWS95} provides an implementation of the IAPWS-95 formulation for properties (including pressure) calculated as a function of temperature and density. To compute the thermodynamic and electrostatic properties of water as a function of temperature and pressure using \code{water.IAPWS95}, \code{water} applies a root-finding function (\code{\link{uniroot}}) to determine the corresponding values of density. Electrostatic properties in this case are derived from values of the static dielectric constant (\code{epsilon}) calculated using equations given by Archer and Wang, 1990 and coded in \code{water.AW90}. Note that the \code{water.AW90} computes the static dielectric constant at given temperatures and pressures, so \code{water} contains routines for calculating its derivatives with respect to temperature and pressure.

  The \code{water.IAPWS95} function returns values of thermodynamic properties in specific units (per gram) which are converted to molar properties by \code{water}. The IAPWS-95 formulation follows the triple point convention used in engineering (values of internal energy and entropy are taken to be zero at the triple point). For compatibility with geochemical modeling conventions, the values of Gibbs energy, enthalpy and entropy output by \code{water.IAPWS95} are converted by \code{water} to the triple point reference state adopted in \code{SUPCRT92} (Johnson and Norton, 1991; Helgeson and Kirkham, 1974). Auxiliary equations to the IAPWS-95 formulation (Wagner and Pruss, 2002) are provided in \code{water.WP02}; the \code{property} for this function can be one of \samp{P.sigma} (saturation vapor pressure in MPa), \samp{dP.sigma.dT} (derivative of saturation vapor pressure with respect to temperature), or \samp{rho.liquid} or \samp{rho.vapor} (density of liquid or vapor in kg m\eqn{^{-3}}{^-3}).

\code{idealgas.IAPWS95} and \code{residual.IAPWS95} are supporting functions to \code{water.IAPWS95} for calculating the ideal-gas and residual parts in the IAPWS-95 formulation.
The value of \code{p} can be one of \samp{phi}, \samp{phi.delta}, \samp{phi.delta.delta}, \samp{phi.tau}, \samp{phi.tau.tau}, or \samp{phi.delta.tau}, to calculate the specific dimensionless Helmholtz free energy (\samp{phi}) or one of its respective derivatives.

The stated temperature limits of validity of calculations in \code{water.SUPCRT92} are from the greater of 0 \eqn{^{\circ}}{°}C or the melting temperature at pressure, to 2250 \eqn{^{\circ}}{°}C (Johnson et al., 1992).
For \code{water.IAPWS} the upper temperature limit of validity is 1000 \eqn{^{\circ}}{°}C, but extrapolation to much higher temperatures is possible (Wagner and Pruss, 2002).
Valid pressures are from the greater of zero bar or the melting pressure at temperature to 30000 bar (\code{water.SUPCRT92}) or 10000 bar (\code{water.IAPWS95}; again, with the provision for extrapolation to more extreme conditions).
The present functions do not check these limits and will attempt calculations for any range of input parameters, but may return \code{NA} for properties that fail to be calculated at given temperatures and pressures and/or produce warnings or even errors when problems are encountered.

}


\value{
  For \code{water}, \code{water.SUPCRT92} and \code{water.IAPWS}, a dataframe the number of rows of which corresponds to the number of input temperature, pressure and/or density values. \code{water.AW90} returns a numeric vector with length corresponding to the number of temperature values.
}

\seealso{
  \code{\link{uniroot}} is the root finder used in \code{water} to back out values of the density (\code{rho}) from those of \code{T} and \code{P} when the \samp{IAPWS} option is set in \code{\link{thermo}$opt$water}. Equations of state for species other than water are coded in \code{\link{hkf}} and \code{\link{cgl}}.
}

\examples{
\dontshow{data(thermo)}
## set temperature, density
T <- 500; rho <- 838.0235
# calculate pressure 				
P <- as.numeric(water.IAPWS95('P',T=T,rho=rho))
# calculate dielectric constant 	
water.AW90(T=T,rho=rho,P=P)	
# find water density for this T, P		
water('rho',T=T,P=convert(P,'bar'))

## density along saturation curve
T <- seq(273.15,623.15,25)
water.WP02(T=T)  # liquid from WP02
water.WP02('rho.vapor',T)  # steam from WP02	
water('rho',T=T,P='Psat')  # liquid from SUPCRT92
# values of the density, Psat, Gibbs energy	
water(c('rho','psat','G'),T=T,P='Psat') 
# derivatives of the dielectric constant (Born functions)
water(c('Q','Y','X','U'),T=T)
# now at constant pressure
water(c('Q','Y','X','U'),T=T,P=2000)

## comparing the computational options
prop <- c("A","G","S","U","H","Cv","Cp","w","epsilon",
  "Y","Q","X","rho","Psat")
thermo$opt$water <<- "SUPCRT"
print(water(prop,T=convert(c(25,100,200,300),"K")))
thermo$opt$water <<- "IAPWS"
print(water(c(prop,"N","UBorn"),T=convert(c(25,100,200,300),"K")))
# fixme: things seem to be working except speed of
# sound in our IAPWS calculations

# calculating Q Born function
# after Table 22 of Johnson and Norton, 1991
thermo$opt$water <<- "SUPCRT"
T <- rep(c(375, 400, 425, 450, 475), each=5)
P <- rep(c(250, 300, 350, 400, 450), 5)
w <- water("Q", T=convert(T, "K"), P=P)
# the rest is to make a readable table
w <- as.data.frame(matrix(w[[1]], nrow=5))
colnames(w) <- T[1:5*5]
rownames(w) <- P[1:5]
print(w)
}

\references{

   Archer, D. G. and Wang, P. M. (1990) The dielectric constant of water and Debye-Huckel limiting law slopes. \emph{J. Phys. Chem. Ref. Data} \bold{19}, 371--411. \url{http://www.nist.gov/srd/PDFfiles/jpcrd383.pdf}

  Haar, L., Gallagher, J. S. and Kell, G. S. (1984) \emph{NBS/NRC Steam Tables}. Hemisphere, Washington, D. C., 320 p. \url{http://www.worldcat.org/oclc/301304139}

 Helgeson, H. C. and Kirkham, D. H. (1974) Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures. I. Summary of the thermodynamic/electrostatic properties of the solvent. \emph{Am. J. Sci.} \bold{274}, 1089--1098. \url{http://www.ajsonline.org/cgi/content/abstract/274/10/1089}

  Johnson, J. W. and Norton, D. (1991) Critical phenomena in hydrothermal systems: state, thermodynamic, electrostatic, and transport properties of H\eqn{_2}{2}O in the critical region. \emph{Am. J. Sci.} \bold{291}, 541--648. \url{http://www.ajsonline.org/cgi/content/abstract/291/6/541}

  Johnson, J. W., Oelkers, E. H. and Helgeson, H. C. (1992) SUPCRT92: A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000\eqn{^{\circ}}{°}C. \emph{Comp. Geosci.} \bold{18}, 899--947. \url{http://dx.doi.org/10.1016/0098-3004(92)90029-Q}

  Levelt-Sengers, J. M. H., Kamgarparsi, B., Balfour, F. W. and Sengers, J. V. (1983) Thermodynamic properties of steam in the critical region. \emph{J. Phys. Chem. Ref. Data} \bold{12}, 1--28. \url{http://www.nist.gov/srd/PDFfiles/jpcrd214.pdf}

  Wagner, W. and Pruss, A. (2002) The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. \emph{J. Phys. Chem. Ref. Data} \bold{31}, 387--535. \url{http://dx.doi.org/10.1063/1.1461829}

}


\keyword{secondary}
