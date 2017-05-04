\encoding{UTF-8}
\name{water}
\alias{water}
\alias{water.props}
\alias{water.SUPCRT92}
\alias{water.IAPWS95}
\title{Properties of Water}
\description{
Calculate thermodynamic and electrostatic properties of water.
}

\usage{
  water(property = NULL, T = get("thermo")$opt$Tr, P = "Psat")
  water.props(formulation = get("thermo")$opt$water)
  water.SUPCRT92(property, T = 298.15, P = 1)
  water.IAPWS95(property, T = 298.15, P = 1)
}

\arguments{
  \item{property}{character, name(s) of property(s) to calculate}
  \item{T}{numeric, temperature (K)}
  \item{P}{numeric, pressure (bar), or \samp{Psat} for vapor-liquid saturation}
  \item{formulation}{character, name of formulation for which to return names of available properties}
}

\details{

These functions compute the thermodynamic (Gibbs energy and it derivatives) and electrostatic (dielectric constant and its derivatives) properties of liquid or supercritical \H2O using equations of state taken from the literature.
The high-level function \code{water} accepts two major computational alternatives.
The default option (i.e., \code{\link{thermo}$opt$water} equal to \samp{SUPCRT92}) is retrieve thermodynamic and electrostatic properties as a function of temperature and pressure using a FORTRAN subroutine taken from the \acronym{SUPCRT92} software package (Johnson et al., 1992).
If \code{thermo$opt$water} is set to \samp{IAPWS95}, the thermodynamic properties are calculated using an implementation in \R code of the IAPWS-95 formulation (Wagner and Pruss, 2002), and electrostatic properties are calculated using the equations of Archer and Wang, 1990.

The allowed \code{property}s for \code{water} are one or more of those given below, depending on the computational option; availability is shown by an asterisk.
The names of properties in the arguments are not case sensitive. Note that some of the properties that can actually be calculated using the different formulations are not implemented here.
Except for \code{rho}, the units are those used by Johnson and Norton, 1991.

  \tabular{lllll}{
     Property \tab Description \tab Units \tab IAPWS95 \tab SUPCRT92 \cr
     \code{A} \tab Helmholtz energy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{G} \tab Gibbs energy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{S} \tab Entropy \tab cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{U} \tab Internal energy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{H} \tab Enthalpy \tab cal mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{Cv} \tab Isochoric heat capacity \tab cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{Cp} \tab Isobaric heat capacity \tab cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{Speed} \tab Speed of sound \tab cm s\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{alpha} \tab Coefficient of isobaric expansivity \tab K\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{beta} \tab Coefficient of isothermal compressibility \tab bar\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{diel} \tab Dielectric constant \tab dimensionless \tab NA \tab * \cr
     \code{visc} \tab Dynamic viscosity \tab g cm\eqn{^{-1}}{^-1} s\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{tcond} \tab Thermal conductivity \tab cal cm\eqn{^{-1}}{^-1} s\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{tdiff} \tab Thermal diffusivity \tab cm\eqn{^2} s\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{Prndtl} \tab Prandtl number \tab dimensionless \tab NA \tab * \cr
     \code{visck} \tab Kinematic viscosity \tab cm\eqn{^2} s\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{albe} \tab Isochoric expansivity \tab bar K\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \tab -compressibility \tab \tab \tab \cr
     \code{ZBorn} \tab Z Born function \tab dimensionless \tab NA \tab * \cr
     \code{YBorn} \tab Y Born function \tab K\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{QBorn} \tab Q Born function \tab bar\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{daldT} \tab Isobaric temperature derivative \tab K\eqn{^{-2}}{^-2} \tab NA \tab * \cr
     \tab of expansibility \tab \tab \tab \cr
     \code{XBorn} \tab X Born function \tab K\eqn{^{-2}}{^-2} \tab * \tab * \cr
     \code{NBorn} \tab N Born function \tab bar\eqn{^{-2}}{^-2} \tab * \tab NA \cr
     \code{UBorn} \tab U Born function \tab bar\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1} \tab * \tab NA \cr
     \code{V} \tab Volume \tab cm\eqn{^3} mol\eqn{^{-1}}{^-1} \tab * \tab * \cr
     \code{rho} \tab Density \tab kg cm\eqn{^3} \tab * \tab * \cr     
     \code{Psat} \tab Saturation vapor pressure \tab bar \tab * \tab * \cr
     \code{E} \tab Isobaric expansivity \tab cm\eqn{^3} K\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{kT} \tab Isothermal compressibility \tab cm\eqn{^3} bar\eqn{^{-1}}{^-1} \tab NA \tab * \cr
     \code{de.dT} \tab Temperature derivative \tab K\eqn{^{-1}}{^-1} \tab * \tab NA \cr
     \tab of dielectric constant \tab \tab \tab \cr
     \code{de.dP} \tab Pressure derivative \tab bar\eqn{^{-1}}{^-1} \tab * \tab NA \cr
     \tab of dielectric constant \tab \tab \tab \cr
     \code{P} \tab Pressure \tab bar \tab * \tab NA \cr
  }

\code{water.props} returns the names of the available properties listed in this table, reflecting the current setting of \code{thermo$opt$water}.

\code{water.SUPCRT92} interfaces to the FORTRAN subroutine taken from the \acronym{SUPCRT92} package (H2O92D.F) for calculating properties of water.
These calculations are based on data and equations of Levelt-Sengers et al., 1983, Haar et al., 1984, and Johnson and Norton, 1991, among others (see Johnson et al., 1992).
A value of \code{P} set to \samp{Psat} refers to one bar below 100 \eqn{^{\circ}}{°}C, otherwise to the vapor-liquid saturation pressure at temperatures below the critical point (\samp{Psat} is not available at temperatures above the critical point).
\code{water.SUPCRT92} provides a limited interface to the FORTRAN subroutine; some functions provided there are not made available here (e.g., using variable density instead of pressure, or calculating the properties of steam).
The properties of steam in CHNOSZ, as in \acronym{SUPCRT92}, are calculated using general equations for crystalline, gaseous and liquid species (\code{\link{cgl}}).

The stated temperature limits of validity of calculations in \code{water.SUPCRT92} are from the greater of 0 \eqn{^{\circ}}{°}C or the melting temperature at pressure, to 2250 \eqn{^{\circ}}{°}C (Johnson et al., 1992).
Valid pressures are from the greater of zero bar or the melting pressure at temperature to 30000 bar.
The present functions do not check these limits and will attempt calculations for any range of input parameters, but may return \code{NA} for properties that fail to be calculated at given temperatures and pressures and/or produce warnings or even errors when problems are encountered.

Starting with version 0.9-9.4, a check for minimum pressure (in valTP function in H2O92D.f) has been bypassed so that properties of H2O can be calculated using \code{water.SUPCRT92} at temperatures below the 0.01 \eqn{^{\circ}}{°}C triple point.
A primary check is still enforced (Tbtm), giving a minimum valid temperature of 253.15 K.

\code{water.IAPWS95} is a wrapper around \code{\link{IAPWS95}}, \code{\link{rho.IAPWS95}} and \code{\link{water.AW90}}.
\code{water.IAPWS95} provides for calculations at specific temperature and pressure; density, needed for \code{IAPWS95}, is inverted from pressure using \code{rho.IAPWS95}.
The function also contains routines for calculating the Born functions as numerical derivatives of the static dielectric constant (from \code{water.AW90}).
For compatibility with geochemical modeling conventions, the values of Gibbs energy, enthalpy and entropy output by \code{IAPWS95} are converted by \code{water.IAPWS95} to the triple point reference state adopted in \code{SUPCRT92} (Johnson and Norton, 1991; Helgeson and Kirkham, 1974).
\code{water.IAPWS95} also accepts setting \code{P} to \samp{Psat}, with the saturation pressure calculated from \code{\link{WP02.auxiliary}}; by default the returned properties are for the liquid, but this can be changed to the vapor in \code{\link{thermo}$opt$IAPWS.sat}.
}


\value{
For \code{water} and \code{water.SUPCRT92} a data frame the number of rows of which corresponds to the number of input temperature, pressure and/or density values.
}

\seealso{
Equations for thermodynamic properties of species other than water are coded in \code{\link{hkf}} and \code{\link{cgl}}.
}

\examples{\dontshow{data(thermo)}
## calculations along saturation curve
T <- seq(273.15, 623.15, 25)
# liquid density, from SUPCRT92
water("rho", T=T, P="Psat")
# values of the saturation pressure, Gibbs energy	
water(c("Psat", "G"), T=T, P="Psat") 
# derivatives of the dielectric constant (Born functions)
water(c("QBorn", "YBorn", "XBorn"), T=T, P="Psat")
# now at constant pressure
water(c("QBorn", "YBorn", "XBorn"), T=T, P=2000)

## comparing the formulations
T <- convert(c(25, 100, 200, 300), "K")
# IAPWS-95
thermo$opt$water <<- "IAPWS95"
water(water.props(), T=T)
# SUPCRT92 (the default)
thermo$opt$water <<- "SUPCRT92"
water(water.props(), T=T)

## calculating Q Born function
# after Table 22 of Johnson and Norton, 1991
thermo$opt$water <<- "SUPCRT92"
T <- rep(c(375, 400, 425, 450, 475), each=5)
P <- rep(c(250, 300, 350, 400, 450), 5)
w <- water("QBorn", T=convert(T, "K"), P=P)
# the rest is to make a neat table
w <- as.data.frame(matrix(w[[1]], nrow=5))
colnames(w) <- T[1:5*5]
rownames(w) <- P[1:5]
print(w)
}

\references{
Archer, D. G. and Wang, P. M. (1990) The dielectric constant of water and Debye-Huckel limiting law slopes. \emph{J. Phys. Chem. Ref. Data} \bold{19}, 371--411. \url{https://doi.org/10.1063/1.555853}

Haar, L., Gallagher, J. S. and Kell, G. S. (1984) \emph{NBS/NRC Steam Tables}. Hemisphere, Washington, D. C., 320 p. \url{http://www.worldcat.org/oclc/301304139}

Helgeson, H. C. and Kirkham, D. H. (1974) Theoretical prediction of the thermodynamic behavior of aqueous electrolytes at high pressures and temperatures. I. Summary of the thermodynamic/electrostatic properties of the solvent. \emph{Am. J. Sci.} \bold{274}, 1089--1098. \url{http://www.ajsonline.org/cgi/content/abstract/274/10/1089}

Johnson, J. W. and Norton, D. (1991) Critical phenomena in hydrothermal systems: state, thermodynamic, electrostatic, and transport properties of H\eqn{_2}{2}O in the critical region. \emph{Am. J. Sci.} \bold{291}, 541--648. \url{http://www.ajsonline.org/cgi/content/abstract/291/6/541}

Johnson, J. W., Oelkers, E. H. and Helgeson, H. C. (1992) SUPCRT92: A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000\eqn{^{\circ}}{°}C. \emph{Comp. Geosci.} \bold{18}, 899--947. \url{https://doi.org/10.1016/0098-3004(92)90029-Q}

Levelt-Sengers, J. M. H., Kamgarparsi, B., Balfour, F. W. and Sengers, J. V. (1983) Thermodynamic properties of steam in the critical region. \emph{J. Phys. Chem. Ref. Data} \bold{12}, 1--28. \url{https://doi.org/10.1063/1.555676}

Wagner, W. and Pruss, A. (2002) The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. \emph{J. Phys. Chem. Ref. Data} \bold{31}, 387--535. \url{https://doi.org/10.1063/1.1461829}
}

\concept{Secondary thermodynamic modeling}
