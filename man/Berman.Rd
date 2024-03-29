\encoding{UTF-8}
\name{Berman}
\alias{Berman}
\title{Thermodynamic Properties of Minerals}
\description{
Calculate thermodynamic properties of minerals using the equations of Berman (1988).
}

\usage{
  Berman(name, T = 298.15, P = 1, check.G = FALSE,
         calc.transition = TRUE, calc.disorder = TRUE)
}

\arguments{
  \item{name}{character, name of mineral}
  \item{T}{numeric, temperature(s) at which to calculate properties (K)}
  \item{P}{numeric, pressure(s) at which to calculate properties (bar)}
  \item{check.G}{logical, check consistency of G, H, and S?}
  \item{calc.transition}{logical, include calculation of polymorphic transition properties?}
  \item{calc.disorder}{logical, include calculation of disordering properties?}
}

\details{
This function calculates the thermodynamic properties of minerals at high \P and \T using equations given by Berman (1988).
These minerals should be listed in \code{thermo()$OBIGT} with the state \samp{cr} and chemical formula, and optionally an abbreviation and references, but all other properties set to NA.

The standard state thermodynamic properties and parameters for the calculations are stored in data files under \code{extdata/Berman}, or can be read from a user-created file specified by \code{\link{thermo}()$opt$Berman}.

The equation used for heat capacity is \Cp = k0 + k1*\T\S{-0.5} + k2*\T\S{-2} + k3*\T\S{-3} + k4*\T\S{-1} + k5*\T + k6*\T\S{2}.
This is an extended form Eq. 4 of Berman (1988) as used in the winTWQ program (Berman, 2007).
The equation used for volume is \I{V}(\P, \T) / \I{V}(1 bar, 298.15 K) =  1 + v1 * (\T - 298.15) + v2 * (\T - 298.15)\S{2} + v3 * (\P - 1) + v4 * (\P - 1)\S{2}  (Berman, 1988, Eq. 5, with terms reordered to follow winTWQ format).
The equations used for lambda transitions follow Eqs. 8-14 of Berman (1988).
The equation used for the disorder contribution between Tmin and Tmax is \Cp[dis] = d0 + d1*\T\S{-0.5} + d2*\T\S{-2} + d3*\T + d4*\T\S{2} (Berman, 1988, Eq. 15).
The parameters correspond to Tables 2 (\code{GfPrTr}, \code{HfPrTr}, \code{SPrTr}, \code{VPrTr}), 3a (\code{k0} to \code{k3}), 4 (\code{v1} to \code{v4}), 3b (transition parameters: \code{Tlambda} to \code{dTH}), and 5 (disorder parameters: \code{Tmax}, \code{Tmin}, \code{d1} to \code{d4} and \code{Vad}) of Berman (1988).
Following the winTWQ data format, multipliers are applied to the volume parameters only (see below).
Note that \code{VPrTr} is tabulated in J bar\S{-1} mol\S{-1}, which is equal to 10 cm\S{3} mol\S{-1}.

A value for \code{GfPrTr} is not required and is only used for optional checks (see below).
Numeric values (possibly 0) should be assigned for all of \code{HfPrTr}, \code{SPrTr}, \code{VPrTr}, \code{k0} to \code{k6} and \code{v1} to \code{v4}.
Missing (or NA) values are permitted for the transition and disorder parameters, for minerals where they are not used.
The data files have the following 30 columns:

    \tabular{lll}{
      \code{name} \tab mineral name (must match an entry with a formula but NA properties in \code{thermo()$OBIGT})\cr
      \code{GfPrTr} \tab standard Gibbs energy at 298.15 K and 1 bar (J mol\S{-1}) (Benson-Helgeson convention)\cr
      \code{HfPrTr} \tab standard enthalpy at 298.15 K and 1 bar (J mol\S{-1})\cr
      \code{SPrTr} \tab standard entropy at 298.15 K and 1 bar (J mol\S{-1} K\S{-1})\cr
      \code{VPrTr} \tab standard volume at 298.15 K and 1 bar (J bar\S{-1}) [1 J bar\S{-1} = 10 cm\S{3}]\cr
      \code{k0} ... \code{k6} \tab k0 (J mol\S{-1} K\S{-1}) to k6 \cr
      \code{v1} \tab v1 (K\S{-1}) * 10\S{5} \cr
      \code{v2} \tab v2 (K\S{-2}) * 10\S{5} \cr
      \code{v3} \tab v3 (bar\S{-1}) * 10\S{5} \cr
      \code{v4} \tab v4 (bar\S{-2}) * 10\S{8} \cr
      \code{Tlambda} \tab \ifelse{latex}{\eqn{T_\lambda}}{\ifelse{html}{\out{<I>T</I><sub>&lambda;</sub>}}{T[lambda]}} (K) \cr
      \code{Tref} \tab \ifelse{latex}{\eqn{T_{ref}}}{\ifelse{html}{\out{<I>T<sub>ref</sub></I>}}{T[ref]}} (K) \cr
      \code{dTdP} \tab d\T / d\P (K bar\S{-1}) \cr
      \code{l1} \tab l1 ((J/mol)\S{0.5} K\S{-1}) \cr
      \code{l2} \tab l2 ((J/mol)\S{0.5} K\S{-2}) \cr
      \code{DtH} \tab \Delta\ifelse{latex}{\eqn{_t}}{\ifelse{html}{\out{<I>T<sub>t</sub></I>}}{[t]}}\I{H} (J mol\S{-1}) \cr
      \code{Tmax} \tab temperature at which phase is fully disordered (\ifelse{latex}{\eqn{T_D}}{\ifelse{html}{\out{<I>T<sub>D</sub></I>}}{T[D]}} in Berman, 1988) (K) \cr
      \code{Tmin} \tab reference temperature for onset of disordering (\ifelse{latex}{\eqn{t}}{\ifelse{html}{\out{<I>t</I>}}{t}} in Berman, 1988) (K) \cr
      \code{d0} ... \code{d4} \tab d0 (J mol\S{-1} K\S{-1}) to d4 \cr
      \code{Vad} \tab constant that scales the disordering enthalpy to volume of disorder (\ifelse{latex}{\eqn{d_5}}{\ifelse{html}{\out{<I>d</I><sub>5</sub>}}{d5}} in Berman, 1988) \cr
    }

The function outputs apparent Gibbs energies according to the Benson-Helgeson convention (\Delta\I{G} = \Delta\I{H} - \T\Delta\I{S}) using the entropies of the elements in the chemical formula of the mineral to calculate \Delta\I{S} (cf. Anderson, 2005).
If \code{check.G} is TRUE, the tabulated value of \code{GfTrPr} (Benson-Helgeson) is compared with that calculated from \code{HfPrTr} - 298.15*\code{DSPrTr} (\code{DSPrTr} is the difference between the entropies of the elements in the formula and \code{SPrTr} in the table).
A warning is produced if the absolute value of the difference between tabulated and calculated \code{GfTrPr} is greater than 1000 J/mol.

If the function is called with missing \code{name}, the parameters for all available minerals are returned.
}

\value{
A data frame with \code{T} (K), \code{P} (bar), \code{G}, \code{H}, \code{S}, and \code{Cp} (energetic units in Joules), and \code{V} (cm\S{3} mol\S{-1}).
}

\examples{
\dontshow{reset()}
# Other than the formula, the parameters aren't stored in
# thermo()$OBIGT, so this shows NAs
info(info("quartz", "cr"))
# Properties of alpha-quartz (aQz) at 298.15 K and 1 bar
Berman("quartz")
# Gibbs energies of aQz and coesite at higher T and P
T <- seq(200, 1300, 100)
P <- seq(22870, 31900, length.out = length(T))
G_aQz <- Berman("quartz", T = T, P = P)$G
G_Cs <- Berman("coesite", T = T, P = P)$G
# That is close to the univariant curve (Ber88 Fig. 4),
# so the difference in G is close to 0
DGrxn <- G_Cs - G_aQz
all(abs(DGrxn) < 100)  # TRUE

# Make a P-T diagram for SiO2 minerals (Ber88 Fig. 4)
basis(c("SiO2", "O2"), c("cr", "gas"))
species(c("quartz", "quartz,beta", "coesite"), "cr")
a <- affinity(T = c(200, 1700, 200), P = c(0, 50000, 200))
diagram(a)

## Getting data from a user-supplied file
## Ol-Opx exchange equilibrium, after Berman and Aranovich, 1996
species <- c("fayalite", "enstatite", "ferrosilite", "forsterite")
coeffs <- c(-1, -2, 2, 1)
T <- seq(600, 1500, 50)
Gex_Ber88 <- subcrt(species, coeffs, T = T, P = 1)$out$G
# Add data from BA96
datadir <- system.file("extdata/Berman/testing", package = "CHNOSZ")
add.OBIGT(file.path(datadir, "BA96_OBIGT.csv"))
thermo("opt$Berman" = file.path(datadir, "BA96_Berman.csv"))
Gex_BA96 <- subcrt(species, coeffs, T = seq(600, 1500, 50), P = 1)$out$G
# Ber88 is lower than BA96 at low T
(Gex_BA96 - Gex_Ber88)[1] > 0  # TRUE
# The curves cross at about 725 deg C (BA96 Fig. 8)
# (actually, in our calculation they cross closer to 800 deg C)
T[which.min(abs(Gex_BA96 - Gex_Ber88))]  # 800
# Reset the database (thermo()$OBIGT, and thermo()$opt$Berman)
reset()
}

\references{
Anderson, G. M. (2005) \emph{Thermodynamics of Natural Systems}, 2nd ed., Cambridge University Press, 648 p. \url{https://www.worldcat.org/oclc/474880901}

Berman, R. G. (1988) Internally-consistent thermodynamic data for minerals in the system Na{\s2}O-K{\s2}O-CaO-MgO-FeO-Fe{\s2}O{\s3}-Al{\s2}O{\s3}-SiO{\s2}-TiO{\s2}-H{\s2}O-CO{\s2}. \emph{J. Petrol.} \bold{29}, 445-522. \doi{10.1093/petrology/29.2.445}

Berman, R. G. and Aranovich, L. Ya. (1996) Optimized standard state and solution properties of minerals. I. Model calibration for olivine, orthopyroxene, cordierite, garnet, and ilmenite in the system FeO-MgO-CaO-Al{\s2}O{\s3}-TiO{\s2}-SiO{\s2}. \emph{Contrib. Mineral. Petrol.} \bold{126}, 1-24. \doi{10.1007/s004100050233}

Berman, R. G. (2007) winTWQ (version 2.3): A software package for performing internally-consistent thermobarometric calculations. \emph{Open File} \bold{5462}, Geological Survey of Canada, 41 p. \doi{10.4095/223425}

Helgeson, H. C., Delany, J. M., Nesbitt, H. W. and Bird, D. K. (1978) Summary and critique of the thermodynamic properties of rock-forming minerals. \emph{Am. J. Sci.} \bold{278-A}, 1--229. \url{https://www.worldcat.org/oclc/13594862}
}

\concept{Thermodynamic calculations}
