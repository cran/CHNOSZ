\encoding{UTF-8}
\name{util.units}
\alias{util.units}
\alias{P.units}
\alias{T.units}
\alias{E.units}
\alias{convert}
\alias{envert}
\alias{outvert}
\title{Functions to Convert Units}
\description{
  These functions to convert values between units and set the user's preferred units.
}

\usage{
  P.units(units = NULL)
  T.units(units = NULL)
  E.units(units = NULL)
  convert(value, units, T = get("thermo")$opt$Tr,
    P = get("thermo")$opt$Pr, pH = 7, logaH2O = 0)
  envert(value, units)
  outvert(value, units)
}

\arguments{
  \item{units}{character, name of units to set or convert to/from}
  \item{value}{numeric, value(s) to be converted}
  \item{T}{numeric, temperature (Kelvin), used in \samp{G}-\samp{logK}, \samp{pe}-\samp{Eh} and \samp{logfO2}-\samp{E0} conversions}
  \item{P}{numeric, pressure (bar), used in \samp{logfO2}-\samp{E0} conversions}
  \item{pH}{numeric, pH, used in \samp{logfO2}-\samp{E0} conversions}
  \item{logaH2O}{numeric, logarithm of activity of water, used in \samp{logfO2}-\samp{E0} conversions}
}

\details{

   The units settings are used by \code{\link{subcrt}}, \code{\link{affinity}}, and \code{\link{diagram}} to accept input in or convert output to the units desired by the user. The settings, which can be queried or changed with \code{T.units}, \code{E.units} and \code{P.units}, refer to the units of temperature (\code{K} or \code{C}), energy (\code{cal} or \code{J}), and pressure (\code{bar}, \code{MPa}). (The first value in each of those pairs refers to the default units).

The actual units conversions are handled by \code{convert}, through which \code{values} are transformed into destination \code{units} (names not case sensitive).
The possible conversions and settings for the \code{units} argument are shown in the following table.
Note that \samp{Eh} and \samp{E0} both stand for the value of Eh (oxidation-reduction potential in volts); they have different names so that one can choose to convert between Eh and either \samp{pe} or \samp{logfO2}. 

  \tabular{lll}{
     property \tab units \tab setting of \code{units} argument \cr
     temperature \tab \eqn{^{\circ}}{°}C, K \tab \code{C}, \code{K} \cr
     pressure \tab bar, MPa \tab \code{bar}, \code{MPa} \cr
     energy \tab cal, J \tab \code{cal}, \code{J} \cr
     energy \tab cal, cm\eqn{^3}{^3} bar \tab \code{calories}, \code{cm3bar} \cr
     energy \tab cal, [none] \tab \code{G}, \code{logK} \cr
     oxidation potential \tab volt, [none] \tab \code{Eh}, \code{pe} \cr
     oxidation potential \tab volt, [none] \tab \code{E0}, \code{logfO2} \cr
  }

  \code{envert} and \code{outvert} are wrappers for \code{convert} that handle the conditional conversion of values from or to the current units settings. \code{envert} converts the value \emph{to} the units in the argument, and \code{outvert} converts the value \emph{from} the units in the argument, only if they are different than the current setting; otherwise the value is returned unchanged.

}

\value{
  Character return for \code{T.units}, \code{P.units} and \code{E.units}; numeric returns by the other functions.
}


\examples{\dontshow{data(thermo)}
## examples using convert
# temperature (Kelvin) to degrees C
convert(273.15, "C")
# temperature (degrees C) to Kelvin
convert(100, "K")
# Gibbs energy (cal mol-1) to/from logK
convert(1000, "logK")		
convert(1000, "logK", T=373.15)
convert(1, "G") 
# Eh (volt) to pe 
convert(-1, "pe")		
convert(-1, "pe", T=373.15)
# logfO2 to E0 (volt)
convert(-80, "E0")             
convert(-80, "E0", pH=5)
convert(-80, "E0", pH=5, logaH2O=-5)
# calorie to/from joule
convert(10, "J")		
convert(10, "cal")
# cm3bar to calories		
convert(10, "calories")

## examples showing unit settings
T.units("C")  
T1in <- envert(25, "C")         # no convertion
T1out <- outvert(313.15, "K")   # K to C
T.units("K") 
T2in <- envert(298.15, "C")     # K to C
T2out <- outvert(-233.15, "C")  # C to K
# these are the same temperature (25 deg C)
stopifnot(all.equal(T1in, T2in))
# these are numerically equivalent (40 deg C / 40 K)
stopifnot(all.equal(T1out, T2out))
T.units("C")  # return to default
}

\keyword{util}
