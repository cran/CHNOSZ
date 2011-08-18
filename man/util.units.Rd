\name{util.units}
\alias{util.units}
\alias{nuts}
\alias{convert}
\alias{envert}
\alias{outvert}
\title{Functions to Convert Units}
\description{
  These functions to convert values between units and set the user's preferred units.
}

\usage{
  nuts(units = NULL)
  convert(value, units, T = thermo$opt$Tr,
    P = thermo$opt$Pr, pH = 7, logaH2O = 0)
  envert(value,units)
  outvert(value,units)
}

\arguments{
  \item{units}{character, name of units to set in global preferences (\code{nuts}), or name of units (see table) to convert to or from.}
  \item{value}{numeric, value(s) to be converted.}
  \item{T}{numeric, temperature (Kelvin), used in \samp{G}-\samp{logK}, \samp{pe}-\samp{Eh} and \samp{logfO2}-\samp{E0} conversions.}
  \item{P}{numeric, pressure (bar), used in \samp{logfO2}-\samp{E0} conversions.}
  \item{pH}{numeric, pH, used in \samp{logfO2}-\samp{E0} conversions.}
  \item{logaH2O}{numeric, logarithm of activity of water, used in \samp{logfO2}-\samp{E0} conversions.}
}

\details{

   The global units settings are used by \code{\link{subcrt}}, \code{\link{affinity}}, and \code{\link{diagram}} to accept input in or convert output to the units desired by the user. The settings, which can be queried or changed with \code{nuts}, refer to the units of temperature (\code{K} or \code{C}), energy (\code{cal} or \code{J}), and pressure (\code{bar}, \code{MPa}). (The first value in each of those pairs refers to the default units). If the \code{units} argument to \code{nuts} is NULL (the default) the current units settings are printed (with a NULL return); if it is one of \samp{T}, \samp{P} or \samp{E}, the units for that property are printed and also returned.

  The actual units conversions are handled by \code{convert}, through which \code{values} are transformed into destination \code{units} (names not case sensitive). The possible unit conversions are shown in the following table. Note that \samp{Eh} and \samp{E0} both stand for the value of Eh (do not confuse \samp{E0} with the standard reduction potential, which is written similarly); they have different names so that one can choose to convert to \samp{pe} or \samp{logfO2}. 

  \tabular{lll}{
     \code{C} \tab \code{K} \tab temperature \cr
     \code{cal} \tab \code{J}  \tab energy \cr
     \code{bar} \tab \code{MPa} \tab pressure \cr
     \code{E0} \tab \code{logfO2} \tab oxidation state \cr
     \code{G} \tab \code{logK} \tab energy \cr
     \code{cm3bar} \tab \code{calories} \tab energy \cr
     \code{pe} \tab \code{Eh}  \tab oxidation state \cr
  }

  \code{envert} and \code{outvert} are wrappers for \code{convert} that handle the conditional conversion of values from or to the user's units. The name of the function \code{\link{nuts}} is such to avoid clashing with the name of \code{\link{units}}, which is present in the \R base package.

}

\value{
  Numeric, except for \code{\link{nuts}}, which returns character or \code{NULL}.
}


\examples{
  \dontshow{data(thermo)}
  ## examples using convert
  # temperature (Kelvin) to degrees C
  convert(273.15,"C")
  # temperature (degrees C) to Kelvin
  convert(100,"K")
  # Gibbs energy (cal mol-1) to/from logK
  convert(1000,"logK")		
  convert(1000,"logK",T=373.15)
  convert(1,"G") 
  # Eh (volt) to pe 
  convert(-1,"pe")		
  convert(-1,"pe",T=373.15)
  # logfO2 to E0 (volt)
  convert(-80,"E0")             
  convert(-80,"E0",pH=5)
  convert(-80,"E0",pH=5,logaH2O=-5)
  # calorie to/from joule
  convert(10,"j")		
  convert(10,"cal")
  # cm3bar to calories		
  convert(10,"calories")

  ## examples showing user unit preference	
  nuts("K")  # set temperature units
  nuts("J")  # set energy units
  # return the name of the (E)nergy units
  nuts("E")	
  # print names of all units (NULL return)		
  nuts()	
  # defaults		
  nuts("K"); nuts("bar"); nuts("cal")	
}

\keyword{util}
