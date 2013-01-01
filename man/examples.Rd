\encoding{UTF-8}
\name{examples}
\alias{examples}
\alias{demos}
\title{Run Examples from the Documentation}
\description{
  Run the examples contained in each of the documentation topics. }

\arguments{
  \item{do.png}{logical, generate PNG files for the plots?}
  \item{which}{character, which example to run.}
}

\usage{
  examples(do.png = FALSE)
  demos(which = c("sources", "NaCl", "cordierite",
    "phosphate", "nucleobase", "orp", "findit", 
    "CO2Ac", "nonideal"))
}

\details{

\code{examples} runs all the examples in the documentation for the package. 
\code{\link{example}} is called for each topic with \code{ask} set to \code{FALSE} (so all of the figures are shown without prompting the user).
If \code{do.png} is TRUE, the plots in the examples are saved as \code{\link{png}} files having names beginning with the name of each of the help topics.

\code{demos} is a function to run other examples that are provided as \code{\link{demo}s}.
\code{\link{demo}} is called with settings to not echo the source code and to not ask before making each plot.
The demo(s) to run is/are specified by \code{which}; the default is to run them in the order of the list below.
See the comments in the source code for more information about each demo.

  \tabular{ll}{
    \code{sources} \tab cross-check the reference list with the thermodynamic database \cr
    \code{NaCl} \tab equilibrium constant for aqueous NaCl dissociation (Shock et al., 1992) \cr
    \code{cordierite} \tab equilibrium constant of hydrous cordierite dehydration \cr
    \code{phosphate} \tab phosphate speciation with pH, temperature and ionic strength \cr
    \code{nucleobase} \tab relative stabilities of nucleobases and some amino acids \cr
    \code{orp} \tab oxidation-reduction potential of redox standards as a function of temperature \cr
    \code{findit} \tab detailed example of usage of \code{\link{findit}} using log-normal distribution as an objective \cr
    \code{CO2Ac} \tab activity of CO2 buffered by acetic acid; comparing \code{\link{affinity}(return.buffer=TRUE)} with \code{\link{diagram}(what="CO2")} \cr
    \code{nonideal} \tab activity coefficient of charged species (Alberty, 2003), using the \code{IS} argument of \code{subcrt} \cr
  }

}

\references{
  Alberty, R. A. (2003) \emph{Thermodynamics of Biochemical Reactions}, John Wiley & Sons, Hoboken, New Jersey, 397 p. \url{http://www.worldcat.org/oclc/51242181}
  
  Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \eqn{^{\circ}}{°}C and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \url{http://dx.doi.org/10.1039/FT9928800803}
}


\examples{
\dontshow{data(thermo)}
demos(c("orp", "NaCl"))
\dontshow{par(thermo$opar)}
\dontrun{
# use the following to run examples in all help topics
examples()}
}


\keyword{extra}
