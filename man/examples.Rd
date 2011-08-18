\name{examples}
\alias{examples}
\alias{longex}
\title{Run examples from the documentation}
\description{
  Run the examples contained in each of the documentation topics. }

\arguments{
  \item{do.png}{logical, generate PNG files for the plots?}
  \item{which}{character, which example to run.}
}

\usage{
  examples(do.png = FALSE)
  longex(which)
}

\details{

  \code{examples} runs all the examples in the documentation for the package. The \code{\link{example}} function is called for each topic with \code{ask} set to \code{FALSE} (so all of the figures are shown without prompting the user). If \code{do.png} is TRUE, the plots in the examples are saved as \code{\link{png}} files having names beginning with the name of each of the help topics.

  \code{longex} is a function that contains code to run various other examples that are not contained in the documentation. The example to run is specified by \code{which}. The available examples are listed below. See the comments in the source code for more information about each one.

  \tabular{ll}{
    \code{sources} \tab cross-check the reference list with the thermodynamic database \cr
    \code{copper} \tab an Eh-pH diagram for the copper-water-glycine system \cr
    \code{cordierite} \tab equilibrium constant of hydrous cordierite dehydration \cr
    \code{phosphate} \tab phosphate speciation with pH, temperature and ionic strength \cr
    \code{nucleobase} \tab relative stabilities of nucleobases and some amino acids \cr
    \code{pie} \tab pie charts comparing relative abundances of organisms (Spear et al., 2005) and model proteins in metastable equilibrium \cr
    \code{orp} \tab oxidation-reduction potential of redox standards as a function of temperature \cr
    \code{findit} \tab detailed example of usage of \code{\link{findit}} using log-normal distribution as an objective \cr
    \code{co2ac} \tab activity of CO2 buffered by acetic acid; comparing \code{\link{affinity}(return.buffer=TRUE)} with \code{\link{diagram}(what="CO2")} \cr
  }

}

\references{
  Spear, J. R., Walker, J. J., McCollom, T. M. and  Pace, N. R. (2005) Hydrogen and bioenergetics in the Yellowstone geothermal ecosystem. \emph{Proc. Natl. Acad. Sci. U. S. A.} \bold{102}, 2555--2560.
}


\examples{
  \dontshow{data(thermo)}
  \donttest{
    longex(c("copper","orp"))
  }
}


\keyword{extra}
