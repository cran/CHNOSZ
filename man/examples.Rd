\name{examples}
\alias{examples}
\alias{longex}
\title{Run examples from the documentation}
\description{
  Run the examples contained each of the documentation topics. }

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
    \code{pie} \tab pie charts for relative abundances of organisms and model proteins \cr
  }

}

\keyword{misc}

\examples{
  \dontshow{data(thermo)}
  \donttest{
    longex("cordierite")
  }
}
