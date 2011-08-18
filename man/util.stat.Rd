\name{util.stat}
\alias{util.stat}
\alias{spearman}
\alias{qqr}
\alias{rmsd}
\alias{cvrmsd}
\alias{lograt}
\title{Functions for Various Statistical Operations}
\description{
  Calculate root mean square deviation (RMSD), coefficient of variation of the RMSD, Spearman's rank correlation coefficient, or correlation coefficient of a q-q plot.
}

\usage{
  qqr(a)
  spearman(a, b)
  rmsd(a, b)
  cvrmsd(a, b)
  lograt(a, b)
}

\arguments{
  \item{a}{numeric values.}
  \item{b}{numeric values.}
}

\details{

  \code{qqr} returns the correlation coefficient of a quantile-quantile plot using the values given in \code{y}. 

  \code{spearman} calculates Spearman's rank correlation coefficient between \code{a} and \code{b}. \code{rmsd} and \code{cvrmsd} return the root mean square deviation and coefficient of variation of the RMSD, repsectively. The coefficient of variation is computed as the RMSD divided by the mean of the values in \code{a}

  \code{lograt} calculates the difference between a and b (i.e. log10 of the ratio of activities b/a) where \code{a} is a list of single values and \code{b} is a list of values, any dimension.

}

\seealso{
  These functions are among the target statistics available in \code{\link{revisit}}.
}

\examples{
  a <- 1:9
  b <- a + 1
  spearman(a, b)       # 1
  spearman(a, rev(b))  # -1
  qqr(a)               # <1
  rmsd(a,a)            # 0
  rmsd(a,b)            # 1
  cvrmsd(a,b)          # 1/5
}

\keyword{util}
