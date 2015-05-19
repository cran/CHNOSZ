\name{util.list}
\alias{util.list}
\alias{lsub}
\alias{lsum}
\alias{pprod}
\alias{which.pmax}
\title{Functions to Work with Lists}

\description{
  Combine lists or perform arithmetic operations on elements of lists.
}

\usage{
  lsub(x, y)
  lsum(x, y)
  pprod(x, y)
  which.pmax(elts, na.rm = FALSE, pmin = FALSE)
}

\arguments{
  \item{x}{list}
  \item{y}{list (\code{lsub}, \code{lsum}), or numeric (\code{pprod})}
  \item{elts}{list, numeric vectors for which to find maximum values (in parallel) (\code{which.pmax}).}
  \item{na.rm}{logical, remove missing values?}
  \item{pmin}{logical, find minimum values instead of maximum ones?}
}

\details{
  \code{lsub} subtracts the elements of list \code{y} from the respective ones in list \code{x}. \code{lsum} sums the respective elements of lists \code{x} and \code{y}. \code{pprod} multiplies each element of list \code{x} by the respective numeric value in \code{y}.

  \code{which.pmax} takes a list of equal-length numeric vectors (or objects that can be coerced to numeric) in \code{elts} and returns the index of the vector holding the maximum value at each position. If \code{na.rm} is \code{TRUE}, values of \code{NA} are removed; if \code{pmin} is \code{TRUE} the function finds locations of the minimum values instead.


}

\value{
  \code{lsub}, \code{lsum} and \code{pprod} return lists. 
}

\keyword{utilities}
