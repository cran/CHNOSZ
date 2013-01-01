\name{util.character}
\alias{util.character}
\alias{c2s}
\alias{s2c}
\alias{can.be.numeric}
\title{Functions to Manipulate Character Objects}

\description{
Convert between strings and character objects. Test for ability to become numeric.
}

\usage{
  c2s(x, sep = " ")
  s2c(x, sep = NULL, keep.sep = TRUE)
  can.be.numeric(x)
}

\arguments{
  \item{x}{character object to convert (\code{s2c}, \code{c2s}, \code{axis.label}), or object to be tested (\code{can.be.numeric}).}
  \item{sep}{character, the separator to insert or separator(s) to match (\code{c2s}, \code{s2c}).}
  \item{keep.sep}{logical, retain the separator in the output (TRUE) or discard it (FALSE) (\code{s2c}).}
}

\details{
  \code{c2s} joins the elements of a character object into a character object of length \eqn{1} (a string), and \code{s2c} splits a string into elements of a character object of length \eqn{n+1}, where \eqn{n} stands for the number of separators in the string. \code{sep} gives the separator to insert between successive items (in \code{c2s}) or the separator(s) to find in a string (in \code{s2c}). The default value of \code{sep} is a space (" ") in \code{c2s}. The default value for \code{sep} is NULL in \code{s2c}, indicating a separator at every position of \code{x} (the result in this case has length equal to \code{nchar(x)}). Argument \code{keep.sep} if \code{TRUE} (the default) instructs \code{s2c} to keep the separating values in the output. The maximum length of the object returned by \code{s2c} is determined by \code{n}; the default value of \code{NULL} indicates an unrestricted length.

  \code{can.be.numeric} returns a value of TRUE or FALSE for each element of \code{x}.
}

\value{
\code{s2c}, \code{c2s} and \code{axis.label} return character values. \code{can.be.numeric} returns logical.
}

\examples{
  \dontshow{data(thermo)}
  ## string to character
  s2c("hello world")
  s2c("hello world",sep=" ",keep.sep=FALSE)
  s2c("3.141592",sep=c(".","9"))
  # character to string
  c2s(aminoacids())
  c2s(aminoacids(),sep=".")
}

\keyword{util}
