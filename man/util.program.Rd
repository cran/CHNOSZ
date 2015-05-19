\name{util.program}
\alias{util.program}
\alias{caller.name}
\alias{msgout}
\alias{palply}

\title{Functions for Various Programming Tasks}
\description{
  Get name of calling function; send messages to standard output; use multiple processors if parallel package is loaded.
}

\usage{
  caller.name(n)
  msgout(..., domain = NULL, appendLF = FALSE)
  palply(X, FUN, ...)
}

\arguments{
  \item{n}{numeric, number of frame to go up}
  \item{...}{equivalent to the same argument in \code{\link{message}} or \code{\link[parallel]{parLapply}}}
  \item{domain}{equivalent to the same argument in \code{message}}
  \item{appendLF}{equivalent to the same argument in \code{message}}
  \item{X}{vector, argument for \code{\link{lapply}} or \code{parLapply}}
  \item{FUN}{function, argument for \code{lapply} or \code{parLapply}}
}

\details{

  \code{caller.name} returns the name of the calling function \code{n} frames up (i.e., for \code{n} equal to 2, the caller of the function that calls this one). If called interactively, returns \code{\link{character}()}.

  \code{msgout} is a variation of \code{\link{message}} in base R. It writes output to \code{stdout} instead of \code{stderr} and has a default setting of \code{appendLF} of FALSE. This function is used throughout CHNOSZ to generat informative messages, that will show up both in an interactive session and in \code{\link{Sweave}} output, but that are suppressed while running the \code{\link[testthat]{test_that}} testing scripts to make it easier to watch their progress.

 \code{palply} is a wrapper function to run \code{\link[parallel]{parLapply}} if length of \code{X} > 100 and package \pkg{parallel} is loaded, otherwise it runs \code{\link{lapply}}. Note that \code{parLapply} is called with \code{methods} set to FALSE. If lots of package startup messages are created when running \code{\link[parallel]{makeCluster}} (which is called by \code{palply}), it can probably be stopped by adding a test for \code{\link{interactive}} sessions around any \code{\link{library}} commands in the \code{\link{Rprofile}}.

}

\seealso{
  \code{palply} is used in \code{\link{read.fasta}} and \code{\link{grep.file}}, \code{\link{count.aa}} and \code{\link{protein.length}}, \code{\link{affinity}} when the \code{iprotein} argument is given, \code{\link{equil.boltzmann}} and \code{\link{equil.reaction}}, and in \code{\link{revisit}}.
}


\examples{
caller.name()  # character(0)
afun <- function() caller.name()
afun()         # character(0)
bfun <- function() afun()
bfun()         # "bfun"

msgout("h3llo w0rld\n")
}

\keyword{utilities}
