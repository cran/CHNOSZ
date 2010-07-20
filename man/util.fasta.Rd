\name{util.fasta}
\alias{util.fasta}
\alias{is.fasta}
\alias{grep.file}
\alias{read.fasta}
\alias{splitline}
\alias{trimfas}
\alias{sum.refseq}
\title{Functions for Accessing FASTA Files}

\description{
  Search the header lines of a FASTA file, read protein sequences from a file and count numbers of amino acids in each sequence.
}

\usage{
  is.fasta(file)
  grep.file(file, x = "", y = NULL, ignore.case = TRUE, startswith = ">", 
    lines = NULL, grep = "grep")
  read.fasta(file, i = NULL, ret = "count", lines = NULL, ihead = NULL)
  splitline(line, length)
  trimfas(file, start, stop)
}

\arguments{
  \item{file}{character, path to FASTA file.}
  \item{x}{character, pattern to search for in header lines.}
  \item{y}{character, term to exclude in searching sequence headers.}
  \item{ignore.case}{logical, ignore differences between upper- and lower-case?}
  \item{startswith}{character, only lines starting with this expression are matched.}
  \item{lines}{list of character, supply the lines here instead of reading them from file.}
  \item{grep}{character, name of system \code{grep} command.}
  \item{i}{numeric, line numbers of sequence headers to read.}
  \item{ret}{character, specification for type of return (count, sequence, or FASTA format).}
  \item{ihead}{numeric, which lines are headers.}
  \item{line}{character, a line to be split into multiple lines.}
  \item{length}{numeric, the maximum length of any line.}
  \item{start}{numeric, starting position to extract from sequences.}
  \item{stop}{numeric, last position to extract from sequences.}
}

\details{
  \code{is.fasta} checks if a file is in FASTA format. A very simple test is performed: if either of the first two lines of the file starts with \samp{>}, then the function returns TRUE, otherwise it returns FALSE.

  \code{grep.file} is used to search for entries in a FASTA file. It returns the line numbers of the matching FASTA headers. It takes a search term in \code{x} and optionally a term to exclude in \code{y}. The \code{ignore.case} option is passed to \code{\link{grep}}, which does the work of finding lines that match. Only lines that start with the expression in \code{startswith} are searched; the default setting reflects the format of the header line for each sequence in a FASTA file. If \code{y} is NULL and a supported operating system is identified (right now, only Linux), the operating system's \samp{grep} function (or other specified in the \code{grep} argument) is applied directly to the file instead of R's \code{\link{grep}}; this mitigates the potential speed penalty of having to read the file into R using \code{\link{readLines}}. If the lines from the file were obtained in a preceding operation, they can be supplied to this function in the \code{lines} argument (this overrides the use of the OS's \samp{egrep}).

  \code{read.fasta} is used to retrieve entries from a FASTA file. The line numbers for the headers of the desired sequences are passed to the function in \code{i} (they can be generated using \code{grep.file}). The function returns various formats depending on the value of \code{ret}; the default \samp{count} returns a dataframe of amino acid counts (the dataframe can be given to \code{\link{add.protein}} in order to add the proteins to \code{\link{thermo}$protein}), \samp{seq} returns a list of sequences, and \samp{fas} returns a list of lines extracted from the FASTA file, including the headers (this can be used e.g. to generate a new FASTA file with only the selected sequences). Similarly to \code{grep.file}, this function utilizes the OS's \samp{grep} on supported operating systems in order to identify the header lines as well as \samp{cat} to read the file, otherwise \code{\link{readLines}} and \R's \code{\link{substr}} are used to read the file and locate the header lines. \code{lines}, if it is given, bypasses the reading of the file and also overrides the use of the OS's tools. If the line numbers of the header lines were previously determined, they can be supplied in \code{ihead}. When computing relative abundances of many proteins that might be found with \code{grep.fasta}, consider using the \code{iprotein} arugment of \code{\link{affinity}} to speed things up. Examples of these operations can be found in the documentation for \code{\link{revisit}}.

  \code{splitline} takes a single character object (the \code{line}) and splits it into multiple lines of the given length (the last line can be shorter than this). It returns a character object that contains the lines. This function is utilized by \code{trimfas}, which extracts the specified positions from a (usually) aligned FASTA file. The length of the lines used by \code{trimfas} is equal to the length of the first sequence line in the given \code{file}.
}

\value{
  \code{grep.file} returns a numeric vector. \code{read.fasta} returns a list of sequences or lines (for \code{ret} equal to \samp{seq} or \samp{fas}, respectively), or a data frame with amino acid compositions of proteins (for \code{ret} equal to \samp{count}) with columns corresponding to those in \code{\link{thermo}$protein}. 
}

\seealso{
  For examples of using \code{grep.file} and \code{read.fasta}, see the help page for \code{\link{revisit}}.
}

\section{Side Effects}{None}

\examples{
  # basic use of splitline
  splitline("abcdefghijklmnopqrstuvwxyz",10)

  # get the first ten positions of each of the sequences
  # in a FASTA file
  f <- system.file("HTCC1062.faa",package="CHNOSZ")
  fnew <- trimfas(f,1,10)
  # just show the first few entries
}

\keyword{misc}
