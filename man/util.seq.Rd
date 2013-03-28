\name{util.seq}
\alias{util.seq}
\alias{aminoacids}

\alias{nucleicacids}

\title{Functions to Work with Sequence Data}

\description{
  Return one- or three-letter abbreviations of amino acids; count nucleotides in nucleic acid sequences, calculate DNA and RNA complements of nucleic acid sequences.
}

\usage{
  aminoacids(nchar=1, which=NULL)
  nucleicacids(seq, type = "DNA", comp = NULL, comp2 = NULL)
}

\arguments{
  \item{nchar}{numeric, \eqn{1} to return one-letter, \eqn{3} to return three-letter abbreviations for amino acids}
  \item{which}{character, which amino acids to name}
  \item{seq}{character, base sequence of a nucleic acid (\code{nucleicacids})}
  \item{type}{character, type of nucleic acid sequence (DNA or RNA) (\code{nucleicads})}
  \item{comp}{character, type of complement sequence}
  \item{comp2}{character, type of second complement sequence}
}

\details{

  Depending on the value of \code{nchar}, \code{aminoacids} returns the one-letter abbreviations (\samp{1}) or the three-letter abbreviations (\samp{3}) or the names of the neutral amino acids (\samp{""}) or the names of the amino acids with ionized side chains (\samp{"Z"}). The default includes 20 amino acids in the order used in \code{thermo$protein}, unless \code{which} is provided, indicating the desired amino acids (either as 1- or 3-letter abbreviations or names of the neutral amino acids).

  \code{aminoacids} takes a character value or list of character values containing a protein sequence. The number of occurrences of each type of amino acid is counted and output in a data frame with 20 columns, ordered in the same way as \code{thermo$protein}.

  \code{nucleicacids} takes a DNA or RNA sequence and counts the numbers of bases of each type. Whether the sequence is DNA or RNA is specified by \code{type}. Setting \code{comp} to \samp{DNA} or \samp{RNA} tells the function to compute the base composition of that type of complement of the sequence. If \code{comp2} is specified, another complement is taken. The two rounds of complementing can be used in a single function call e.g. to go from a sequence on DNA minus strand (given in \code{seq}) to the plus strand (with \code{comp="DNA"}) and then from the DNA plus strand to RNA (with \code{comp2="RNA"}). The value returned by the function is a dataframe of base composition, which can be passed back to the function to obtain the overall chemical formula for the bases.

}

\examples{\dontshow{data(thermo)}
## count nucleobases in a sequence
nucleicacids("ACCGGGTTT")
# the DNA complement of that sequence
nucleicacids("ACCGGGTTT", comp="DNA")
# the RNA complement of the DNA complement
n <- nucleicacids("ACCGGGTTT", comp="DNA", comp2="RNA")
# the formula of the RNA complement
(nf <- nucleicacids(n, type="RNA"))
stopifnot(identical(nf, "C40H42N32O11"))
}

\keyword{util}
