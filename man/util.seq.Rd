\name{util.seq}
\alias{util.seq}
\alias{aminoacids}
\alias{nucleicacids}

\title{Functions to Work with Sequence Data}

\description{
  Count amino acids in protein sequences, return one- or three-letter abbreviations of amino acids; count nucleotides in nucleic acid sequences, calculate DNA and RNA complements of nucleic acid sequences.
}

\usage{
  aminoacids(seq, nchar=1)
  nucleicacids(seq, type = "DNA", comp = NULL, comp2 = NULL)
}

\arguments{
  \item{seq}{character, amino acid sequence of a protein (\code{aminoacids}) or base sequence of a nucleic acid (\code{nucleicacids}).}
  \item{nchar}{numeric, \eqn{1} to return one-letter, \eqn{3} to return three-letter abbreviations for amino acids (\code{aminoacids}).}
  \item{type}{character, type of nucleic acid sequence (DNA or RNA) (\code{nucleicads}).}
  \item{comp}{character, type of complement sequence.}
  \item{comp2}{character, type of second complement sequence.}
}

\details{

  \code{aminoacids} takes a character argument containing a protein sequence and counts the number of occurrences of each type of amino acid. The output is a dataframe with 20 columns, each corresponding to an amino acid, ordered in the same way as \code{thermo$protein}. If the first argument is NULL, the function returns the one-letter abbreviations (for \code{nchar} equal to 1) or the three-letter ones (if \code{nchar} is equal to 3) or the names of the amino acids (if \code{nchar} is NA) of twenty amino acids in the order used in \code{thermo$protein}.

  \code{nucleicacids} takes a DNA or RNA sequence and counts the numbers of bases of each type. Whether the sequence is DNA or RNA is specified by \code{type}. Setting \code{comp} to \samp{DNA} or \samp{RNA} tells the function to compute the base composition of that type of complement of the sequence. If \code{comp2} is specified, another complement is taken. The two rounds of complementing can be used in a single function call e.g. to go from a sequence on DNA minus strand (given in \code{seq}) to the plus strand (with \code{comp="DNA"}) and then from the DNA plus strand to RNA (with \code{comp2="RNA"}). The value returned by the function is a dataframe of base composition, which can be passed back to the function to obtain the overall chemical formula for the bases.

}

\value{
  An object of type character or dataframe.
}

\examples{
  \dontshow{data(thermo)}
  ## count amino acids in a sequence
  aminoacids("GGSGG")
  aminoacids("WhatAmIMadeOf?")

  ## count nucleobases in a sequence
  nucleicacids("ACCGGGTTT")
  # the DNA complement of that sequence
  nucleicacids("ACCGGGTTT",comp="DNA")
  # the RNA complement of the DNA complement
  n <- nucleicacids("ACCGGGTTT",comp="DNA",comp2="RNA")
  # the formula of the RNA complement
  nucleicacids(n,type="RNA")
}

\keyword{misc}



