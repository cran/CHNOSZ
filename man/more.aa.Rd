\name{more.aa}
\alias{more.aa}
\title{Proteins from Model Organisms}
\description{
  Retrieve the amino acid compositions of one or more proteins from \emph{Escherichia coli} or \emph{Saccharomyces cerevisiae}.
}

\usage{
  more.aa(protein = NULL, organism)
}

\arguments{
  \item{protein}{character, name of protein}
  \item{organism}{character, name of organism (\samp{Eco} or \samp{Sce})}
}

\details{
\code{more.aa} retrieves the amino acid composition(s) of the indicated proteins in either \emph{Escherichia coli} or \emph{Saccharomyces cerevisiae}.
The value of \code{organism} can be one of \samp{Eco} or \samp{Sce}.
The calculation depends on the data files \code{\link{extdata}/protein/Eco.csv.xz} and \code{\link{extdata}/protein/Sce.csv.xz}, which contain the amino acid compositions of the proteins.
The \code{protein} argument should be a vector or a list of vectors of one or more Ordered Locus Names (OLN) or Open Reading Frame (ORF) names that are found in these files.
The output data frame contains rows with NA compositions for names that are not matched.
}

\value{
A data frame, or list of data frames, containing the amino acid composition(s) of the specified protein(s) in the format of \code{\link{thermo}$protein}.
}

\seealso{
\code{\link{extdata}} describes the sources of compositional data for the proteins. Other examples of usage of \code{more.aa} are shown for \code{\link{read.expr}}.
}


\examples{
\dontshow{data(thermo)}
# the first 13 names in UniProt for "aminotransferase ecoli"
at.ecoli <- c("BIOA", "ARNB", "SERC", "AAT", "TYRB", "ARGD", 
  "ILVE", "ALAA", "ALAC", "YBDL", "AVTA", "GLMS", "PUUE")
# get the amino acid compositions
# note that ALAA and ALAC are not matched
at.aa <- more.aa(at.ecoli, "Eco")
# what are their formulas and oxidation states
protein.formula(at.aa)
ZC(protein.formula(at.aa))
}

\concept{Protein thermodynamic modeling}
