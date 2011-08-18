\name{get.expr}
\alias{get.expr}
\title{Protein Expression Data}
\description{
Get abundance data from a protein expression experiment and add the proteins to the working instance of CHNOSZ.}

\usage{
  get.expr(file, idcol, abundcol, seqfile, filter=NULL, 
    is.log=FALSE, loga.total = 0)
}

\arguments{
  \item{file}{character, name of file with sequence IDs and abundance data.}
  \item{idcol}{character, name of the column with sequence IDs.}
  \item{abundcol}{character, name of the column with abundances.}
  \item{seqfile}{character, name of the FASTA file with protein sequences.}
  \item{filter}{list, optional filters to apply.}
  \item{is.log}{logical, are the abundances in the file in logarithmic (base 10) units?}
  \item{loga.total}{numeric, logarithm of total activity of residues.}
}

\details{

  This function reads a CSV \code{file} that contains protein sequence IDs and protein abundance data. The header (first line) of this file contains the column names; the names of the columns holding the sequence IDs and protein abundances are indicated by \code{idcol} and \code{abundcol}, respectively. The sequence IDs are searched for in the accession lines in the FASTA file indicated by \code{seqfile} (using \code{\link{grep}}); a match can occur in any part of an accession line, and the first such match is used. Any IDs that are NA or can not be found in \code{seqfile} are excluded from further consideration. The amino acid compositions of the matched proteins are computed (using \code{\link{read.fasta}}) and are added to the inventory of proteins in CHNOSZ (\code{\link{thermo}$protein}). 

  The function returns values of the logarithms of activities of the proteins. We associate molality with activity (i.e., activity coefficients are implicitly unity).  If \code{loga.total} is not NULL, the abundances of the proteins from the data file are scaled to give a logarithm of total activity of amino acid residues equal to the value in \code{loga.total}, usually set to zero (see \code{\link{unitize}}). This operation preserves the relative abundances of the proteins. If the abundances of the proteins in the file are already in logarithmic units, set \code{is.log} to TRUE.

  If \code{seqfile} is one of \samp{SGD}, \samp{ECO} or \samp{HUM} it refers to the database of amino acid compositions of proteins packaged with CHNOSZ for either \emph{Saccharomyces cerevisiae}, \emph{Escherichia coli} or \emph{Homo sapiens}. In this case, the search for matching IDs is performed using \code{\link{get.protein}}.

  The data file can be filtered by using \code{filter}. This argument should be a list with one element, the name of which indicates the column to apply the filter to, and the value of which is a search term. 

}

\value{
  Returns a list with objects \code{iprotein} (the indices of the proteins in \code{\link{thermo}$protein}) and \code{loga.ref} (the logarithms of activities of the proteins).
}

\seealso{
  \code{\link{findit}} for finding combinations of chemical activities that optimize the fit of metastable protein assemblages to experimental protein abundances.
}

\examples{
  \dontshow{data(thermo)}
  # let's use a sample data file
  file <- system.file("extdata/abundance/ISR+08.csv",package="CHNOSZ")
  # read the abundances and get the proteins from ECO.csv
  expr <- get.expr(file,"ID","emPAI","ECO")
  # what if we just wanted kinases?
  expr <- get.expr(file,"ID","emPAI","ECO",list(description="kinase"))
  # the abundances were scaled so that the total activity of residues is unity
  pl <- protein.length(-expr$iprotein)
  stopifnot(all.equal(sum(pl*10^expr$loga),1))
  # see the 'protactiv' vignette for comparison with equilibrium calculations

  # if you want to read the protein sequences from a FASTA file...
  # e <- get.expr(file,"ID","emPAI","ECOLI.fasta")
}

\keyword{primary}
