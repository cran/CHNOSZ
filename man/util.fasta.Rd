\name{util.fasta}
\alias{util.fasta}
\alias{grep.file}
\alias{read.fasta}
\alias{uniprot.aa}
\alias{count.aa}
\title{Functions for Reading FASTA Files and Downloading from UniProt}

\description{
  Search the header lines of a FASTA file, read protein sequences from a file, count numbers of amino acids in each sequence, and download sequences from UniProt.
}

\usage{
  grep.file(file, pattern = "", y = NULL, ignore.case = TRUE, 
    startswith = ">", lines = NULL, grep = "grep")
  read.fasta(file, i = NULL, ret = "count", lines = NULL, 
    ihead = NULL, start=NULL, stop=NULL, type="protein")
  count.aa(seq, start=NULL, stop=NULL, type="protein")
  uniprot.aa(protein, start=NULL, stop=NULL)
}

\arguments{
  \item{file}{character, path to FASTA file}
  \item{pattern}{character, pattern to search for in header lines}
  \item{y}{character, term to exclude in searching sequence headers}
  \item{ignore.case}{logical, ignore differences between upper- and lower-case?}
  \item{startswith}{character, only lines starting with this expression are matched}
  \item{lines}{list of character, supply the lines here instead of reading them from file}
  \item{grep}{character, name of system \code{grep} command}
  \item{i}{numeric, line numbers of sequence headers to read}
  \item{ret}{character, specification for type of return (count, sequence, or FASTA format)}
  \item{ihead}{numeric, which lines are headers}
  \item{start}{numeric, position in sequence to start counting}
  \item{stop}{numeric, position in sequence to stop counting}
  \item{type}{character, sequence type (protein or DNA)}
  \item{seq}{character, amino acid sequence of a protein}
  \item{protein}{character, entry name for protein in UniProt}
}

\details{
\code{grep.file} returns the line numbers of header lines in a FASTA file.
Matching header lines are identified having the search term \code{pattern} and optionally a term to exclude in \code{y}.
The \code{ignore.case} option is passed to \code{\link{grep}}, which does the work of finding lines that match.
Only lines that start with the expression in \code{startswith} are searched; the default setting reflects the format of the header lines in a FASTA file. 
If \code{y} is NULL and a supported operating system is identified, the operating system's \samp{grep} function (or other specified in the \code{grep} argument) is applied directly to the file instead of R's \code{\link{grep}}.
This avoids having to read the file into R using \code{\link{readLines}}.
If the lines from the file were obtained in a preceding operation, they can be supplied to this function in the \code{lines} argument.

\code{read.fasta} is used to retrieve entries from a FASTA file.
To read only selected sequences pass the line numbers of the header lines to the function in \code{i} (they can be identified using e.g. \code{grep.file}).
The function returns various formats depending on the value of \code{ret}; the default \samp{count} returns a dataframe of amino acid counts (the data frame can be given to \code{\link{add.protein}} in order to add the proteins to \code{\link{thermo}$protein}), \samp{seq} returns a list of sequences, and \samp{fas} returns a list of lines extracted from the FASTA file, including the headers (this can be used e.g. to generate a new FASTA file with only the selected sequences).
Similarly to \code{grep.file}, this function utilizes the OS's \samp{grep} on supported operating systems in order to identify the header lines as well as \samp{cat} to read the file, otherwise \code{\link{readLines}} and \R's \code{\link{substr}} are used to read the file and locate the header lines.
If the line numbers of the header lines were previously determined, they can be supplied in \code{ihead}. 
Optionally, the lines of a previously read file may be supplied in \code{lines} (in this case no file is needed so \code{file} should be set to "").

\code{count.aa} counts the occurrences of each amino acid or nucleic-acid base in a sequence (\code{seq}).
For amino acids, the columns in the returned data frame are in the same order as \code{thermo$protein}.
Letters are matched without regard for case.
A warning is generated if any character in \code{seq}, excluding spaces, is not one of the single-letter amino acid or nucleobase abbreviations.
\code{start} and/or \code{stop} can be provided to count a fragment of the sequence (extracted using \code{\link{substr}}).
If only one of \code{start} or \code{stop} is present, the other defaults to 1 (\code{start}) or the length of the sequence (\code{stop}).

\code{uniprot.aa} returns a data frame of amino acid composition, in the format of \code{thermo$protein}, retrieved from the protein sequence if it is available from UniProt (\url{http://uniprot.org}; The UniProt Consortium, 2012).
The \code{protein} argument corresponds to the \samp{Entry name} on the UniProt search pages.


}

\value{
\code{grep.file} returns a numeric vector.
\code{read.fasta} returns a list of sequences or lines (for \code{ret} equal to \samp{seq} or \samp{fas}, respectively), or a data frame with amino acid compositions of proteins (for \code{ret} equal to \samp{count}) with columns corresponding to those in \code{\link{thermo}$protein}. 
}

\seealso{
\code{\link{nucleic.formula}} for an example of counting nucleobases in a DNA sequence.
When computing relative abundances of many proteins that might be found with \code{grep.file} and \code{read.fasta}, consider using the \code{iprotein} arugment of \code{\link{affinity}} to speed things up; for an example see the help page for \code{\link{revisit}}. 
}

\examples{
## reading a protein FASTA file
# the path to the file
file <- system.file("extdata/fasta/EF-Tu.aln", package="CHNOSZ")
# read the sequences, and print the first one
read.fasta(file, ret="seq")[[1]]
# count the amino acids in the sequences
aa <- read.fasta(file)
# compute lengths (number of amino acids)
protein.length(aa)

\dontrun{
# download amino acid composition of a protein
# start at position 2 to remove the initiator methionine
aa <- uniprot.aa("ALAT1_HUMAN", start=2)
# add it to thermo$protein
ip <- add.protein(aa)
# now it's possible to calculate some properties
protein.length(ip)
protein.formula(ip)
subcrt("ALAT1_HUMAN", c("cr", "aq"), c(-1, 1))
# the amino acid composition can be saved for future use
write.csv(aa, "saved.aa.csv", row.names=FALSE)
# in another R session, the protein can be loaded without using uniprot.aa()
aa <- read.aa("saved.aa.csv")
add.protein(aa)

## count amino acids in a sequence
count.aa("GGSGG")
# warnings are issued for unrecognized characters
atest <- count.aa("WhatAmIMadeOf?")
# there are 3 "A" (alanine)
stopifnot(atest[, "A"]==3)
}
}

\references{
  The UniProt Consortium (2012) Reorganizing the protein space at the Universal Protein Resource (UniProt). \emph{Nucleic Acids Res.} \bold{40}, D71--D75. \url{http://dx.doi.org/10.1093/nar/gkr981}
}

\keyword{util}
