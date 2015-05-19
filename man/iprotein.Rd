\name{iprotein}
\alias{iprotein}
\alias{ip2aa}
\alias{aa2eos}
\alias{seq2aa}
\alias{aasum}
\alias{read.aa}
\alias{add.protein}
\title{Amino Acid Compositions of Proteins}
\description{
  Functions to identify proteins, get and set amino acid compositions, and calculate thermodynamic properties from group additivity.
}

\usage{
  iprotein(protein, organism=NULL)
  ip2aa(protein, organism=NULL, residue=FALSE)
  aa2eos(aa, state=get("thermo")$opt$state)
  seq2aa(protein, sequence)
  aasum(aa, abundance = 1, average = FALSE, protein = NULL, organism = NULL)
  read.aa(file = "protein.csv")
  add.protein(aa)
}

\arguments{
  \item{protein}{character, name of protein; numeric, indices of proteins (rownumbers of \code{\link{thermo}$protein})}
  \item{organism}{character, name of organism}
  \item{residue}{logical, compute per-residue counts?}
  \item{aa}{data frame, amino acid composition in the format of \code{thermo$protein}}
  \item{state}{character, physical state}
  \item{sequence}{character, protein sequence}
  \item{abundance}{numeric, abundances of proteins}
  \item{average}{logical, return the weighted average of amino acid counts?}
  \item{file}{character, path to file with amino acid compositions}
}

\details{
  A \samp{protein} in CHNOSZ is defined by a name and by the counts of amino acids, stored in \code{\link{thermo}$protein}. The purpose of the functions described here is to identify proteins and work with their amino acid compositions. From the amino acid compositions, the thermodynamic properties of the proteins can be estimated (Dick et al., 2006) for use in other functions in the package. 

  Often, the names of proteins are sufficient to set up calculations using functions such as \code{\link{subcrt}} or \code{\link{species}}. The names of proteins in CHNOSZ are distinguished from those of other chemical species by having an underscore character ("_") that separates two identifiers, referred to as the \code{protein} and \code{organism} (but any other meaning can be attached to these names). An example is \samp{LYSC_CHICK}. 

  The first few functions provide low-level operations:

  \code{iprotein} returns the rownumber(s) of \code{thermo$protein} that match the protein names. The names can be supplied in the single \code{protein} argument or as separated \code{protein}s and \code{organism}s (without the underscore). Any protein not matched returns an NA and generates a message.

  \code{ip2aa} returns the row(s) of \code{thermo$protein} that match the supplied protein names, OR the protein indices found by \code{iprotin}. Set \code{residue} to TRUE to return the per-residue composition (i.e. amino acid composition of the protein divided by total number of residues). For this function only, if the \code{protein} argument is a data frame, it is returned unchanged, except for possibly the per-residue calculation.

  \code{aa2eos} calculates the thermodynamic properties and equations-of-state parameters for the completely nonionized proteins using group additivity with parameters taken from Dick et al., 2006 (aqueous proteins) and LaRowe and Dick, 2012 (crystalline proteins and revised aqueous methionine sidechain group). The return value is a data frame in the same format as \code{thermo$obigt}. \code{state} indicates the physical state for the parameters used in the calculation (\samp{aq} or \samp{cr}).

  The remaining functions are more likely to be called directly by the user:

  \code{seq2aa} returns a data frame of amino acid composition, in the format of \code{thermo$protein}, corresponding to the provided \code{sequence}. Here, the \code{protein} argument indicates the name of the protein with an underscore (e.g. \samp{LYSC_CHICK}).

  \code{aasum} returns a data frame representing the sum of amino acid compositions in the rows of the input \code{aa} data frame. The amino acid compositions are multiplied by the indicated \code{abundance}; that argument is recycled to match the number of rows of \code{aa}. If \code{average} is TRUE the final sum is divided by the number of input compositions. The name used in the output is taken from the first row of \code{aa} or from \code{protein} and \code{organism} if they are specified. This function is useful for calculating bulk amino acid compositions in stress response experiments or localization studies; see \code{\link{read.expr}} for examples of its use.

\code{read.aa} returns a data frame of amino acid composition based on the contents of the indicated \code{file}, which should be a CSV file with the same column names as \code{thermo$protein}.

\code{add.protein} completes the loop; any amino acid composition returned by the \code{*aa} functions described above can be added to \code{thermo$protein} using this function to be made available to other functions in the package.
The amino acid compositions of proteins in \code{aa} with the same name as one in \code{thermo$protein} are replaced.
The value returned by this function is the rownumbers of \code{thermo$protein} that are added and/or replaced.

}

\examples{
\dontshow{data(thermo)}
# search by name in thermo$protein
ip1 <- iprotein("LYSC_CHICK")
ip2 <- iprotein("LYSC", "CHICK")
# these are the same
stopifnot(all.equal(ip1, ip2))
# two organisms with the same protein name
ip3 <- iprotein("MYG", c("HORSE", "PHYCA"))
# their amino acid compositions
ip2aa(ip3)
# their thermodynamic properties by group additivity
aa2eos(ip2aa(ip3))

# an example of an unrecognized protein name
ip4 <- iprotein("MYGPHYCA")
stopifnot(is.na(ip4))

# manually adding a new protein
# Human Gastric juice peptide 1
aa <- seq2aa("GAJU_HUMAN", "LAAGKVEDSD")
ip <- add.protein(aa)
stopifnot(protein.length(ip)==10)
# the chemical formula of this peptide
stopifnot(as.chemical.formula(protein.formula(ip))=="C41H69N11O18")
# we can also calculate a formula without using add.protein
as.chemical.formula(protein.formula(seq2aa("pentapeptide_test", "ANLSG")))

# read a fasta file, calculate H/C and O/C for all proteins
# and averages by polypeptide chain, residue and mass
ffile <- system.file("extdata/fasta/HTCC1062.faa.xz", package="CHNOSZ")
aa <- read.fasta(ffile)
pf <- as.data.frame(protein.formula(aa))
plot(pf$H/pf$C, pf$O/pf$C, pch=NA)
points(pf$H/pf$C, pf$O/pf$C, pch=20, cex=0.5, col="grey")
# average composition per polypeptide chain
chainaa <- aasum(aa, average=TRUE)
chainpf <- as.data.frame(protein.formula(chainaa))
points(chainpf$H/chainpf$C, chainpf$O/chainpf$C, pch=15, col="red")
# average by amino acid residue
pl <- protein.length(aa)
resaa <- aasum(aa, abundance=pl, average=TRUE)
respf <- as.data.frame(protein.formula(resaa))
points(respf$H/respf$C, respf$O/respf$C, pch=16, col="red")
# average by mass
pm <- mass(pf)
massaa <- aasum(aa, abundance=pm, average=TRUE)
masspf <- as.data.frame(protein.formula(massaa))
points(masspf$H/masspf$C, masspf$O/masspf$C, pch=17, col="red")
# add a legend and title
legend("topright", pch=c(20, 15, 16, 17), col=c("grey", rep("red", 3)),
  legend=c("protein", "chain average", "residue average", "mass average"))
title(main=paste("O/C vs H/C for HTCC1062 proteins\n",
  "n =", nrow(aa)))
}

\seealso{
\code{\link{read.fasta}} and \code{\link{uniprot.aa}} for getting amino acid compositions from a FASTA file or downloading them from UniProt.

\code{\link{protein.info}} for higher-level functions (chemical formulas, summaries of reaction coefficients and energies), \code{\link{more.aa}} for getting amino acid compositions for model organisms from additional data files in the \code{\link{extdata}/protein} directory, and \code{\link{read.expr}} for working with protein abundance and subcellular localization data.

Examples of stability calculations for proteins are in \code{\link{protein}}.
}

\references{
  Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \url{http://dx.doi.org/10.5194/bg-3-311-2006}

}

\concept{Protein thermodynamic modeling}
