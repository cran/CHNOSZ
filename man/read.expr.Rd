\encoding{UTF-8}
\name{read.expr}
\alias{yeastgfp}
\alias{read.expr}

\title{Experimental Data for Protein Abundances and Localizations}

\description{
Get abundance data from a protein expression experiment and add the proteins to the current list of proteins.
Retrieve the amino acid compositions of proteins with localizations and abundances taken from the YeastGFP project.
}

\usage{
  yeastgfp(location, exclusive = TRUE)
  read.expr(file, idcol, abundcol, filter=NULL)
}

\arguments{
  \item{location}{character, name of subcellular location (compartment)}
  \item{exclusive}{logical, report only proteins exclusively localized to a compartment?}
  \item{file}{character, name of file with sequence IDs and abundance data}
  \item{idcol}{character, name of the column with sequence IDs}
  \item{abundcol}{character, name of the column with abundances}
  \item{filter}{list, optional filters to apply}
}

\details{
\code{read.expr} and \code{yeastgfp} read data files stored in \code{\link{extdata}/abundance} to retrieve identities and possibly abundances of proteins in certain conditions.

\code{yeastgfp} returns the identities and abundances of proteins with the requested subcellular localization(s) (specified in \code{location}) using data from the YeastGFP project that is stored in \code{\link{extdata}/abundance/yeastgfp.csv.xz}.
If \code{exclusive} is \code{FALSE}, the function grabs all proteins that are localized to a compartment even if they are also localized to other compartments.
If \code{exclusive} is \code{TRUE} (the default), only those proteins that are localized exclusively to the requested compartments are identified, unless there are no such proteins, then the non-exclusive localizations are used (applies to the \samp{bud} localization).

\code{read.expr} reads a \code{file} (CSV format) that contains protein sequence names or IDs and protein abundance data.
\code{idcol} and \code{abundcol} are either the names of the columns holding the sequence IDs and protein abundances, or numeric values indicating the column numbers where these data are found.
The column indicated by \code{abundcol} might not actually be abundance (it is likely to be abundance ratios).
The data can be filtered to only include records that contain the term in the named argument \code{filter}, the name of which indicates the column to apply the filter to.
}

\value{
Each of these functions returns a list with elements named \code{protein} (names of proteins) and \code{abundance} (counts or concentrations without any conversion from the units in the data file).
For \code{yeastgfp}, if \code{location} is NULL, the function returns the names of all known locations, and if the length of \code{location} is >1, the \code{protein} and \code{abundance} values are lists of the results for each location.
}

\seealso{ \code{\link{more.aa}} for getting the amino acid compositions of the proteins. }

\examples{\dontshow{data(thermo)}
## overall oxidation state of proteins exclusively localized 
## to cytoplasm of S. cerevisiae with/without abundance weighting
y <- yeastgfp("cytoplasm")
aa <- more.aa(y$protein, "Sce")
aaavg <- aasum(aa, average=TRUE)
ZC(protein.formula(aaavg))
# the average composition weighted by abundance
waaavg <- aasum(aa, abundance=y$abundance, average=TRUE)
ZC(protein.formula(waaavg))

## read.expr using one of the provided data files,
## from Ishihama et al., 2008
file <- system.file("extdata/abundance/ISR+08.csv.xz", package="CHNOSZ")
# read all protein names and abundances in ID and emPAI columns
# (emPAI - exponentially modified protein abundance index)
expr <- read.expr(file, "ID", "emPAI")
# scatter plot of average oxidation state and emPAI
aa <- more.aa(expr$protein, "Eco")
pf <- protein.formula(aa)
zc <- ZC(pf)
# note we specify ylim here that excludes some high-emPAI values
plot(zc, expr$abundance, xlab=expr.property("ZC"), ylim=c(0, 90), ylab="emPAI",
  main="Proteins in E. coli cytosol\nAbundance vs oxidation state of carbon")
legend("topleft", pch=1, legend="Ishihama et al., 2008")
# what if we just want kinases?
# "description" is the name of the column where we search for "kinase"
expr.kinase <- read.expr(file, "ID", "emPAI", list(description="kinase"))

## potential fields for overall protein compositions 
## transcriptionally induced and repressed in aerobic
## and anaerobic carbon limitation
## (experiments of Tai et al., 2005)
# the activities of ammonium and sulfate used here
# are similar to the non-growth-limiting concentrations
# used by Boer et al., 2003
basis(c("glucose", "H2O", "NH4+", "hydrogen", "SO4-2", "H+"),
  c(-1, 0, -1.3, 999, -1.4, -7))
# the names of the experiments in TBD+05.csv
expt <- c("Clim.aerobic.down", "Clim.aerobic.up",
  "Clim.anaerobic.down", "Clim.anaerobic.up")
file <- system.file("extdata/abundance/TBD+05.csv", package="CHNOSZ")
dat <- read.csv(file, as.is=TRUE)
# more.aa: get the amino acid compositions
# aasum: average them together
for(thisexpt in expt) {
  p <- dat$protein[dat[, thisexpt]]
  aa <- more.aa(p, "Sce")
  aa <- aasum(aa, average=TRUE, protein=thisexpt)
  add.protein(aa)
}
species(expt, "Sce")
a <- affinity(C6H12O6=c(-30, 0), H2=c(-20, 0))
d <- diagram(a, normalize=TRUE, fill=NULL)
title(main=paste("Formation potential of proteins associated with\n",
  "transcriptional response to carbon limitation in yeast"))
# the affinity of formation favors the proteins upregulated 
# by carbon limitation at low chemical potentials of C6H12O6 ...
stopifnot(c(d$predominant[1,1], d$predominant[1,128])==grep("up", expt))
# ... and favors proteins downregulated by aerobic conditions
# at high hydrogen fugacities
stopifnot(c(d$predominant[128, 128], d$predominant[128, 1])==grep("down", expt))
}

\references{
Boer, V. M., de Winde, J. H., Pronk, J. T. and Piper, M. D. W. (2003) The genome-wide transcriptional responses of \emph{Saccharomyces cerevisiae} grown on glucose in aerobic chemostat cultures limited for carbon, nitrogen, phosphorus, or sulfur. \emph{J. Biol. Chem.} \bold{278}, 3265--3274. \url{https://doi.org/10.1074/jbc.M209759200}

Ishihama, Y., Schmidt, T., Rappsilber, J., Mann, M., Hartl, F. U., Kerner, M. J. and Frishman, D. (2008) Protein abundance profiling of the \emph{Escherichia coli} cytosol. \emph{BMC Genomics} \bold{9}:102. \url{https://doi.org/10.1186/1471-2164-9-102}

Tai, S. L., Boer, V. M., Daran-Lapujade, P., Walsh, M. C., de Winde, J. H., Daran, J.-M. and Pronk, J. T. (2005) Two-dimensional transcriptome analysis in chemostat cultures: Combinatorial effects of oxygen availability and macronutrient limitation in \emph{Saccharomyces cerevisiae}. \emph{J. Biol. Chem.} \bold{280}, 437--447. \url{https://doi.org/10.1074/jbc.M410573200}
}

\concept{Protein thermodynamic modeling}
