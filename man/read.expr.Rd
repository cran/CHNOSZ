\name{read.expr}
\alias{stress}
\alias{yeastgfp}
\alias{read.expr}

\title{Experimental Data for Protein Abundances and Localizations}

\description{
  Get abundance data from a protein expression experiment and add the proteins to the working instance of CHNOSZ.
  Retrieve the amino acid compositions of proteins with localizations and abundances taken from the YeastGFP project.
}

\usage{
  stress(condition, organism)
  yeastgfp(location, exclusive = TRUE)
  read.expr(file, idcol, abundcol, filter=NULL)
}

\arguments{
  \item{condition}{character, name of condition of stress response experiment}
  \item{organism}{character, organism in stress response experiment}
  \item{location}{character, name of subcellular location (compartment)}
  \item{exclusive}{logical, report only proteins exclusively localized to a compartment?}
  \item{file}{character, name of file with sequence IDs and abundance data}
  \item{idcol}{character, name of the column with sequence IDs}
  \item{abundcol}{character, name of the column with abundances}
  \item{filter}{list, optional filters to apply}
}

\details{

  \code{read.expr}, \code{yeastgfp} and \code{stress} all interact with data files stored in \code{\link{extdata}/abundance} to retrieve identities and possibly abundances of proteins in certain conditions.

  \code{stress} is the simplest of these functions since the source of its data, \code{stress.csv}, lists proteins without any abundance data. \code{condition} indicates the name of the stress response experiment (column name of \code{stress.csv}, e.g. \samp{low.C}) and \code{organism} indicates the organism (\samp{Eco} or \samp{Sce}).

The \code{yeastgfp} function returns the identities and abundances of proteins with the requested subcellular localization(s) (specified in \code{location}) using data from the YeastGFP project that is stored in \code{\link{extdata}/abundance/yeastgfp.csv.xz}.
The default value of \code{exclusive} (\code{FALSE}) tells the function to grab all proteins that are localized to a compartment even if they are also localized to other compartments.
If \code{exclusive} is \code{TRUE}, only those proteins that are localized exclusively to the requested compartments are identified, unless there are no such proteins, then the non-exclusive localizations are used (applies to the \samp{bud} localization).

  \code{read.expr} reads a \code{file} (CSV format) that contains protein sequence names or IDs and protein abundance data. \code{idcol} and \code{abundcol} are either the names of the columns holding the sequence IDs and protein abundances, or numeric values indicating the column numbers where these data are found. The column indicated by \code{abundcol} might not actually be abundance (it is likely to be abundance ratios). The data can be filtered to only include records that contain the term in the named argument \code{filter}, the name of which indicates the column to apply the filter to.

}

\value{
Each of these functions returns a list with elements named \code{protein} (names of proteins) and \code{abundance} (counts or concentrations without any conversion from the units in the data file).
For \code{stress}, the \code{abundance} value is all 1's.
For \code{yeastgfp}, if \code{location} is NULL, the function returns the names of all known locations, and if the length of \code{location} is >1, the \code{protein} and \code{abundance} values are lists of the results for each location.
}

\seealso{
  \code{\link{more.aa}} for getting the amino acid compositions of proteins whose names are returned by these functions.
}

\examples{
\dontshow{data(thermo)}

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

## read.expr using a different data file,
## from Anderson and Anderson, 2003
file <- system.file("extdata/abundance/AA03.csv", package="CHNOSZ")
# look for proteins described as "Complement"
read.expr(file, "name", "log10(pg/ml)", list(description="Complement"))

\donttest{## speciation diagram for ER.to.Golgi proteins (COPII coat 
## proteins) as a function of logfO2, after Dick, 2009
# add old parameters for [Met] sidechain to database
add.obigt()
y <- yeastgfp("ER.to.Golgi")
# don't use those with NA abundance
ina <- is.na(y$abundance)
# get the amino acid compositions of the proteins
aa <- more.aa(y$protein[!ina], "Sce")
# add proteins to thermo$protein
ip <- add.protein(aa)
# use logarithms of activities of proteins such
# that total activity of residues is unity
pl <- protein.length(ip)
logact <- unitize(rep(1, length(ip)), pl)
# load the proteins
basis("CHNOS+")
a <- affinity(O2=c(-80, -73), iprotein=ip, loga.protein=logact)
# make a speciation diagram
e <- equilibrate(a, normalize=TRUE)
diagram(e, ylim=c(-4.9, -2.9))
# where we are closest to experimental log activity
logfO2 <- rep(-78, length(ip))
abline(v=logfO2[1], lty=3)
# scale experimental abundances such that
# total activity of residues is unity
logact.expt <- unitize(log10(y$abundance[!ina]), pl)
# plot experimental log activity
points(logfO2, logact.expt, pch=16)
text(logfO2+0.5, logact.expt, y$protein[!ina])
# add title
title(main=paste("ER.to.Golgi; points - relative abundances",
  "from YeastGFP. Figure after Dick, 2009",sep="\n"))
# restore default thermodynamic database
data(thermo)

#############################
## examples using stress() ##
#############################

## predominance fields for overall protein compositions 
## induced and repressed in an/aerobic carbon limitation
## (experiments of Tai et al., 2005)
# the activities of ammonium and sulfate
# are similar to the non-growth-limiting concentrations
# used by Boer et al., 2003
basis(c("glucose", "H2O", "NH4+", "hydrogen", "SO4-2", "H+"),
  c(-1, 0, -1.3, 999, -1.4, -7))
# the names of the experiments in thermo$stress
expt <- c("Clim.aerobic.down", "Clim.aerobic.up",
  "Clim.anaerobic.down", "Clim.anaerobic.up")
# here we use abundance to indicate that the protein
# compositions should be summed together in equal amounts
for(i in 1:length(expt)) {
  p <- stress(expt[i], "Sce")
  aa <- more.aa(p$protein, "Sce")
  aa <- aasum(aa, average=TRUE, protein=expt[i])
  add.protein(aa)
}
species(expt, "Sce")
a <- affinity(C6H12O6=c(-30, 0), H2=c(-20, 0))
d <- diagram(a, normalize=TRUE, fill=NULL)
title(main=paste("Relative stabilities of proteins observed in\n",
  "an/aerobic carbon limitation in yeast"))
# the equilibrium distribution favors the proteins upregulated 
# by carbon limitation at low chemical potentials of C6H12O6 ...
stopifnot(c(d$predominant[1,1], d$predominant[1,128])==grep("up", expt))
# ... and favors proteins downregulated by aerobic conditions
# at high hydrogen fugacities
stopifnot(c(d$predominant[128, 128], d$predominant[128, 1])==grep("down", expt))

#################################
##  examples using yeastgfp()  ##
#################################

## Oxygen fugacity - activity of H2O predominance 
## diagrams for proteologs for 23 YeastGFP localizations
# need old properties of [Met] to reproduce this example
data(thermo)
add.obigt()
# arranged by decreasing metastability:
# order of this list of locations is based on the 
# (dis)appearance of species on the current set of diagrams
names <- c("vacuole", "early.Golgi", "ER", "lipid.particle",
  "cell.periphery", "ambiguous", "Golgi", "mitochondrion",
  "bud", "actin", "cytoplasm", "late.Golgi",
  "endosome", "nucleus", "vacuolar.membrane", "punctate.composite",
  "peroxisome", "ER.to.Golgi", "nucleolus", "spindle.pole",
  "nuclear.periphery", "bud.neck", "microtubule")
nloc <- c(4, 5, 3, 4, 4, 3)
# define the system
basis("CHNOS+")
# get protein names and abundances in each location
gfp <- yeastgfp(names)
# get amino acid compositions of proteins
aa <- more.aa(gfp$protein, "Sce")
# calculate average amino acid compositions 
for(i in 1:length(names)) {
  avgaa <- aasum(aa[[i]], gfp$abundance[[i]], average=TRUE, protein=names[i])
  add.protein(avgaa)
}
species(names, "Sce")
a <- affinity(H2O=c(-5, 0, 256), O2=c(-80, -66, 256))
# setup the plot
layout(matrix(c(1, 1,2:7), byrow=TRUE, nrow=4), heights=c(0.7, 3, 3, 3))
par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, paste("Subcellular proteins of S. cerevisiae,",
 "after Dick, 2009\n"), cex=1.5)
text(0.5, 0.2, describe.basis(ibasis=c(1, 3, 4, 6), oneline=TRUE), cex=1.5)
opar <- par(mar=c(3, 4, 1, 1), xpd=TRUE)
fill <- heat.colors(length(names))
inames <- 1:length(names)
for(i in 1:length(nloc)) {
  diagram(a, normalize=TRUE, names=names[inames], groups=as.list(inames),
    fill=fill[inames], cex.axis=0.75, cex.names=1)
  label.plot(letters[i])
  title(main=paste(length(inames), "locations"))
  # take out the stable species
  inames <- inames[-(1:nloc[i])]
}
# return to plot defaults
layout(matrix(1))
par(opar)
# reset thermodynamic database
data(thermo)

}
}

\references{
  Anderson, N. L. and Anderson, N. G. (2003) The human plasma proteome: History, character and diagnostic prospects (Vol. 1 (2002) 845-867). \emph{Molecular and Cellular Proteomics} \bold{2}, 50. \url{http://dx.doi.org/10.1074/mcp.A300001-MCP200}

  Boer, V. M., de Winde, J. H., Pronk, J. T. and Piper, M. D. W. (2003) The genome-wide transcriptional responses of \emph{Saccharomyces cerevisiae} grown on glucose in aerobic chemostat cultures limited for carbon, nitrogen, phosphorus, or sulfur. \emph{J. Biol. Chem.} \bold{278}, 3265--3274. \url{http://dx.doi.org/10.1074/jbc.M209759200}

  Dick, J. M. (2009) Calculation of the relative metastabilities of proteins in subcellular compartments of \emph{Saccharomyces cerevisiae}. \emph{BMC Syst. Biol.} \bold{3}:75. \url{http://dx.doi.org/10.1186/1752-0509-3-75}

  Ishihama, Y., Schmidt, T., Rappsilber, J., Mann, M., Hartl, F. U., Kerner, M. J. and Frishman, D. (2008) Protein abundance profiling of the \emph{Escherichia coli} cytosol. \emph{BMC Genomics} \bold{9}:102. \url{http://dx.doi.org/10.1186/1471-2164-9-102}

  Richmond, C. S., Glasner, J. D., Mau, R., Jin, H. F. and Blattner, F. R. (1999) Genome-wide expression profiling in \emph{Escherichia coli} K-12. \emph{Nucleic Acids Res.} \bold{27}, 3821--3835. \url{http://nar.oxfordjournals.org/cgi/content/abstract/27/19/3821}

  Tai, S. L., Boer, V. M., Daran-Lapujade, P., Walsh, M. C., de Winde, J. H., Daran, J.-M. and Pronk, J. T. (2005) Two-dimensional transcriptome analysis in chemostat cultures: Combinatorial effects of oxygen availability and macronutrient limitation in \emph{Saccharomyces cerevisiae}. \emph{J. Biol. Chem.} \bold{280}, 437--447. \url{http://dx.doi.org/10.1074/jbc.M410573200}
}

\keyword{protein}
