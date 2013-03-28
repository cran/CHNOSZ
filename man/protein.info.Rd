\encoding{UTF-8}
\name{protein.info}
\alias{group.formulas}
\alias{protein.formula}
\alias{protein.length}
\alias{protein.info}
\alias{protein.basis}
\alias{protein.equil}
\alias{MP90.cp}
\title{Summaries of Thermodynamic Properties of Proteins}

\description{
  Calculate chemical formulas, lengths, standard Gibbs energies and net charges, stoichiometric coefficients of basis species in reactions to form proteins (possibly per residue), and show steps in calculation of chemical activities of proteins in metastable equilibrium.
}

\usage{
  protein.formula(protein, organism = NULL, residue = FALSE)
  protein.length(protein, organism = NULL)
  protein.info(protein, T = 25, residue = FALSE, round.it = FALSE)
  protein.basis(protein, T = 25, normalize = FALSE)
  protein.equil(protein, T=25, loga.protein = 0)
  MP90.cp(protein, T)
  group.formulas()
}

\arguments{
  \item{protein}{character, names of proteins; numeric, species index of proteins; data frame; amino acid composition of proteins}
  \item{organism}{character, names of organisms}
  \item{residue}{logical, return per-residue values (those of the proteins divided by their lengths)?}
  \item{normalize}{logical, return per-residue values (those of the proteins divided by their lengths)?}
  \item{T}{numeric, temperature in \eqn{^{\circ}}{°}C}
  \item{round.it}{logical, round the values in the output?}
  \item{loga.protein}{numeric, decimal logarithms of reference activities of proteins}
}

\details{
  These functions accept \code{protein} (and optionally \code{organism}) in the same way as \code{\link{ip2aa}}, that is, as a protein name (optionally with the \code{organism} part separated), one or more row numbers in \code{thermo$protein} that can be identified using \code{\link{iprotein}}, or a data frame in the format of \code{thermo$protein}.

  \code{protein.formula} returns a stoichiometrix matrix representing the chemical formulas of the proteins that can be pased to e.g. \code{\link{mass}} or \code{\link{ZC}}. The amino acid compositions are multiplied by the output of \code{group.formulas} to generate the result. \code{group.formulas} returns the chemical formulas of each of the 20 common amino acid residues in proteins, as well as the terminal -H and -H (treated as the [H2O] group).

  \code{protein.length} returns the lengths (number of amino acids) of the proteins.

  \code{protein.info} tabulates some properties of proteins. A data frame is returned with a row for each protein, and columns named \samp{protein}, \samp{length}, \samp{formula}, \samp{G}, \samp{Z}, \samp{G.Z} and \samp{ZC}, indicating the names of the proteins, their lengths, chemical formulas, and values of the standard molal Gibbs energy of the neutral (nonionized) proteins, net charges and standard molal Gibbs energy of the ionized proteins, and average oxidation states of carbon. \samp{Z} and \samp{G.Z} are calculated using \code{\link{ionize.aa}} with values of pH taken from \code{thermo$basis}; \samp{Z} and \samp{G.Z} become NA if the \code{\link{basis}} species are not loaded or \samp{H+} is not in the basis definition. \samp{ZC} is calculated using \code{\link{ZC}}. The value of \code{T} indicates the temperature at which to calculate the Gibbs energies and net charge. The values of standard Gibbs energy are shown in cal/mol; these and other numeric values are rounded at a set number of digits if \code{round.it} is TRUE. The values (including chemical formula but not \samp{ZC}) are divided by the lengths of the proteins if \code{residue} is TRUE.

  The following two functions depend on an existing definition of the basis species:

\code{protein.basis} calculates the numbers of the basis species (i.e. opposite of the coefficients in the formation reactions) that can be combined to form the composition of each of the proteins.
The basis species must be present in \code{thermo$basis}, and if \samp{H+} is among the basis species, the ionization states of the proteins are included.
As with \code{protein.info}, the ionization state of the protein is calculated at the pH defined in \code{thermo$basis} and at the temperature specified by the \code{T} argument.
If \code{normalize} is TRUE, the coefficients on the basis species are divided by the lengths of the proteins. 

  \code{protein.equil} produces a series of messages showing step-by-step a calculation of the chemical activities of proteins in metastable equilibrium. For the first protein, it shows the standard Gibbs energies of the reaction to form the nonionized protein from the basis species and of the ionization reaction of the protein (if \samp{H+} is in the basis), then the standard Gibbs energy/RT of the reaction to form the (possibly ionized) protein per residue. The per-residue values of \samp{logQstar} and \samp{Astar/RT} are also shown for the first protein. Equilibrium calculations are then performed, only if more than one protein is specified. This calculation applies the Boltzmann distribution to the calculation of the equilibrium degrees of formation of the residue equivalents of the proteins, then converts them to activities of proteins taking account of \code{loga.protein} and protein length. If the \code{protein} argument is numeric (indicating rownumbers in \code{thermo$protein}), the values of \samp{Astar/RT} are compared with the output of \code{\link{affinity}}, and those of the equilibrium degrees of formation of the residues and the chemical activities of the proteins with the output of \code{\link{diagram}}. If the values in any of these tests are are not \code{\link{all.equal}} an error is produced indicating a bug. 

  \code{MP90.cp} takes \code{protein} (name of protein) and \code{T} (one or more temperatures in \eqn{^{\circ}}{°}C) and returns the additive heat capacity (J mol \eqn{^{-1}}{^-1}) of the unfolded protein using values of heat capacities of the residues taken from Makhatadze and Privalov, 1990. Those authors provided values of heat capacity at six points between 5 and 125 \eqn{^{\circ}}{°}C; this function interpolates (using \code{\link{splinefun}}) values at other temperatures.
}

\seealso{
  \code{\link{ionize.aa}} for an example that compares \code{MP90.cp} with heat capacities calculated in CHNOSZ at different temperatures and pHs. The functions for interacting with the database of amino acid compositions of proteins are documented at \code{\link{iprotein}}, and examples of relative stability calculations can be found on the \code{\link{protein}} help page.
}

\examples{\dontshow{data(thermo)}
## example for chicken lysozyme C
# index in thermo$protein
ip <- iprotein("LYSC_CHICK")
# amino acid composition
ip2aa(ip)
# length and chemical formula
protein.length(ip)
protein.formula(ip)
# formula, Gibbs energy, average oxidation state of carbon
protein.info(ip)
# as above, now with charge and Gibbs energy of ionized protein at pH 7
basis("CHNOS+")
protein.info(ip)
# group additivity for thermodynamic properties and HKF equation-of-state
# parameters of non-ionized protein
aa2eos(ip2aa(ip))
# calculation of standard thermodynamic properties
# (subcrt uses the species name, not ip)
subcrt("LYSC_CHICK")
# affinity calculation, protein identified by ip
affinity(iprotein=ip)
# affinity calculation, protein loaded as a species
species("LYSC_CHICK")
affinity()
# NB: subcrt() only shows the properties of the non-ionized
# protein, but affinity() uses the properties of the ionized
# protein if the basis species have H+

## these are all the same
protein.formula("P53_PIG")
protein.formula(iprotein("P53_PIG"))
protein.formula(ip2aa(iprotein("P53_PIG")))

## steps in calculation of chemical activities of two proteins
## in metastable equilibrium, after Dick and Shock, 2011
protein <- iprotein(c("CSG_METVO", "CSG_METJA"))
# clear out amino acid residues loaded by the example above
# ( in affinity(iprotein=ip) )
data(thermo)
# load supplemental database to use "old" [Met] sidechain group
add.obigt()
# set up the basis species to those used in DS11
basis("CHNOS+")
# note this yields logaH2 = -4.657486
swap.basis("O2", "H2")
# demonstrate the steps of the equilibrium calculation
protein.equil(protein, loga.protein=-3)
## we can also look at the affinities
# (Reaction 7, Dick and Shock, 2011)
# A/2.303RT for protein at unit activity (A-star for the protein)
a <- affinity(iprotein=protein[1], loga.protein=0)
Astar.protein <- a$values[[1]]
# divide affinity by protein length (A-star for the residue)
pl <- protein.length(protein[1])
Astar.residue <- a$values[[1]]/pl  # 0.1893, Eq. 11
# A/2.303RT per residue corresponding to protein activity of 10^-3
loga.residue <- log10(pl*10^-3)
Aref.residue <- Astar.residue - loga.residue  # 0.446, after Eq. 16
# A-star of the residue in natural log units (A/RT)
log(10) * Astar.residue  # 0.4359, after Eq. 23

\donttest{## using protein.formula: average oxidation state of 
## carbon of proteins from different organisms
# get amino acid compositions of microbial proteins 
# generated from the RefSeq database 
file <- system.file("extdata/refseq/protein_refseq.csv.xz", package="CHNOSZ")
ip <- add.protein(read.aa(file))
# only use those organisms with a certain
# number of sequenced bases
ip <- ip[as.numeric(thermo$protein$abbrv[ip]) > 100000]
pf <- protein.formula(thermo$protein[ip, ])
zc <- ZC(pf)
# the organism names we search for
# "" matches all organisms
terms <- c("Natr", "Halo", "Rhodo", "Acido", "Methylo",
  "Nitro", "Desulfo", "Chloro", "Geo", "Methano",
  "Thermo", "Pyro", "Sulfo", "Buchner", "")
tps <- thermo$protein$ref[ip]
plot(0, 0, xlim=c(1, 15), ylim=c(-0.3, -0.05), pch="",
  ylab="average oxidation state of carbon in proteins",
  xlab="", xaxt="n", mar=c(6, 3, 1, 1))
for(i in 1:length(terms)) {
  it <- grep(terms[i], tps)
  zct <- zc[it]
  points(jitter(rep(i, length(zct))), zct, pch=20)
}
terms[15] <- paste("all", length(ip))
axis(1, 1:15, terms, las=2)
title(main=paste("Average Oxidation State of Carbon:",
  "Total Protein per taxID in NCBI RefSeq", sep="\n"))
}
}

\references{
  Dick, J. M. and Shock, E. L. (2011) Calculation of the relative chemical stabilities of proteins as a function of temperature and redox chemistry in a hot spring. \emph{PLoS ONE} \bold{6}, e22782. \url{http://dx.doi.org/10.1371/journal.pone.0022782}

  Makhatadze, G. I. and Privalov, P. L. (1990) Heat capacity of proteins. 1. Partial molar heat capacity of individual amino acid residues in aqueous solution: Hydration effect \emph{J. Mol. Biol.} \bold{213}, 375--384. \url{http://dx.doi.org/10.1016/S0022-2836(05)80197-4}
}

\keyword{protein}
