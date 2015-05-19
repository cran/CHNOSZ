\name{protein}
\alias{protein}
\title{Examples of Calculations for Proteins}
\description{
  This page contains some examples of using the functions in CHNOSZ to calculate thermodynamic properties and relative stabilities of proteins.
}

\examples{\dontshow{data(thermo)}
## thermodynamic properties and activity diagrams
prot4 <- c("LYSC_CHICK", "BPT1_BOVIN", "CYC_BOVIN", "MYG_PHYCA")
# aqueous protein properties (nonionized)
subcrt(prot4, T=seq(0, 200, 10))$out
# T-logfO2 activity diagram
basis("CHNOS")
species(prot4)
a <- affinity(T=c(0, 200), O2=c(-80, -40))
diagram(a)

## one way to calculate the standard Gibbs energy of a 
## reaction to form an ionized protein at 100 degrees and pH 8
basis("CHNOS+")   # do this to auto-balance the formation reaction
Gr.nonionized <- subcrt("LYSC_CHICK", 1, T=100)$out$G
basis("pH", 8)
pinfo <- protein.info("LYSC_CHICK", round.it=FALSE, T=100)
Gr.ionization <- pinfo$G.Z - pinfo$G
# standard Gibbs energy of the reaction
# in cal/mol ionized protein:
Gr.ionized <- Gr.nonionized + Gr.ionization
# in cal/mol ionized residue:
Gr.ionized_residue <- Gr.ionized/protein.length("LYSC_CHICK")

## Standard molal entropy of a protein reaction
basis("CHNOS")
# here we provide the reaction coefficients of the 
# proteins (per protein backbone); 'subcrt' function calculates 
# the coefficients of the basis species in the reaction
s <- subcrt(c("CSG_METTL","CSG_METJA"), c(-1/530,1/530),
  T=seq(0, 350, length.out=50))
thermo.plot.new(xlim=range(s$out$T), ylim=range(s$out$S),
  xlab=axis.label("T"), ylab=axis.label("DS0r"))
lines(s$out$T, s$out$S)
# do it at high pressure as well
s <- subcrt(c("CSG_METTL","CSG_METJA"), c(-1/530,1/530),
  T=seq(0,350,length.out=50), P=3000)
lines(s$out$T, s$out$S, lty=2)
# label the plot
title(main=paste("Standard molal entropy\n",
  "P = Psat (solid), P = 3000 bar (dashed)"))
s$reaction$coeff <- round(s$reaction$coeff, 3)
dsr <- describe.reaction(s$reaction, iname=c(1,2))
text(170, -3, dsr, cex=0.8)


### Equilibrium activity diagrams

## surface-layer proteins from Methanococcus and others
## as a function of oxygen fugacity, after Dick, 2008, Fig. 5b
# use old properties of [Met] to reproduce this example
data(thermo)
add.obigt()
# make our protein list
organisms <- c("METSC", "METJA", "METFE", "HALJP", "METVO",
  "METBU", "ACEKI", "GEOSE", "BACLI", "AERSA")
proteins <- c(rep("CSG", 6), rep("SLAP", 4))
proteins <- paste(proteins, organisms, sep="_")
# load the basis species and proteins
basis("CHNOS+")
species(proteins)
# calculate affinities; we go to lower logfO2 than Dick, 2008
# and find an interesting convergence of stabilities there
a <- affinity(O2=c(-100, -65))
# try normalize=FALSE to make Fig. 5a in the paper
e <- equilibrate(a, normalize=TRUE)
d <- diagram(e, ylim=c(-5, -1), legend.x=NA, names=organisms)
# add water stability line
abline(v=-83.1, lty=2)
title(main="Surface-layer proteins, after Dick, 2008")
# checking the geometry of the diagram
# most preominant along the x-axis
stopifnot(organisms[unique(which.pmax(e$loga.equil))] ==
  c("METFE", "METJA", "METVO", "HALJP"))
# stability order close to logfO2=-83.1
stopifnot(order(as.data.frame(e$loga.equil)[62,],
  decreasing=TRUE)==c(2, 6, 7, 5, 3, 1, 9, 8, 10, 4))
# reset thermodynamic database
data(thermo)

## relative stabilities of bovine proteins
## as a function of temperature along a glutathione redox buffer
mod.buffer("GSH-GSSG", c("GSH","GSSG"), logact=c(-3, -7))   
basis(c("CO2", "H2O", "NH4+", "SO4-2", "H2", "H+"),
  c(-1, 0, -4, -4, 999, -7)) 
basis("H2", "GSH-GSSG")
basis("CO2", "gas")
prot <- c("CYC", "RNAS1", "BPT1", "ALBU", "INS", "PRIO")
species(prot, "BOVIN")
a <- affinity(T=c(0, 200))
# set line colors according to oxidation state of carbon
ZC <- ZC(species()$ispecies)
col <- rep("red", length(prot))
col[ZC > 0] <- "blue"
e <- equilibrate(a, normalize=TRUE)
d <- diagram(e, col=col, legend.x=NA)
title(main="Bovine proteins, GSH/GSSG redox buffer")
}

\seealso{
  For accessing, updating, and downloading amino acid compositions of proteins, see \code{\link{iprotein}}. 
  For getting chemical formulas and stoichiometric coefficients in reactions of proteins, see \code{\link{protein.info}}.
For more examples of metastable equilibrium calculations for proteins, see \code{\link{read.expr}}, \code{\link{more.aa}}, \code{\link{ionize.aa}}, and \code{\link{apc}} for reaction-path calculations.
}

\concept{Protein thermodynamic modeling}
