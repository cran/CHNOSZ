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

\donttest{## Standard molal entropy of a protein reaction
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

## relative stabilities of plasma proteins,
## using chemical activities of H2 and O2
# find that insulin is very stable in oxidizing conditions
basis(c("CO2", "NH3", "H2S", "H2", "O2"), "aq", c(-3, -3, -10))
f <- system.file("extdata/abundance/AA03.csv", package="CHNOSZ")
pdat <- read.csv(f, as.is=TRUE)
species(pdat$name, "HUMAN")
basis("H2", -3)
a <- affinity(O2=c(-90, -75))
e <- equilibrate(a, normalize=TRUE)
diagram(e, col="green", legend.x=NA, ylim=c(-10, 1))
basis("H2", -20)
a <- affinity(O2=c(-90, -75))
e <- equilibrate(a, normalize=TRUE)
diagram(e, add=TRUE, col="blue")
legend("bottomleft", col=c(NA, "green", "blue"), lty=1,
  legend=c(axis.label("H2"), -3, -20))
title(main="Equilibrium activities of human plasma proteins")


### Calculations showing effects of ionization of proteins

## Eh-pH diagrams for intra/extracellular proteins
organism <- c("PYRFU", "ECOLI", "YEAST")
intracellular <- c("AMPM", "AMPM", "AMPM1")
extracellular <- c("O08452", "AMY1", "PST1")
basis("CHNOSe")  # for Eh we need electrons
mycol <- c("red", "green", "blue")
for(i in 1:3) {
  species(delete=TRUE)
  species(c(intracellular[i], extracellular[i]), organism[i])
  a <- affinity(pH=c(0, 14), Eh=c(-1, 0))
  diagram(a, add=i>1, fill=NULL, names=species()$name,
    col=mycol[i], col.names=mycol[i], normalize=TRUE)
}
title(main="Intracellular (AMPM) and extracellular proteins")
text(9, -0.05, describe.basis(ibasis=1:4, oneline=TRUE))


## Buffer + ionization: Metastablilities of
## thiol peroxidases from model bactera
## (ECOLI, BACSU mesophile; AQUAE thermophile,
## THIDA acidophile, BACHD alkaliphile)
basis("CHNOS+")
organisms <- c("ECOLI", "AQUAE", "BACSU", "BACHD", "THIDA")
species("TPX", organisms)
# create a buffer with our proteins in it
mod.buffer("TPX", paste("TPX", organisms, sep="_"))
# set up the buffered activities
basis(c("CO2", "H2O", "NH3", "O2"), "TPX")
a <- affinity(return.buffer=TRUE, T=50)
basis(c("CO2", "H2O", "NH3", "O2"), as.numeric(a[1:4]))
a <- affinity(pH=c(0, 14, 200), T=c(25, 70, 200))
diagram(a, fill=NULL)
title(main="Thiol peroxidases from bacteria")
text(0.5, 66, describe.basis(thermo$basis[-6,], oneline=TRUE), adj=0)


## Buffer + ionization: relative stabilities
## of E. coli sigma factors on a T-pH diagram
# (sigma factors 24, 32, 38, 54, 70, i.e.
# RpoE, RpoH, RpoS, RpoN, RpoD)
proteins <- c("RPOE", "RP32", "RPOS", "RP54", "RPOD")
basis("CHNOS+")
basis("pH", 7.4)
# define and set the buffer
mod.buffer("sigma", paste(proteins, "ECOLI", sep="_"))
basis(c("CO2", "NH3", "H2S", "O2"), "sigma")
logact <- affinity(return.buffer=TRUE, T=25)
# Set the activities of the basis species to constants 
# corresponding to the buffer, and diagram the relative
# stabilities as a function of T and pH
basis(c("CO2", "NH3", "H2S", "O2"), as.numeric(logact))
species(paste(proteins, "ECOLI", sep="_"))
a <- affinity(pH=c(5, 10), T=c(10, 40))
diagram(a, normalize=FALSE, fill="heat")
title(main="Relative stabilities of sigma factors in E. coli")
ptext <- c(describe.property("T", 25), 
  describe.basis(ibasis=c(2, 6), oneline=TRUE))
btext <- describe.basis(ibasis=c(1, 3, 4, 5), oneline=TRUE)
legend("bottomleft", legend=c("preset (input values):",
  ptext, "buffered (results):", btext), bty="n")
}
}

\seealso{
  For accessing, updating, and downloading amino acid compositions of proteins, see \code{\link{iprotein}}. 
  For getting chemical formulas and stoichiometric coefficients in reactions of proteins, see \code{\link{protein.info}}.
For more examples of metastable equilibrium calculations for proteins, see \code{\link{read.expr}}, \code{\link{more.aa}}, \code{\link{ionize.aa}}, and \code{\link{apc}} for reaction-path calculations.
}

\keyword{protein}
