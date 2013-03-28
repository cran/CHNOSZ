### R code from vignette source 'equilibrium.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: equilibrium.Rnw:23-24
###################################################
  if(exists(".orig.enc")) options(encoding = .orig.enc)


###################################################
### code chunk number 2: add_obigt
###################################################
library(CHNOSZ)
data(thermo)
add.obigt()


###################################################
### code chunk number 3: ProteinFormation
###################################################
basis("CHNOS+")
species("CSG",c("METVO", "METJA"))


###################################################
### code chunk number 4: ProteinInfo
###################################################
protein.info(species()$name)


###################################################
### code chunk number 5: ProteinAffinity
###################################################
a <- affinity()
a$values


###################################################
### code chunk number 6: ProteinActivities
###################################################
e <- equilibrate(a)
e$loga.equil


###################################################
### code chunk number 7: equilibrium.Rnw:334-335
###################################################
protein.basis(species()$name, normalize=TRUE)


###################################################
### code chunk number 8: protein_equil
###################################################
# get an error if we don't data(thermo), only in the re-building vignettes of R CMD check
data(thermo)
protein <- iprotein(c("CSG_METVO", "CSG_METJA"))
basis("CHNOS+")
swap.basis("O2", "H2")
protein.equil(protein, loga.protein=-3)


###################################################
### code chunk number 9: ProteinSpeciation
###################################################
organisms <- c("METSC", "METJA", "METFE", "HALJP",
  "METVO", "METBU", "ACEKI", "GEOSE", "BACLI", "AERSA")
proteins <- c(rep("CSG", 6), rep("SLAP", 4))
basis("CHNOS+")
species(proteins, organisms)
a <- affinity(O2=c(-100, -65))
par(mfrow=c(2, 1))
e <- equilibrate(a)
diagram(e, ylim=c(-5, -1), legend.x=NA)
title(main="Equilibrium activities of proteins, whole formulas")
e <- equilibrate(a, normalize=TRUE)
diagram(e, ylim=c(-5, -1), legend.x=NA)
title(main="Equilibrium activities of proteins, normalized formulas")


###################################################
### code chunk number 10: SulfurSpeciation
###################################################
basis("CHNOS+")
basis("pH",5)
species(c("H2S", "S2-2", "S3-2", "S2O3-2", "S2O4-2",
  "S3O6-2", "S5O6-2", "S2O6-2", "HSO3-", "SO2", "HSO4-"))
a <- affinity(O2=c(-50, -15), T=325, P=350)
par(mfrow=c(2, 1))
e <- equilibrate(a, loga.balance=-2)
diagram(e, ylim=c(-30, 0), legend.x="topleft", cex.names=0.8)
title(main="Aqueous sulfur speciation, whole formulas")
e <- equilibrate(a, loga.balance=-2, normalize=TRUE)
diagram(e, ylim=c(-30, 0), legend.x="topleft", cex.names=0.8)
title(main="Aqueous sulfur speciation, normalized formulas")


###################################################
### code chunk number 11: Plasma
###################################################
data(thermo)  # cleanup from previous plot
basis(c("CO2", "NH3", "H2S", "H2O", "O2"), c(-3, -3, -10))
f <- system.file("extdata/abundance/AA03.csv", package="CHNOSZ")
pdat <- read.csv(f, as.is=TRUE)
iil <- grep("^IL", pdat$name)
species(pdat$name[iil], "HUMAN")
a <- affinity(O2=c(-82, -78), H2O=c(-12, -2))
par(mfrow=c(1, 2))
dA <- diagram(a, normalize=TRUE, main="maximum affinity")
e <- equilibrate(a, normalize=TRUE)
dE <- diagram(e, main="equilibrium activities")
stopifnot(identical(dA$predominant, dE$predominant))


###################################################
### code chunk number 12: Amino
###################################################
basis("CHNOS+")
species(aminoacids(""))
a <- affinity(O2=c(-71, -66), H2O=c(-8, 4))
par(mfrow=c(1, 2))
dA <- diagram(a, main="maximum affinity")
e <- equilibrate(a)
dE <- diagram(e, main="equilibrium activities")


###################################################
### code chunk number 13: SessionInfo
###################################################
sessionInfo()


