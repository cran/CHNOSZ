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
e <- equilibrate(a, normalize=TRUE)
diagram(e, fill=NULL)
title(main="Thiol peroxidases from bacteria")
text(0.5, 66, describe.basis(thermo$basis[-6,], oneline=TRUE), adj=0)
