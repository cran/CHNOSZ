## a case where the maximum affinity method doesn't
## reproduce an equilibrium predominance diagram
basis("CHNOS+")
# this adds data for some metabolites in the TCA cycle
# from Dalla-Betta and Schulte, 2010
add.obigt()
species(c("oxaloacetate-2", "malate-2", "fumarate-2",
  "a-ketoglutarate-2", "citrate-3"))
a <- affinity(O2=c(-80, -60), H2O=c(-5, 5))
diagram(a, dotted=1, fill="heat")
e <- equilibrate(a)
diagram(e, add=TRUE, names=NULL, col="purple")
e <- equilibrate(a, normalize=TRUE)
diagram(e, add=TRUE, names=NULL)
title(main=paste("maximum affinity method (fields)\n",
  "equilibrium calculations (lines)"))
data(thermo)
