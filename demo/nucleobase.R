## Nucleobase - Amino Acid Interaction Eh-H2O
# for this example we try a unique basis definition
basis(c("CO2","H2O","glutamine","e-","H+"),c(-3,0,-3,0,-7))
species(c("uracil","cytosine","adenine","guanine",
  "phenylalanine","proline","lysine","glycine"),"aq")
# this loaded four nucleobases and four related amino acids
# (coded for by the homocodon triplets)
# check out the predominance diagrams
a.1 <- affinity(H2O=c(-5,0),Eh=c(-0.5,0))
diagram(a.1,fill=NULL)
# overlay a different temperature
a.2 <- affinity(H2O=c(-5,0),Eh=c(-0.5,0),T=100)
diagram(a.2,col="red",add=TRUE,names=NULL)
# start make a title for the plot
tb <- thermo$basis   # includes activities of basis species
# exclude those that are on the axes
tb <- tb[!((rownames(tb) %in% c("e-","H2O"))),]
title(main="Nucleobases and amino acids; P=Psat")
dp <- describe.property(c("T", "T"), c(25, 100))
db <- describe.basis(tb)
legend("bottomleft", lty=c(1, 1, NA, NA, NA), col=c("black","red"), legend=c(dp, db))
