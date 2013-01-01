opar <- par(mfrow=c(2, 2))
# find chemical activities where metastable equilibrium activities of
# selected proteins in P. ubique have high correlation
# with a lognormal distribution (i.e., maximize r of q-q plot)
f <- system.file("extdata/fasta/HTCC1062.faa.xz", package="CHNOSZ")
# search for ribosomal proteins
g <- grep.file(f, "ribosomal")
# read their amino acid compositions from the file
aa <- read.fasta(f, g)
# add these proteins to thermo$protein
ip <- add.protein(aa)
# load a predefined set of uncharged basis species
# (speeds things up as we won't model protein ionization)
basis("CHNOS")
# plot 1: calculated logarithms of chemical activity
# as a function of logfO2 ... a bundle of curves near logfO2 = -77
a <- affinity(O2=c(-90, -60), iprotein=ip)
e <- equilibrate(a, loga.balance=0, normalize=TRUE)
d <- diagram(e, names=NULL)
title(as.expression("Ribosomal proteins in"~italic("Pelagibacter ubique")))
db <- describe.basis(ibasis=c(2, 1, 3))
legend("bottomright", legend=db, bg="white")
# plot 2: calculate q-q correlation coefficient
# the lognormal distribution is favored near logfO2 = -77
r <- revisit(d, "qqr")
title(main="correlation with a normal distribution")
text(r$xopt, r$optimum, paste(" qqr", round(r$optimum, 3), sep="="), adj=0)
# plot 3: findit... maximize qqr as a function of activities of O2-H2O-NH3-CO2
f1 <- findit(list(O2=c(-106, -75), H2O=c(-40, -20), CO2=c(-20, 10), NH3=c(-15, 0)),
  "qqr", iprotein=ip, niter=12, normalize=TRUE)
title(main="searching 4-D chemical activity space")
# plot 4: q-q plot at the final loga O2, H2O, CO2, NH3
# higher correlation coefficient than plot 3
a <- affinity(iprotein=ip)
e <- equilibrate(a, loga.balance=0, normalize=TRUE)
qqr5 <- revisit(e, "qqr")$H
db <- describe.basis(ibasis=c(5, 2, 1, 3))
legend("bottomright", legend=db)
# plot 5: trajectory of O2, H2O, CO2, NH3, and the
# q-q correlation coefficient in the search
#plot(f1,mar=c(2,5,1,1),mgp=c(4,1,0))
par(opar)
