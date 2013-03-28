### non-ideality calculations -- activity coefficients of 
### aqueous species as a function of charge, temperature,
### and ionic strength -- after Alberty, 2003 
## p. 16 Table 1.3  apparent pKa of acetic acid with
## changing ionic strength
subcrt(c("acetic acid","acetate","H+"),c(-1,1,1),
  IS=c(0,0.1,0.25),T=25,property="logK")
# note that these *apparent* values of G and logK approach
# their *standard* counterparts as IS goes to zero.
## p. 95: basis and elemental stoichiometries of species 
## (a digression here from the nonideality calculations) 
# note coefficient of O2 and NH3 will be zero for these species
basis(c("ATP-4","H+","H2O","HPO4-2","O2","NH3"))
# cf Eq. 5.1-33: (basis composition) 
species(c("ATP-4","H+","H2O","HPO4-2","ADP-3","HATP-3","HADP-2","H2PO4-"))
lb <- nrow(basis())
# cf Eq. 5.1-32: (elemental composition)
as.matrix(species()[,1:lb]) %*% as.matrix(basis()[,1:lb]) 
## p. 273-275: activity coefficient (gamma)
## as a function of ionic strength and temperature
## (as of 20080304, these do look quantitatively different 
## from the plots in Alberty's book.)
iplotfun <- function(T,col,add=TRUE) {
  IS <- seq(0,0.25,0.0025)
  s <- subcrt(c("H2PO4-","HADP-2","HATP-3","ATP-4"),IS=IS,grid="IS",T=T)
  if(!add) thermo.plot.new(xlim=range(IS),ylim=c(0,1),
    xlab=axis.label("IS"),ylab="gamma")
  for(i in 1:4) lines(IS,10^s$out[[i]]$loggam,col=col)
}
iplotfun(0,"blue",add=FALSE)
text(0.1, 0.62, "Z = -1")
iplotfun(25,"black")
text(0.075, 0.18, "Z = -2")
iplotfun(40,"red")
text(0.05, 0.06, "Z = -3")
title(main=paste("activity coefficient (gamma) of -1,-2,-3,-4",
  "charged species at 0, 25, 40 deg C, after Alberty, 2003",
  sep="\n"),cex.main=0.95)
legend("topright", lty=c(NA, 1, 1, 1), col=c(NA, "blue", "black", "red"),
  legend=c(as.expression(axis.label("T")), 0, 25, 40))
