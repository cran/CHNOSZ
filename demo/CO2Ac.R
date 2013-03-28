# one can solve for the logarithm of activity of a 
# basis species using the 'what' argument of diagram
basis("CHNOS")
basis("CO2", 999)
species("acetic acid", -3)
a <- affinity(O2=c(-85, -70, 4), T=c(25, 100, 4))
# hacking to write a title with formulas and subscripts
lCO2 <- axis.label("CO2")
main <- substitute(a~~b~~c,list(a=lCO2, b="buffered by",
  c="acetic acid"))
d <- diagram(a, what="CO2", main=main)
species(1, -10)
a <- affinity(O2=c(-85, -70, 4), T=c(25, 100, 4))
d <- diagram(a, what="CO2", add=TRUE, lty=2)
# add a legend
lAC <- expr.species("CH3COOH", log="aq")
ltext <- c(as.expression(lAC), -3, -10)
lty <- c(NA, 1, 2)
legend("topright", legend=ltext, lty=lty, bg="white")
# do return.buffer and diagram(what) give the same results?
and <- as.numeric(d$plotvals[[1]])
basis("CO2", "AC")
mod.buffer("AC", logact=-10)
a.buffer <- affinity(O2=c(-85, -70, 4), T=c(25, 100, 4), return.buffer=TRUE)
ana <- as.numeric(unlist(a.buffer[[1]]))
stopifnot(all.equal(ana, and))
