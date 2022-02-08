## ----setup, include=FALSE-----------------------------------------------------
## use pngquant to optimize PNG images
library(knitr)
knit_hooks$set(pngquant = hook_pngquant)
pngquant <- "--speed=1 --quality=0-25"
if (!nzchar(Sys.which("pngquant"))) pngquant <- NULL

## ----CHNOSZ_reset, include=FALSE----------------------------------------------
library(CHNOSZ)
reset()

## ----AAsetup, results = "hide", message = FALSE-------------------------------
library(CHNOSZ)
reset()
basis("CHNOS")
species(aminoacids(""))

## ----AAfunctions--------------------------------------------------------------
aaA <- function() {
  a <- affinity(O2 = c(-90, -70), H2O = c(-20, 10))
  diagram(a, balance = 1, names = aa)
}

aaB <- function() {
  a <- affinity(O2 = c(-90, -70), H2O = c(-20, 10))
  e <- equilibrate(a, balance = 1)
  diagram(e, names = aa)
}

aaC <- function() {
  a <- affinity(O2 = c(-71, -66), H2O = c(-8, 4))
  diagram(a, balance = "CO2", names = aa)
}

aaD <- function() {
  a <- affinity(O2 = c(-71, -66), H2O = c(-8, 4))
  e <- equilibrate(a, balance = "CO2")
  diagram(e, names = aa)
}

aaE <- function() {
  basis("O2", -66)
  a <- affinity(H2O = c(-8, 4))
  e <- equilibrate(a, balance = "CO2")
  diagram(e, ylim = c(-5, -1), names = aa)
}

aaF <- function() {
  species(1:20, -7)
  a <- affinity(H2O = c(-8, 4))
  e <- equilibrate(a, balance = "CO2")
  diagram(e, ylim = c(-8, -4), names = aa)
}

## ----AAabbrv------------------------------------------------------------------
aa <- aminoacids()
aa

## ----AAplot, eval = FALSE-----------------------------------------------------
#  showtime <- function(st) {
#    # plot time in lower-right of figure region
#    f <- usrfig()
#    par(xpd=TRUE)
#    if(st[3] > 2) col <- "red" else col <- "black"
#    text(f$x[2], f$y[1], paste(round(st[3], 1), "s\n"), adj=1, col=col)
#    par(xpd=FALSE)
#  }
#  
#  layout(t(matrix(c(1:7, 11, 8:10, 12), nrow=4)), widths=c(1, 4, 4, 4), heights=c(0.7, 4, 4))
#  
#  ## row 0 (column titles)
#  opar <- par(mar=c(0, 0, 0, 0))
#  plot.new()
#  plot.new()
#  text(0.58, 0.5, "maximum affinity", cex=1.4)
#  plot.new()
#  text(0.58, 0.5, "equilibration", cex=1.4)
#  plot.new()
#  text(0.58, 0.5, "equilibration", cex=1.4)
#  par(opar)
#  
#  ## row 1 (balance = 1)
#  opar <- par(mar=c(0, 0, 0, 0))
#  plot.new()
#  text(0.5, 0.5, "balance = 1", srt=90, cex=1.4)
#  par(opar)
#  # figure A
#  st <- system.time(dA <- aaA())
#  showtime(st)
#  title(main="loga(species) = -3", cex.main=1)
#  label.figure("A", yfrac=0.92, xfrac=0.1, font = 2)
#  # figure B
#  st <- system.time(dB <- aaB())
#  showtime(st)
#  title(main=paste("loga(total species) =", round(dB$loga.balance[1], 2)), cex.main=1)
#  label.figure("C", yfrac=0.92, xfrac=0.1, font = 2)
#  
#  ## row 2 (balance = nCO2)
#  opar <- par(mar=c(0, 0, 0, 0))
#  plot.new()
#  text(0.5, 0.5, 'balance = "CO2"', srt=90, cex=1.4)
#  par(opar)
#  # figure C
#  st <- system.time(dC <- aaC())
#  showtime(st)
#  title(main="loga(species) = -3", cex.main=1)
#  label.figure("B", yfrac=0.92, xfrac=0.1, font = 2)
#  # figure D
#  st <- system.time(dD <- aaD())
#  showtime(st)
#  title(main=paste("loga(total CO2) =", round(dD$loga.balance[1], 2)), cex.main=1)
#  label.figure("D", yfrac=0.92, xfrac=0.1, font = 2)
#  
#  ## right (speciation at different total activity of CO2)
#  par(xpd=NA)
#  lines(c(-66, -64.5), c(4, 9), lty=2)
#  lines(c(-66, -64.5), c(-8, -8.5), lty=2)
#  par(xpd=FALSE)
#  # figure E
#  st <- system.time(dE <- aaE())
#  showtime(st)
#  title(main=paste("loga(total CO2) =", round(dE$loga.balance[1], 2)), cex.main=1)
#  label.figure("E", yfrac=0.92, xfrac=0.1, font = 2)
#  # figure F
#  st <- system.time(dF <- aaF())
#  showtime(st)
#  title(main=paste("loga(total CO2) =", round(dF$loga.balance[1], 2)), cex.main=1)
#  label.figure("F", yfrac=0.92, xfrac=0.1, font = 2)

## ----AAplot, echo = FALSE, results = "hide", message = FALSE, fig.width = 13/2, fig.height = 8.7/2, out.width = "100%", pngquant = pngquant----
showtime <- function(st) {
  # plot time in lower-right of figure region
  f <- usrfig()
  par(xpd=TRUE)
  if(st[3] > 2) col <- "red" else col <- "black"
  text(f$x[2], f$y[1], paste(round(st[3], 1), "s\n"), adj=1, col=col)
  par(xpd=FALSE)
}

layout(t(matrix(c(1:7, 11, 8:10, 12), nrow=4)), widths=c(1, 4, 4, 4), heights=c(0.7, 4, 4))

## row 0 (column titles)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
plot.new()
text(0.58, 0.5, "maximum affinity", cex=1.4)
plot.new()
text(0.58, 0.5, "equilibration", cex=1.4)
plot.new()
text(0.58, 0.5, "equilibration", cex=1.4)
par(opar)

## row 1 (balance = 1)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "balance = 1", srt=90, cex=1.4)
par(opar)
# figure A
st <- system.time(dA <- aaA())
showtime(st)
title(main="loga(species) = -3", cex.main=1)
label.figure("A", yfrac=0.92, xfrac=0.1, font = 2)
# figure B
st <- system.time(dB <- aaB())
showtime(st)
title(main=paste("loga(total species) =", round(dB$loga.balance[1], 2)), cex.main=1)
label.figure("C", yfrac=0.92, xfrac=0.1, font = 2)

## row 2 (balance = nCO2)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, 'balance = "CO2"', srt=90, cex=1.4)
par(opar)
# figure C
st <- system.time(dC <- aaC())
showtime(st)
title(main="loga(species) = -3", cex.main=1)
label.figure("B", yfrac=0.92, xfrac=0.1, font = 2)
# figure D
st <- system.time(dD <- aaD())
showtime(st)
title(main=paste("loga(total CO2) =", round(dD$loga.balance[1], 2)), cex.main=1)
label.figure("D", yfrac=0.92, xfrac=0.1, font = 2)

## right (speciation at different total activity of CO2)
par(xpd=NA)
lines(c(-66, -64.5), c(4, 9), lty=2)
lines(c(-66, -64.5), c(-8, -8.5), lty=2)
par(xpd=FALSE)
# figure E
st <- system.time(dE <- aaE())
showtime(st)
title(main=paste("loga(total CO2) =", round(dE$loga.balance[1], 2)), cex.main=1)
label.figure("E", yfrac=0.92, xfrac=0.1, font = 2)
# figure F
st <- system.time(dF <- aaF())
showtime(st)
title(main=paste("loga(total CO2) =", round(dF$loga.balance[1], 2)), cex.main=1)
label.figure("F", yfrac=0.92, xfrac=0.1, font = 2)

## ----PRsetup, results = "hide", message = FALSE-------------------------------
basis("CHNOS+")
organisms <- c("METJA", "HALJP", "METVO", "ACEKI", "GEOSE", "BACLI")
proteins <- c(rep("CSG", 3), rep("SLAP", 3))
species(proteins, organisms)

## ----PRlength-----------------------------------------------------------------
protein.length(species()$name)

## ----PRfunctions--------------------------------------------------------------
prA <- function() {
  a <- affinity(O2 = c(-90, -70), H2O = c(-20, 0))
  diagram(a, balance = "length", names = organisms)
}

prB <- function() {
  a <- affinity(O2 = c(-90, -70))
  e <- equilibrate(a, balance = "length")
  ylab <- quote(log~italic(a)~protein)
  diagram(e, names = organisms, ylim = c(-5, -1), ylab = ylab)
}

prC <- function() {
  a <- affinity(O2 = c(-90, -70), H2O = c(-20, 0))
  diagram(a, normalize = TRUE, names = organisms)
}

prD <- function() {
  a <- affinity(O2 = c(-90, -70))
  e <- equilibrate(a, normalize = TRUE)
  ylab <- quote(log~italic(a)~protein)
  diagram(e, names = organisms, ylim = c(-5, -1), ylab = ylab)
}

prE <- function() {
  a <- affinity(O2 = c(-90, -70), H2O = c(-20, 0))
  diagram(a, as.residue = TRUE, names = organisms)
}

prF <- function() {
  a <- affinity(O2 = c(-90, -70))
  e <- equilibrate(a, as.residue = TRUE, loga.balance = 0)
  ylab <- quote(log~italic(a)~residue)
  diagram(e, names = organisms, ylim = c(-3, 1), ylab = ylab)
}

## ----PRplot, eval = FALSE-----------------------------------------------------
#  layout(t(matrix(1:12, nrow=4)), widths=c(1, 4, 4, 4), heights=c(0.7, 4, 4))
#  
#  ## row 0 (column titles)
#  opar <- par(mar=c(0, 0, 0, 0))
#  plot.new()
#  plot.new()
#  text(0.58, 0.5, 'balance = "length"', cex=1.4)
#  plot.new()
#  text(0.58, 0.5, "normalize = TRUE\n(balance = 1)", cex=1.4)
#  plot.new()
#  text(0.58, 0.5, "as.residue = TRUE\n(balance = 1)", cex=1.4)
#  par(opar)
#  
#  ## row 1 (maximum affinity 2D)
#  opar <- par(mar=c(0, 0, 0, 0))
#  plot.new()
#  text(0.5, 0.5, "maximum affinity", srt=90, cex=1.4)
#  par(opar)
#  # figure A (balance = "length")
#  st <- system.time(dA <- prA())
#  showtime(st)
#  label.figure("A", yfrac=0.9, xfrac=0.1, font = 2)
#  # figure C (normalize = TRUE)
#  st <- system.time(dC <- prC())
#  showtime(st)
#  label.figure("C", yfrac=0.9, xfrac=0.1, font = 2)
#  # figure E (as.residue = TRUE)
#  st <- system.time(dE <- prE())
#  showtime(st)
#  label.figure("E", yfrac=0.9, xfrac=0.1, font = 2)
#  
#  ## row 2 (equilibrate 1D)
#  opar <- par(mar=c(0, 0, 0, 0))
#  plot.new()
#  text(0.5, 0.5, "equilibration", srt=90, cex=1.4)
#  par(opar)
#  # figure B (balance = "length")
#  st <- system.time(prB())
#  showtime(st)
#  label.figure("B", yfrac=0.9, xfrac=0.1, font = 2)
#  # figure D (normalize = TRUE)
#  st <- system.time(prD())
#  showtime(st)
#  label.figure("D", yfrac=0.9, xfrac=0.1, font = 2)
#  # figure F (as.residue = TRUE)
#  st <- system.time(prF())
#  showtime(st)
#  label.figure("F", yfrac=0.9, xfrac=0.1, font = 2)

## ----PRplot, echo = FALSE, results = "hide", message = FALSE, fig.width = 13/2, fig.height = 8.7/2, out.width = "100%", pngquant = pngquant----
layout(t(matrix(1:12, nrow=4)), widths=c(1, 4, 4, 4), heights=c(0.7, 4, 4))

## row 0 (column titles)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
plot.new()
text(0.58, 0.5, 'balance = "length"', cex=1.4)
plot.new()
text(0.58, 0.5, "normalize = TRUE\n(balance = 1)", cex=1.4)
plot.new()
text(0.58, 0.5, "as.residue = TRUE\n(balance = 1)", cex=1.4)
par(opar)

## row 1 (maximum affinity 2D)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "maximum affinity", srt=90, cex=1.4)
par(opar)
# figure A (balance = "length")
st <- system.time(dA <- prA())
showtime(st)
label.figure("A", yfrac=0.9, xfrac=0.1, font = 2)
# figure C (normalize = TRUE)
st <- system.time(dC <- prC())
showtime(st)
label.figure("C", yfrac=0.9, xfrac=0.1, font = 2)
# figure E (as.residue = TRUE)
st <- system.time(dE <- prE())
showtime(st)
label.figure("E", yfrac=0.9, xfrac=0.1, font = 2)

## row 2 (equilibrate 1D)
opar <- par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "equilibration", srt=90, cex=1.4)
par(opar)
# figure B (balance = "length")
st <- system.time(prB())
showtime(st)
label.figure("B", yfrac=0.9, xfrac=0.1, font = 2)
# figure D (normalize = TRUE)
st <- system.time(prD())
showtime(st)
label.figure("D", yfrac=0.9, xfrac=0.1, font = 2)
# figure F (as.residue = TRUE)
st <- system.time(prF())
showtime(st)
label.figure("F", yfrac=0.9, xfrac=0.1, font = 2)

## ----ProteinSpeciation, results = "hide", message = FALSE, fig.width = 8, fig.height = 5.5, out.width = "100%", pngquant = pngquant----
organisms <- c("METSC", "METJA", "METFE",  "METVO", "METBU",
               "HALJP", "ACEKI", "GEOSE", "BACLI", "AERSA")
proteins <- c(rep("CSG", 6), rep("SLAP", 4))
# use red for Methano* genera
col <- c(rep(2, 5), rep(1, 5))
basis("CHNOS+")
species(proteins, organisms)
a <- affinity(O2 = c(-100, -65))
layout(matrix(1:2), heights = c(1, 2))
e <- equilibrate(a)
diagram(e, ylim = c(-2.8, -1.6), names = organisms, col = col)
water.lines(e, col = 4)
title(main="Equilibrium activities of proteins, normalize = FALSE")
e <- equilibrate(a, normalize = TRUE)
diagram(e, ylim = c(-4, -2), names = organisms, col = col)
water.lines(e, col = 4)
title(main = "Equilibrium activities of proteins, normalize = TRUE")

