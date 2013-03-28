### R code from vignette source 'hotspring.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: hotspring.Rnw:20-21
###################################################
  if(exists(".orig.enc")) options(encoding = .orig.enc)


###################################################
### code chunk number 2: options
###################################################
options(width=85,continue=" ")


###################################################
### code chunk number 3: libraryCHNOSZ
###################################################
library(CHNOSZ)
data(thermo)


###################################################
### code chunk number 4: add_obigt
###################################################
add.obigt()


###################################################
### code chunk number 5: classes
###################################################
tp <- thermo$protein
classes <- tp$protein[tp$organism=="bisonN"]
classes


###################################################
### code chunk number 6: sites
###################################################
sites <- c("N", "S", "R", "Q", "P")
sitenames <- paste("bison", sites, sep="")


###################################################
### code chunk number 7: ZCplot
###################################################
plot(0, 0, xlim=c(-1, 6), ylim=c(-0.33, -0.11), xlab="site", 
  ylab=expression(bar(italic(Z))[C]), xaxt="n")
axis(1, at=1:5)
col <- c("green", rep("black", 20))
lwd <- c(3, rep(1, 20))
clab <- c("hydrolase", "overall", "protease", 
  "oxidoreductase", "transport", "membrane", "permease")
ip <- iprotein(rep(classes, each=5), rep(sitenames, 20))
p <- thermo$protein[ip, ]
pf <- protein.formula(p)
ZC.p <- ZC(pf)
for(i in 1:length(classes)) {
  lines(1:5, ZC.p[(1:5)+5*(i-1)], col=col[i], lwd=lwd[i])
  if(classes[i] %in% clab) text(0.8, ZC.p[1+5*(i-1)], classes[i], adj=1)
}


###################################################
### code chunk number 8: basis
###################################################
basis(c("HCO3-", "H2O", "NH3", "HS-", "H2", "H+"))


###################################################
### code chunk number 9: TpH
###################################################
distance <- c(0, 6, 11, 14, 22)
T.bison <- c(93.3, 79.4, 67.5, 65.3, 57.1)
pH.bison <- c(7.350, 7.678, 7.933, 7.995, 8.257)


###################################################
### code chunk number 10: TpHplot
###################################################
Tfun <- splinefun(distance, T.bison, method="mono")
pHfun <- splinefun(distance, pH.bison, method="mono")
xpoints <- seq(0, 22, 0.25)
par(mfrow=c(1, 2), mar=c(4, 4, 3, 2))
plot(distance, T.bison,xlab="distance, m", ylab=axis.label("T"))
lines(xpoints, Tfun(xpoints))
plot(distance, pH.bison,xlab="distance, m", ylab="pH")
lines(xpoints, pHfun(xpoints))


###################################################
### code chunk number 11: basis
###################################################
basis(c("HCO3-", "NH3", "HS-", "H+"), c(-3, -4, -7, -7.933))


###################################################
### code chunk number 12: formation_protein
###################################################
species("overall", sitenames)


###################################################
### code chunk number 13: formation_residue
###################################################
ip <- iprotein("overall", sitenames)
pl <- protein.length(ip)
mys <- species()
mys[, 1:6]/pl


###################################################
### code chunk number 14: T_H2
###################################################
Tlim <- c(50, 100)
H2prot <- function(T) -11 + T * 3/40


###################################################
### code chunk number 15: TlogaH2plot
###################################################
a <- affinity(T=Tlim, H2=c(-7, -4))
diagram(a, fill=NULL, names=1:5, normalize=TRUE)
lines(Tlim, H2prot(Tlim), lty=2)
title(main="overall")


###################################################
### code chunk number 16: loadclass
###################################################
loadclass <- function(class) {
  species(delete=TRUE)
  species(rep(class, each=5), rep(sitenames, length(class)))
}


###################################################
### code chunk number 17: TlogaH2plots
###################################################
forclasses <- c("transferase", "transport", "dehydrogenase", "synthase")
par(mfrow=c(2, 2))
for(i in 1:4) {
  loadclass(forclasses[i])
  a <- affinity(T=Tlim, H2=c(-7, -4))
  diagram(a, fill=NULL, names=1:5, normalize=TRUE)
  title(main=forclasses[i])
  lines(Tlim, H2prot(Tlim), lty=2)
}
loadclass("overall")


###################################################
### code chunk number 18: affinity
###################################################
species(1:5,0)
a <- affinity(T=T.bison, pH=pH.bison, H2=H2prot(T.bison))
a.res <- t(as.data.frame(a$values))/pl
a.res


###################################################
### code chunk number 19: affinitymax
###################################################
sapply(1:5, function(x) which.max(a.res[, x]))


###################################################
### code chunk number 20: affinityadj
###################################################
a.res <- a.res - log10(pl)
a.res
sapply(1:5, function(x) which.max(a.res[, x]))


###################################################
### code chunk number 21: logaplot
###################################################
species(1:5, -3)
xT <- Tfun(xpoints)
xpH <- pHfun(xpoints)
xH2 <- H2prot(xT)
a <- affinity(T=xT, pH=xpH, H2=xH2)


###################################################
### code chunk number 22: logaplot
###################################################
a$vars[1] <- "distance, m"
a$vals[[1]] <- xpoints
e <- equilibrate(a, normalize=TRUE)
diagram(e, legend.x=NULL)
legend("bottom", lty=1:5, legend=paste(1:5,"       "), bty="n")


###################################################
### code chunk number 23: stripplot
###################################################
xclasses <- c("overall", "transferase", "transport", "synthetase", "membrane", "permease")
loadclass(xclasses)
a <- affinity(T=xT, pH=xpH, H2=xH2)
a$vars[1] <- "distance, m"
a$vals[[1]] <- xpoints
col <- c("red", "orange", "yellow", "green", "blue")
par(mfrow=c(1, 2), mar=c(4, 4, 1, 1))
for(i in 1:2) {
  ispecies <- lapply((1:3)+(i-1)*3, function(x) {1:5+(x-1)*5} )
  names(ispecies) <- xclasses[(1:3)+(i-1)*3]
  strip(a = a, ispecies = ispecies, col = col, xticks = distance, cex.names = 1)
}


###################################################
### code chunk number 24: E.AgAgCl
###################################################
E.AgAgCl <- function(T) {
  0.23737 - 5.3783e-4 * T - 2.3728e-6 * T^2 - 2.2671e-9 * (T+273)
}


###################################################
### code chunk number 25: meters
###################################################
T.ORP <- c(93.9, 87.7, 75.7, 70.1, 66.4, 66.2)
pH.ORP <- c(8.28, 8.31, 7.82, 7.96, 8.76, 8.06)
ORP <- c(-258, -227, -55, -58, -98, -41)


###################################################
### code chunk number 26: ORP2Eh
###################################################
Eh <- ORP/1000 + E.AgAgCl(T.ORP)
pe <- convert(Eh, "pe", T=convert(T.ORP, "K"))
logK.ORP <- subcrt(c("e-", "H+", "H2"), c(-2, -2, 1), c("aq", "aq", "aq"), T=T.ORP)$out$logK
logaH2.ORP <- logK.ORP - 2*pe - 2*pH.ORP


###################################################
### code chunk number 27: sulfur
###################################################
loga.HS <- log10(c(4.77e-6, 2.03e-6, 3.12e-7, 4.68e-7, 2.18e-7))
loga.SO4 <- log10(c(2.10e-4, 2.03e-4, 1.98e-4, 2.01e-4, 1.89e-4))
logK.S <- subcrt(c("HS-", "H2O", "SO4-2", "H+", "H2"), c(-1, -4, 1, 1, 4), 
  state=c("aq", "liq", "aq", "aq", "aq"), T=T.bison)$out$logK
logaH2.S <- (logK.S + pH.bison - loga.SO4 + loga.HS) / 4


###################################################
### code chunk number 28: oxygen
###################################################
DO <- c(0.173, 0.776, 0.9, 1.6, 2.8)
logaO2 <- log10(DO/1000/32)
logK <- subcrt(c("O2", "H2", "H2O"), c(-0.5, -1, 1), c("aq", "aq", "liq"), T=T.bison)$out$logK
logaH2.O <- 0 - 0.5*logaO2 - logK


###################################################
### code chunk number 29: logaH2plot
###################################################
plot(Tlim, H2prot(Tlim), xlim=Tlim, ylim=c(-45,0),
  xlab=axis.label("T"), ylab=axis.label("H2"), type="l")
points(T.ORP, logaH2.ORP, pch=15)
lines(T.ORP, logaH2.ORP, lty=2)
points(T.bison, logaH2.O, pch=16)
lines(T.bison, logaH2.O, lty=2)
points(T.bison, logaH2.S, pch=17)
lines(T.bison, logaH2.S, lty=2)
llab <- c("H2prot", "ORP", "oxygen", "sulfur")
text(c(72, 80, 72, 74), c(-4, -25, -34, -10.5), llab)
# legend("bottomright", lty=c(1, NA, NA, NA), pch=c(NA, 15, 16, 17), legend=llab)


