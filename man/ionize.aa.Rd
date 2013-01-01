\encoding{UTF-8}
\name{ionize.aa}
\alias{ionize.aa}
\title{Properties of Ionization of Proteins}
\description{
  Calculate the charges of proteins and contributions of ionization to the thermodynamic properties of proteins.
}

\usage{
  ionize.aa(aa, property = "Z", T = 25, P = "Psat", pH = 7,
    ret.val = NULL, suppress.Cys = FALSE)
}

\arguments{
  \item{aa}{data frame, amino acid composition in the format of \code{thermo$protein}}
  \item{property}{character, property to calculate}
  \item{T}{numeric, temperature in \eqn{^{\circ}}{°}C}
  \item{P}{numeric, pressure in bar, or \samp{Psat} for vapor pressure of \eqn{\mathrm{H_2O}}{H2O} above 100 \eqn{^{\circ}}{°}C}
  \item{pH}{numeric, pH}
  \item{ret.val}{character, return the indicated value from intermediate calculations}
  \item{suppress.Cys}{logical, suppress (ignore) the ionization of the cysteine groups?}
}

\details{
  The properties of ionization of proteins calculated by this function take account of the standard molal thermodynamic properties of ionizable amino acid sidechain groups and the terminal groups in proteins ([AABB]) and their equations of state parameters taken from Dick et al., 2006. The values of the ionization constants (pK) are calculated as a function of temperature, and the charges and the ionization contributions of other thermodynamic properties to the proteins are calculated additively, without consideration of electrostatic interactions, so they are best applied to the unfolded protein reference state.

  For each amino acid composition in \code{aa}, the additive value of the \code{property} is calculated as a function of \code{T}, \code{P} and \code{pH}. \code{property} can be NULL to denote net charge, or if not NULL is one of the properties available in \code{\link{subcrt}}, or is \samp{A} to calculate the dimensionless chemical affinity (A/2.303RT) of the ionization reaction for the protein. If \code{ret.val} is one of \samp{pK}, \samp{alpha}, or \samp{aavals} it indicates to return the value of the ionization constant, degree of formation, or the values of the \code{property} for each ionizable group rather than taking their sums for the amino acid compositions in \code{aa}.
}

\value{
  The function returns a matrix (possibly with only one row or column) with number of rows corresponding to the longest of \code{T}, \code{P} or \code{pH} (values of any of these with shorter length are recycled) and a column for each of the amino acid compositions in \code{aa}.
}

\seealso{
  The amino acid composition in \code{aa} can be generated using e.g. \code{\link{ip2aa}}. 

  \code{\link{protein.info}} and \code{\link{protein.basis}} use this function to compute the properties of ionized proteins.

  This function is called by \code{\link{A.ionization}} as part of a calculation of \code{\link{affinity}} if proteins are among the species of interest, \samp{H+} is in the basis, and \code{thermo$opt$ionize} is TRUE. 
}

\examples{\dontshow{data(thermo)}
## Charge of LYSC_CHICK as a function of pH and T
# After Fig. 10 of Dick et al., 2006
# the rownumber of the protein in thermo$protein
ip <- iprotein("LYSC_CHICK")
# its amino acid composition
aa <- ip2aa(ip)
# additive charges of unfolded protein at 25, 100, 150 degrees C
# as a function of pH
pH <- seq(0, 14, 0.1)
Z.25 <- ionize.aa(aa, T=25, pH=pH)
plot(pH, Z.25[, 1], type="l", xlab="pH", ylab="net charge (Z)")
lines(pH, ionize.aa(aa, T=100, pH=pH)[, 1], col="red")
lines(pH, ionize.aa(aa, T=150, pH=pH)[, 1], col="orange")
text(c(12, 10, 9), c(-15, -16, -18),
  labels=paste("T=", c(25, 100, 150), sep=""))
# suppress ionization of cysteine as if it was oxidized 
# (to form non-ionizable cystine disulfide bonds)
lines(pH, ionize.aa(aa, T=25, pH=pH, suppress.Cys=TRUE)[, 1], lty=2)
text(12, -7, "T=25, oxidized")
# add experimental points
RT71 <- read.csv(system.file("extdata/cpetc/RT71.csv", package="CHNOSZ"))
points(RT71$pH, RT71$Z)
legend("topright", pch=1, legend="Roxby and Tanford, 1971")
title(main=paste("Ionization of unfolded LYSC_CHICK\n",
  "After Dick et al., 2006"))

\donttest{## Heat capacity of LYSC_CHICK as a function of T
pH <- c(5, 9, 3)
T <- seq(0, 100)
# Cp of non-ionized protein
Cp.nonion <- subcrt("LYSC_CHICK", T=T)$out[[1]]$Cp
plot(T, Cp.nonion, xlab=axis.label("T"), type="l",
  ylab=axis.label("Cp"), ylim=c(5000, 8000))
# Cp of ionization and ionized protein
aa <- ip2aa("LYSC_CHICK")
for(pH in c(5, 9, 3)) {
  Cp.ionized <- Cp.nonion + ionize.aa(aa, "Cp", T=T, pH=pH)[, 1]
  lines(T, Cp.ionized, lty=2)
  text(80, Cp.ionized[70], paste("pH =",pH) )
}
# Makhatadze and Privalov's group contributions
T <- c(5, 25, 50, 75, 100, 125)
points(T, convert(MP90.cp("LYSC_CHICK", T), "cal"))
# Privalov and Makhatadze's experimental values
e <- read.csv(system.file("extdata/cpetc/PM90.csv", package="CHNOSZ"))
e <- e[e$protein=="LYSC_CHICK",]
points(e$T, convert(e$Cp, "cal"), pch=16)
legend("bottomright", pch=c(16, 1, NA, NA), lty=c(NA, NA, 1, 2),
  legend=c("PM90 experiment", "MP90 groups", 
  "DLH06 groups no ion", "DLH06 groups ionized"))
title("Heat capacity of unfolded LYSC_CHICK")

## Contour plots of net charge and ionization properties of LYSC_CHICK
\dontshow{par(thermo$opar)}
aa <- ip2aa("LYSC_CHICK")
pH <- seq(0, 14, 0.2)
T <- seq(0, 200, 2)
val <- expand.grid(pH=pH, T=T)
par(mfrow=c(2, 2))
for(X in c("Z", "A", "Cp", "V")) {
  Y <- ionize.aa(aa, property=X,  pH=val$pH, T=val$T)
  contour(pH, T, matrix(Y[, 1], ncol=length(T)),
    xlab="pH", ylab=axis.label("T"))
  title(main=axis.label(X))
}
par(mfrow=c(1, 1))
pu <- par("usr")
text(mean(pu[1:2]), sum(pu[3:4])*0.45, 
  "additive properties of ionization of LYSC_CHICK")
\dontshow{par(thermo$opar)}
}
}

\references{
  Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \url{http://www.biogeosciences.net/3/311/2006/bg-3-311-2006.html}

  Makhatadze, G. I. and Privalov, P. L. (1990) Heat capacity of proteins. 1. Partial molar heat capacity of individual amino acid residues in aqueous solution: Hydration effect. \emph{J. Mol. Biol.} \bold{213}, 375--384. \url{http://dx.doi.org/10.1016/S0022-2836(05)80197-4}

  Privalov, P. L. and Makhatadze, G. I. (1990) Heat capacity of proteins. II. Partial molar heat capacity of the unfolded polypeptide chain of proteins: Protein unfolding effects. \emph{J. Mol. Biol.} \bold{213}, 385--391. \url{http://dx.doi.org/10.1016/S0022-2836(05)80198-6}
}

\keyword{protein}
