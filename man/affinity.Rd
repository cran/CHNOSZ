\encoding{UTF-8}
\name{affinity}
\alias{affinity}
\title{Chemical Affinities of Formation Reactions}
\description{
  Calculate the chemical affinities of formation reactions of species. Do it for single values of temperature, pressure, ionic strength and chemical activities of the basis species, or as a function of one or more of these variables. Or, return other properties including standard molal Gibbs energies of basis species and species of interest, and standard molal Gibbs energies, equilibrium constants and activity products of formation reactions. 
}

\usage{
  affinity(..., property=NULL, sout=NULL, exceed.Ttr=FALSE,
    return.buffer=FALSE, balance="PBB", iprotein=NULL, loga.protein=-3)
}

\arguments{
  \item{...}{numeric, zero or more named arguments, used to identify the variables of interest in the calculations}
  \item{property}{character, denoting the property to be calculated. Default is \samp{A}, for chemical affinity of formation reactions of species of interest}
  \item{sout}{list, output from subcrt function}
  \item{exceed.Ttr}{logical, allow subcrt to compute properties for phases beyond their transition temperature?}
  \item{return.buffer}{logical. If \code{TRUE}, and a \code{\link{buffer}} has been associated with one or more basis species in the system, return the values of the activities of the basis species calculated using the buffer (it is not necessary in this case to have defined any species of interest). Default is \code{FALSE}}
  \item{balance}{character. This argument is passed to \code{\link{buffer}} to identify a conserved basis species (or \samp{PBB}) in a chemical activity buffer. Default is \samp{PBB}}
  \item{iprotein}{numeric, indices of proteins in \code{\link{thermo}$protein} for which to calculate properties}
  \item{loga.protein}{numeric, logarithms of activities of proteins identified in \code{iprotein}}
}

\details{
\code{affinity} calculates the chemical affinities of reactions to form the species of interest from the basis species.
The calculation of chemical affinities relies on the current definitions of the \code{\link{basis}} species and \code{\link{species}} of interest.
It is possible to use the results of \code{affinity} to generate equilibrium activity diagrams using \code{\link{diagram}}.

The equation used to calculate chemical affinity \emph{\bold{A}} is \emph{\bold{A}}=\eqn{RT\ln (K/Q)}{RT*ln(K/Q)} (Kondepudi and Prigogine, 1998), where \eqn{K} denotes the equilibrium constant of the reaction and \eqn{Q} stands for the activity product of the species in the reaction.
(The equilibrium constant is related to standard Gibbs energy of reaction, \eqn{\Delta G^{\circ}_r}{DeltaG0r}, by \eqn{\Delta G^{\circ}_r = -2.303RT\log K}{DeltaG0r = -2.303*RT*logK}, where \eqn{R} and \eqn{T} stand for, respectively, the gas constant and temperature).
With the approach of a given reaction to a state of equilibrium, the chemical affinity tends toward zero, or \eqn{K = Q}. 

Valid properties are \samp{A} or NULL for chemical affinity, \samp{logK} or \samp{logQ} for logarithm of equilibrium constant and reaction activity product, or any of the properties available in \code{\link{subcrt}} except for \samp{rho}.
The properties returned are those of the formation reactions of the species of interest from the basis species.
It is also possible to calculate the properties of the species of interest themselves (not their formation reactions) by setting the \code{property} to \samp{G.species}, \samp{Cp.species}, etc.
Except for \samp{A}, the properties of proteins or their reactions calculated in this manner are restricted to nonionized proteins.

Zero, one, or more leading arguments to the function identify which of the chemical activities of basis species, temperature, pressure and/or ionic strength to vary.
The names of each of these arguments may be the formula of any of the basis species of the system, or \samp{T}, \samp{P}, \samp{pe}, \samp{pH}, \samp{Eh}, or \samp{IS} (but names may not be repeated).
To use the names of charged basis species such as \samp{K+} and \samp{SO4-2} as the arguments, they should be enclosed in quotes (see the example for aluminum speciation in \code{\link{diagram}}).
The value of each argument is of the form \code{c(min,max)} or \code{c(min,max,res)} where \code{min} and \code{max} refer to the minimimum and maximum values of variable identified by the name of the argument, and \code{res} denotes the resolution, or number of points along which to do the calculations (this is assigned a default value of 128 if it is missing).
For any arguments that refer to basis species, the numerical values are the logarithms of the activities of that basis species, or logarithms of fugacities if it is a gas.
Unlike the \code{energy} function, the units of \samp{T} and \samp{P} in \code{affinity} are those the user has set using \code{\link{T.units}} and \code{\link{P.units} }(on program start-up these are \eqn{^{\circ}}{°}C and bar, respectively). 

If one or more buffers are assigned to the definition of \code{\link{basis}} species, \code{affinity} calls \code{\link{buffer}} to calculate the logarithms of activities of these basis species from the buffer.

The \code{iprotein} and \code{loga.protein} arguments can be used to compute the chemical affinities of formation reactions of proteins that are not in the current \code{\link{species}} definition.
This approach can be utilized in order to calculate the properties of many proteins in a fraction of the time it would take to calculate them individually.
The appropriate \code{\link{basis}} species still must be defined prior to calling \code{affinity}.
\code{iprotein} contains indices of desired proteins in \code{\link{thermo}$protein}; \code{affinity} adds to the species list the amino acid residues and and terminal H2O group (indicated by \dQuote{RESIDUE} in \code{thermo$protein}) then calculates the properties of the reactions for the residues and terminal group, including ionization effects, and adds them together to get those of the indicated proteins.

In CHNOSZ version 0.9, \code{energy} gained a new argument \samp{transect} which is set to TRUE by \code{energy.args} when the length(s) of the variables is(are) greater than three.
In this mode of operation, instead of performing the calculations on an \eqn{n}{n}-dimensional grid, the affinities are calculated on an \eqn{n}{n}-dimensional transect through chemical potential (possibly including T and/or P) space. 

}

\value{
For \code{affinity}, a list, elements of which are \code{sout} output from \code{\link{subcrt}}, \code{property} name of the calculated property (\samp{A} for chemical affinity), \code{basis} and \code{species} definition of basis species and species of interest in effect at runtime, \code{T} and \code{P} temperature and pressure, in the system units of Kelvin and bar, set to \code{numeric()} (length=0) if either one is a variable, \code{vars} the names of the variables, \code{vals} the values of the variables (a list, one element for each variable), \code{values} the result of the calculation (a list, one element for each species, with names taken from the species index in \code{\link{thermo}$obigt}).
The elements of the lists in \code{vals} and \code{values} are arrays of \eqn{n}{n} dimensions, where \eqn{n}{n} is the number of variables.
The values of chemical affinity of formation reactions of the species are returned in dimensionless units (for use with decimal logarithms, i.e., \emph{\strong{A}}/\eqn{2.303RT}).

Names other than \samp{T} or \samp{P} in \code{vars} generally refer to basis species, and the corresponding \code{vals} are the logarithms of activity or fugacity. However, if one or more of \code{pe}, \code{Eh} or \code{pH} is among the variables of interest, \code{vals} holds the values of the those variables as indicated.
}

\seealso{
Usually, \code{\link{equilibrate}} is the next step in calculations of chemical equilibrium.
The help for \code{\link{buffer}} has some examples of using chemical activity buffers.
The guts of the calculations provided by \code{affinity} involve \code{\link{energy.args}} and \code{\link{energy}}, which normally are not part of the user interaction.
}

\examples{
\dontshow{data(thermo)}## set up a system and calculate
## chemical affinities of formation reactions
basis(c("SiO2","MgO","H2O","O2"),c(-5,-5,0,999))
species(c("quartz","enstatite","forsterite"))
# chemical affinities (A/2.303RT) at 25 deg C and 1 bar
affinity()
# at higher temperature and pressure
affinity(T=500,P=2000)
# ten different temperatures at one pressure
affinity(T=c(500,1000,10),P=2000)
# at 25 temperatures and pressures
affinity(T=c(500,1000,5),P=c(1000,5000,5))
# as a function of logarithm of activity of MgO
affinity(MgO=c(-10,-5,10))
## equilibrium constants of formation reactions
affinity(property="logK")
# Standard molal Gibbs energies of species,
# user units (default: cal/mol)
affinity(property="G.species")
# Standard molal Gibbs energies of reactions
affinity(property="G")

## affinities of metabolic reactions
## after Amend and Shock, 2001, Fig. 7
# use aq state for all basis species (including O2)
basis(c("CO2", "H2", "NH3", "O2", "H2S", "H+"), "aq")
# we're going to make H2O
species("H2O")
# a function to create the plots
doplot <- function(T) {
  res <- 20
  # calculate affinity/2.303RT as a function of loga(H2) and loga(O2)
  a <- affinity(H2=c(-10, 0, res), O2=c(-10, 0, res), T=T)
  T.K <- convert(T, "K")                   # temperature in Kelvin
  acal <- convert(a$values[[1]], "G", T.K) # affinity (cal/mol)
  akJ <- convert(acal, "J")/1000           # affinity (kJ/mol)
  # now contour the values
  xyvals <- seq(-10, 0, length.out=res)
  contour(x=xyvals, y=xyvals, z=t(akJ), levels=seq(-150, -250, -20),
    labcex=1, xlab=axis.label("H2"), ylab=axis.label("O2"))
  # show the temperature
  legend("topleft", bg="white", cex=1,
    legend=describe.property("T", T, digits=0, ret.val=TRUE) )
}
# plot layout with space for title at top
layout(matrix(c(1, 1, 2, 3, 4, 5), ncol=2, byrow=TRUE), heights=c(1, 4, 4))
par(mar=c(0, 0, 0, 0))
plot.new()
# we use subcrt() to generate a reaction for titling the plot
rxnexpr <- describe.reaction(subcrt("H2O", 1)$reaction, states="all")
# also in the title is the property with its units
E.units("J")
Gexpr <- axis.label("DGr", prefix="k")[[2]]
text(0.5, 0.6, substitute(paste(G~~"for"~~r), list(G=Gexpr, r=rxnexpr)), cex=2)
text(0.5, 0.2, "after Amend and Shock, 2001 Figure 7", cex=2)
# now make the plots
par(mar=c(3, 3, 0.5, 0.5), cex=1.3, mgp=c(2, 1, 0))
sapply(c(25, 55, 100, 150), doplot)
# affinity() can handle the three dimensions simultaneously
print(affinity(H2=c(-10, 0, 3), O2=c(-10, 0, 3), T=c(25, 150, 4))$values)
# this is so the plots in the next examples show up OK
E.units("cal")
layout(matrix(1))
par(mar=c(5.1, 4.1, 4.1, 2.1))

## amino acid synthesis at low and high temperatures
## after Amend and Shock, 1998
# select the basis species and species of interest
# and set their activities, first for the 18 degree C case
basis(c("H2O", "CO2", "NH4+", "H2", "H+", "H2S"),
  log10(c(1, 1e-4, 5e-8, 2e-9, 5e-9, 1e-15)))
species(sort(aminoacids("Z")),
  log10(c(3.9, 0.7, 1.1, 3.3, 0.5, 3.8, 1.0, 5.8, 1.2, 0.7,
  0.8, 1.0, 2.8, 0.5, 0.5, 4.6, 5.8, 0.6, 0.9, 2.8)/1e9))
T <- 18
TK <- convert(T, "K")
# calculate A/2.303RT (dimensionless), convert to G of reaction (cal/mol)
a <- affinity(T=T)
G.18.cal <- convert(unlist(a$values), "G", T=TK)
# covvert to kJ/mol
G.18.kJ <- convert(G.18.cal, "J")/1000
# the 100 degree C case
basis(c("H2O", "CO2", "NH4+", "H2", "H+", "H2S"),
  log10(c(1, 2.2e-3, 2.9e-6, 3.4e-4, 1.9e-6, 1.6e-3)))
species(1:20, log10(c(2.8e-9, 5.0e-10, 7.9e-10, 2.4e-9, 3.6e-10,
  2.7e-9, 7.2e-10, 4.2e-9, 8.6e-10, 5.0e-10, 5.7e-10, 7.2e-10, 2.0e-9,
  3.6e-10,3.6e-10, 3.3e-9, 4.2e-9, 4.3e-10, 6.5e-10, 2.0e-9)))
T <- 100
TK <- convert(T, "K")
a <- affinity(T=T)
G.100.cal <- convert(unlist(a$values), "G", T=TK)
G.100.kJ <- convert(G.100.cal, "J")/1000
# the average oxidation states of carbon
Z.C <- ZC(thermo$obigt$formula[thermo$species$ispecies])
# put everything together a la Table 3 in the paper
print(out <- data.frame(G.18=G.18.kJ, G.100=G.100.kJ, Z.C=Z.C))
# make a plot; set units to get correct label
E.units("J")
plot(out$Z.C, out$G.18, pch=20, xlim=c(-1.1, 1.1), ylim=c(-200, 500), 
  xlab=axis.label("ZC"), ylab=axis.label("DGr"))
points(out$Z.C, out$G.100, col="red", pch=20)
legend("topleft", pch=c(20, 20), col=c("black", "red"),
  legend=describe.property(c("T", "T"), c(18, 100)))
title(main="Amino acid synthesis, after Amend and Shock, 1998")
# 9 amino acids have negative delta Gr under hydrothermal conditions
# (cf. AS98 with 11; we are using more recent thermodynamic data)
stopifnot(sum(out$G.100 < 0)==9)
# reset units and species to run next examples
E.units("cal")
species(delete=TRUE)

## calculations along a transect: methanogenesis and biosynthetic 
## reactions in hydrothermal systems, after Shock and Canovas, 2010
# this file has their mixing path results for Rainbow hydrothermal field
file <- system.file("extdata/cpetc/SC10_Rainbow.csv", package="CHNOSZ")
rb <- read.csv(file, check.names=FALSE)
# write all synthesis reactions in terms of these basis species
# it's okay not to set the activities of the basis species now
# because they'll be changing along with temperature
basis(c("CO2", "H2", "NH4+", "H2O", "H2S", "H+"))
# now a selection of the species from SC10, with activities equal to 1e-6
species(c("CH4", "formaldehyde", "ethylene", "glycolic acid", 
  "n-nonanoic acid", "leucine", "aspartic acid", "tryptophan", "deoxyribose", 
  "adenine", "cytosine"), -6)
# the exception is methane; unlike SC10 we use a constant activity 1e-3
# (accounting for variable activities of the species of interest here 
# is possible but would require longer code ....)
species("CH4", -3)
# synchronized change of temperature and five basis activities
a <- affinity(T=rb$T, CO2=rb$CO2, H2=rb$H2, `NH4+`=rb$`NH4+`, H2S=rb$H2S, pH=rb$pH)
# the tricky part: affinity() uses dimensionless values (A/2.303RT)
# but we want to show the values in cal/mol
a$values <- lapply(a$values, function(val) {
  -convert(val, "G", T=convert(a$vals[[1]], "K")) })
# if we didn't have balance=1 here the values would be 
# divided by the number of moles of CO2 in the reactions ...
diagram(a, balance=1, ylim=c(-100000, 100000), ylab=axis.label("A"),
  col=topo.colors(4), lwd=2)
# add a zero-affinity line and a title
abline(h=0, lty=2, lwd=2)
title(main="Affinities of organic synthesis, after Shock and Canovas, 2010")
}

\references{
  Amend, J. P. and Shock, E. L. (1998) Energetics of amino acid synthesis in hydrothermal ecosystems. \emph{Science} \bold{281}, 1659--1662. \url{http://dx.doi.org/10.1126/science.281.5383.1659}

  Amend, J. P. and Shock, E. L. (2001) Energetics of overall metabolic reactions of thermophilic and hyperthermophilic Archaea and Bacteria. \emph{FEMS Microbiol. Rev.} \bold{25}, 175--243. \url{http://dx.doi.org/10.1016/S0168-6445(00)00062-0}

  Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \url{http://www.biogeosciences.net/3/311/2006/bg-3-311-2006.html}

  Kondepudi, D. K. and Prigogine, I. (1998) \emph{Modern Thermodynamics: From Heat Engines to Dissipative Structures}, John Wiley & Sons, New York, 486 p. \url{http://www.worldcat.org/oclc/38055900}

  Shock, E. and Canovas, P. (2010) The potential for abiotic organic synthesis and biosynthesis at seafloor hydrothermal systems. \emph{Geofluids} \bold{10}, 161--192. \url{http://dx.doi.org/10.1111/j.1468-8123.2010.00277.x}
}


\concept{Primary thermodynamic modeling}
