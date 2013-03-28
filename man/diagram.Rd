\name{diagram}
\alias{diagram}
\alias{strip}
\alias{find.tp}
\title{Equilibrium Chemical Activity Diagrams}
\description{
Plot equilibrium chemical activity (1-D speciation) or equal-activity (2-D predominance) diagrams as a function of chemical activities of basis speecies, temperature and/or pressure.
}

\usage{
  diagram(eout, what = "loga.equil", alpha = FALSE, normalize = FALSE,
    balance=NULL, groups=as.list(1:length(eout$values)), xrange=NULL,
    mar=NULL, yline=par("mgp")[1]+0.7, side=1:4,
    ylog=TRUE, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
    cex=par("cex"), cex.names=1, cex.axis=par("cex"),
    lty=NULL, lwd=par("lwd"), dotted=0,
    bg=par("bg"), col=par("col"), col.names=par("col"), fill=NULL,
    names=NULL, main=NULL, legend.x="topright", add=FALSE, plot.it=TRUE)
  strip(affinity, ispecies = NULL, col = NULL, ns = NULL, 
    xticks = NULL, ymin = -0.2, xpad = 1, cex.names = 0.7)
  find.tp(x)
}

\arguments{
  \item{eout}{list, object returned by \code{\link{equilibrate}}} or \code{\link{affinity}}
  \item{what}{character, what property to calculate and plot}
  \item{alpha}{logical, for speciation diagrams, plot degree of formation instead of activities?}
  \item{normalize}{logical, normalize chemical formulas or reaction properties by the balance vector?}
  \item{balance}{character, balancing constraint; see \code{\link{balance}}}
  \item{groups}{list of numeric, groups of species to consider as a single effective species}
  \item{xrange}{numeric, range of x-values between which predominance field boundaries are plotted}
  \item{mar}{numeric, margins of plot frame}
  \item{yline}{numeric, margin line on which to plot the y-axis name}
  \item{side}{numeric, which sides of plot to draw axes}
  \item{xlim}{numeric, limits of x-axis}
  \item{ylim}{numeric, limits of y-axis}
  \item{xlab}{character, label to use for x-axis}
  \item{ylab}{character, label to use for y-axis}
  \item{ylog}{logical, use a logarithmic y-axis (on 1D degree diagrams)?}
  \item{cex}{numeric, character expansion (scaling relative to current)}
  \item{cex.names}{numeric, character expansion factor to be used for names of species on plots}
  \item{cex.axis}{numeric, character expansion factor for names of axes}
  \item{lty}{numeric, line types to be used in plots}
  \item{lwd}{numeric, line width}
  \item{dotted}{numeric, how often to skip plotting points on predominance field boundaries (to gain the effect of dotted or dashed boundary lines)}
  \item{bg}{character, background color for legend}
  \item{col}{character, color of activity lines (1D diagram) or predominance field boundaries (2D diagram), or colors of bars in a strip diagram (\code{strip})}
  \item{col.names}{character, colors for labels of species}
  \item{fill}{character, colors used to fill predominance fields}
  \item{names}{character, names of species for activity lines or predominance fields}
  \item{main}{character, a main \code{\link{title}} for the plot; \code{NULL} means to plot no title}
  \item{legend.x}{character, description of legend placement passed to \code{\link{legend}}}
  \item{add}{logical, add to current plot?}
  \item{plot.it}{logical, make a plot?}
  \item{affinity}{list, object returned by \code{\link{affinity}}}
  \item{ispecies}{numeric, which species to consider (default of \code{NULL} is to consider all species)}
  \item{ns}{numeric, numbers of species, used to make inset plots for strip diagrams}
  \item{xticks}{numeric, location of supplemental tick marks on x-axis}
  \item{ymin}{numeric, lower limit of y-axis}
  \item{xpad}{numeric, amount to extend x-axis on each side}
  \item{x}{matrix, value of the \code{predominant} list element from \code{diagram}}
}

\details{

\code{diagram} takes as its primary input the results from \code{\link{equilibrate}} and displays diagrams representing the equilibrium chemical activities of the species.
0-D diagrams, at a single point, are shown as barcharts.
1-D diagrams, for a single variable on the x-axis, are plotted as lines.
2-D diagrams, for two variables, are plotted as predominance fields.
The allowed variables are any that \code{\link{affinity}} accepts: temperature, pressure, or the chemical activities of the basis species

If \code{alpha} is TRUE, the fractional degrees of formation (ratios of activities to total activity) are plotted.
This setting is useful for visualizing the molalities (activities for ideal interactions) of monomer groups in a system of biopolymers (e.g. proteins).
If \code{groups} is supplied, the activities of the species identified in each numeric vector of this list are summed together and subsequently treated as a single species; the names of the new species are taken from the names of this list.

For 1-D diagrams, the default setting for the y-axis is a logarithmic scale (unless \code{alpha} is TRUE) with limits corresponding to the range of logarithms of activities (or 0,1 if \code{alpha} is TRUE); these actions can be overridden by \code{ylog} and \code{ylim}.
A \code{\link{legend}} is placed at the location identified by \code{legend.x}, or omitted if \code{legend.x} is \code{FALSE}.
If \code{legend.x} is NA, instead of a legend, the lines are labeled with the names of the species near the maximum value.
The line type and line width can be controlled with \code{lty} and \code{lwd}, respectively.

On 2-D diagrams, the fields represent the species with the highest equilibrium activity, after any \code{normalize} or \code{alpha} operations.
\code{fill} determines the color of the predominance fields, \code{col} that of the boundary lines.
By default, \code{\link{heat.colors}} are used to fill the predominance fields in diagrams on the screen plot device.
The style of the boundary lines on 2-D diagrams can be altered by supplying one or more non-zero integers in \code{dotted}, which indicates the fraction of line segments to omit; a value of \samp{1} or NULL for \code{dotted} has the effect of not drawing the boundary lines.

For all diagrams, the \code{names} of the species and their colors in \code{col.names} can be changed, as can \code{cex}, \code{cex.names}, and \code{cex.axis} to adjust the overall character expansion factors (see \code{\link{par}}) and those of the species names and axis labels.
The x- and y-axis labels are automatically generated unless they are supplied in \code{xlab} and \code{ylab}. 
A new plot is started unless \code{add} is TRUE.
If \code{plot.it} is FALSE, no plot will be generated but all the intermediate computations will be performed and the results returned.

\code{diagram} also accepts the output from \code{\link{affinity}}, for which three actions are possible:
1) plot a property of a reaction, such as the equilibrium constant (\samp{logK}) (0-D or 1-D);
2) plot the equilibrium activity of a selected basis species in all of the formation reactions (0-D, 1-D or 2-D);
3) plot predominance fields, based on the relative magnitudes of the affinities of the formation reactions (2-D only).

Some of the arguments have different effects when the output from \code{affinity} is being used instead of the equilibrium activities from \code{equilibrate}.
If \code{what} is missing, option (1) is selected if the number of dimensions is 0 or 1, and option (3) is selected if the number of dimensions is 2.
The latter is referred to as the maximum affinity method.
In cases where it applies (see Warning), the maximum affinity method is much faster than an equilibrium calculation.
\code{balance} is the option, sent to \code{\link{balance}}, that determines the balancing coefficients used in the normalization (this argument has no effect on the calculations of equilibrium activities.)
If \code{what} is the name of a basis species, it refers to option (2) above.
A contour plot is made in the case of 2-D diagrams of the equilibrium activity of a basis species (see \code{\link{demos}("CO2Ac")}, and only the first species of interest is used in the calculation; a warning is produced if there is more than one.

A different incarnation of 1-D speciation diagrams is provided by \code{strip}.
This function generates any number of strip diagrams in a single plot.
The diagrams are made up of colors bars whose heights represent the relative abundances of species; the color bars are arranged in order of abundance and the total height of the stack of colors bars is constant.
If \code{ispecies} is a list, the number of strip diagrams is equal to the number of elements of the list, and the elements of this list are numeric vectors that identify the species to consider for each diagram.
The strips are labeled with the \code{\link{names}} of \code{ispecies}.
If \code{col} is NULL, the colors of the bars are generated using \code{\link{rainbow}}.
Supplemental ticks can be added to the x-axis at the locations specified in \code{xtick}; they are larger than the standard ticks and have colors corresponding to those of the color bars.
\code{ymin} can be decreased in order to add more space at the bottom of the plot, and \code{xpad} can be changed in order to increase or decrease the size of the x-axis relative to the width of the strips.
An inset dot-and-line plot is created below each strip if \code{ns} is given.
This argument has the same format as \code{ispecies}, and can be used e.g. to display the relative numbers of species for comparison with the stability calculations.

\code{find.tp} finds the locations in a matrix of integers that are surrounded by the greatest number of different values.
The function counts the unique values in a 3x3 grid around each point and returns a matrix of indices (similar to \code{\link{which}(..., arr.ind = TRUE)}) for the maximum count (ties result in more than one pair of indices).
It can be used with the output from \code{diagram} for calculations in 2 dimensions to approximately locate the triple points on the diagram.

}

\value{
  For speciation diagrams, an \code{\link{invisible}} list of the chemical activities of the species, or their degrees of formation (if \code{alpha} is \code{TRUE}), at each point. For predominance diagrams, an invisible list with elements \code{species}, the dataframe describing the species, \code{out}, which species predominates at each grid point, and \code{A}, a list of the calculated values of the chemical affinity (per balanced quantity) (log10 dimensionless) at each point.
}

\section{Warning}{
The maximum affinity method is premised on the calculation of the affinities of formation reactions at equal activities of the species of interest.
Then, the species with the highest affinity of formation, after normalization by the balancing coefficients, corresponds to the predominant species in an equilibrium calculation.
The examples below \dQuote{work} because they are relatively simple - the balancing coefficients are unity or all the same value (aqueous aluminum example), or the species are solids with unit activities (the mineral examples).
The examples shown for \code{\link{protein}s} elsewhere also take the balancing coefficients to unity, after normalizing by protein length.
However, if aqueous species are present with different balancing coefficients, the maximum affinity method is not dependable, as shown in the TCA metabolite example below.
}

\seealso{ 
Other examples are present in the help for \code{\link{protein}} and \code{\link{buffer}}, and even more can be found in \code{\link{demos}}.
}

\examples{
\dontshow{data(thermo)}
### 1-D diagrams: logarithms of activities

## Aqueous sulfur species (after Seewald, 1997 and 2001)
basis("CHNOS+")
basis("pH", 5)
species(c("H2S", "S2-2", "S3-2", "S2O3-2", "S2O4-2", 
  "S3O6-2", "S5O6-2", "S2O6-2", "HSO3-", "SO2", "HSO4-"))
a <- affinity(O2=c(-50, -15), T=325, P=350)
e <- equilibrate(a, loga.balance=-2)
diagram(e, ylim=c(-30, 0), legend.x="topleft")
title(main=paste("Aqueous sulfur speciation, 325 degC, 350 bar\n",
  "After Seewald, 1997"))
# try it with and without the loga.balance argument (total activity of 
# the balanced quantity, in this case H2S aka sulfur) 

\donttest{## Degrees of formation of ionized forms of glycine
## After Fig. 1 of Aksu and Doyle, 2001 
basis("CHNOS+")
species(ispecies <- info(c("glycinium", "glycine", "glycinate")))
a <- affinity(pH=c(0, 14))
e <- equilibrate(a)
diagram(e, alpha=TRUE, lwd=1)
title(main=paste("Degrees of formation of aqueous glycine species\n",
  "after Aksu and Doyle, 2001"))

## Degrees of formation of ATP species as a function of 
## temperature, after LaRowe and Helgeson, 2007, Fig. 10b
# to make a similar diagram, activity of Mg+2 here is set to
# 10^-4, which is different from LH07, who used 10^-3 total molality
basis(c("CO2", "NH3", "H2O", "H3PO4", "O2", "H+", "Mg+2"),
  c(999, 999, 999, 999, 999, -5, -4))
species(c("HATP-3", "H2ATP-2", "MgATP-2", "MgHATP-"))
a <- affinity(T=c(0, 120, 25))
e <- equilibrate(a)
diagram(e, alpha=TRUE)  
title(main=paste("Degrees of formation of ATP species,\n",
  "pH=5, log(aMg+2)=-3. After LaRowe and Helgeson, 2007"),
  cex.main=0.9)

### 2-D diagrams: predominance diagrams
### these use the maximum affinity method

## Fe-S-O at 200 deg C, after Helgeson, 1970
basis(c("Fe", "O2", "S2"))
species(c("iron", "ferrous-oxide", "magnetite",
  "hematite", "pyrite", "pyrrhotite"))
# include the high-temperature phase of pyrrhotite
species("pyrrhotite", "cr2")
a <- affinity(S2=c(-50, 0), O2=c(-90, -10), T=200)
diagram(a, fill="heat")
title(main=paste("Fe-S-O, 200 degrees C, 1 bar",
  "After Helgeson, 1970", sep="\n"))

## pe-pH diagram for hydrated iron sulfides,
## goethite and pyrite, after Majzlan et al., 2006
# add some of these species to the database
add.obigt()
basis(c("Fe+2", "SO4-2", "H2O", "H+", "e-"), 
  c(0, log10(3), log10(0.75), 999, 999))
species(c("rhomboclase", "ferricopiapite", "hydronium jarosite",
  "goethite", "melanterite", "pyrite"))
a <- affinity(pH=c(-1, 4, 256), pe=c(-5, 23, 256))
d <- diagram(a, main="Fe-S-O-H, after Majzlan et al., 2006")
# the first four species show up along the top of the diagram
stopifnot(all.equal(unique(t(d$predominant)[256,]), 1:4))
water.lines(yaxis="pe")
text(3, 22, describe.basis(thermo$basis[2:3,], digits=2, oneline=TRUE))
text(3, 21, describe.property(c("T", "P"), c(25, 1), oneline=TRUE))
# reset the database
data(thermo)

## Aqueous Aluminum Species F-/OH-, after Tagirov and Schott, 2001
# some of the species are not in the default databse
add.obigt()
# the 999s have no effect on the diagram:
# pH and log_a(F-) are plotting variables
# O2 is not in the formation reactions
# Al+3 is the balanced quantity
basis(c("Al+3", "F-", "H+", "O2", "H2O"), c(rep(999, 4), 0))
species(c("Al+3", "Al(OH)4-", "AlOH+2", "Al(OH)2+", "Al(OH)3",
  "AlF+2", "AlF2+", "AlF3", "AlF4-", "Al(OH)2F2-", "Al(OH)2F",
  "AlOHF2"), "aq")
a <- affinity(pH=c(0, 10), "F-"=c(-1, -9), T=200)
diagram(a, fill="heat")
title(main=paste("Aqueous aluminium species, T=200 C, P=Psat\n",
  "after Tagirov and Schott, 2001"))
# restore thermodynamic database to default
data(thermo)

## Temperature-Pressure: kayanite-sillimanite-andalusite
# cf. Fig. 49 of Helgeson et al., 1978
# this is a system of one component (Al2SiO5), but we need the same 
# number of basis species as elements; and avoid using H2O or other 
# aqueous species because of T/P limits of the water() calculations;
basis(c("corundum", "quartz", "oxygen"))
species(c("kyanite", "sillimanite", "andalusite"))
# database has transition temperatures of kyanite and andalusite
# at 1 bar only, so we permit calculation at higher temperatures
a <- affinity(T=c(200, 900, 99), P=c(0, 9000, 101), exceed.Ttr=TRUE)
d <- diagram(a, fill=NULL)
bexpr <- sapply(c("Al2O3", "SiO2", "H2O"), expr.species, simplify=FALSE)
btext <- substitute(Al2O3 - SiO2 - H2O, unlist(bexpr))
mtitle(c(as.expression(btext), "after Helgeson et al., 1978"))
# find the approximate position of the triple point
tp <- find.tp(d$predominant)
Ttp <- a$vals[[1]][tp[1, 2]]
Ptp <- rev(a$vals[[2]])[tp[1, 1]]
points(Ttp, Ptp, pch=10, cex=5)
# some testing of the overall geometry
stopifnot(species()$name[d$predominant[1, 1]]=="andalusite")
stopifnot(species()$name[d$predominant[1, 101]]=="kyanite")
stopifnot(species()$name[d$predominant[99, 101]]=="sillimanite")

### other examples

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

## calculate the equilibrium logarithm of activity of a 
## basis species in different reactions
basis("CHNOS")
species(c("ethanol", "lactic acid", "deoxyribose", "ribose"))
a <- affinity(T=c(0, 150))
diagram(a, what="O2", legend.x="topleft", col=rev(rainbow(4)), lwd=2)
title(main="equilibrium logfO2 for 1e-3 mol/kg of CO2 and ... ")

## using strip()
## proteins from different mammals
organisms <- c("BOVIN", "CANFA", "HUMAN", "MOUSE", "PIG")
proteins <- c("AQP1", "MYG", "P53")
basis("CHNOS+")
species(rep(proteins, each=5), organisms)
a <- affinity(O2=c(-85, -65, 128))
ispecies <- list(1:5, 6:10, 11:15)
desc <- c("(Aquaporin-1)", "(Myoglobin)",
  "(Cellular tumor antigen p53)")
names(ispecies) <- paste(proteins, desc)
col <- rainbow(5)
strip(a, ispecies=ispecies, ymin=-0.7, col=col)
legend("bottomright", legend=organisms, col=col,
  lty=1, lwd=4, bty="n")
title(main=paste("Equilibrium degrees of formation of",
  "proteins from different mammals", sep="\n"))
}
}

\references{

  Aksu, S. and Doyle, F. M. (2001) Electrochemistry of copper in aqueous glycine solutions. \emph{J. Electrochem. Soc.} \bold{148}, B51--B57. \url{http://dx.doi.org/10.1149/1.1344532}

  Helgeson, H. C. (1970) A chemical and thermodynamic model of ore deposition in hydrothermal systems. \emph{Mineral. Soc. Amer. Spec. Pap.} \bold{3}, 155--186. \url{http://www.worldcat.org/oclc/583263}

  Helgeson, H. C., Delany, J. M., Nesbitt, H. W. and Bird, D. K. (1978) Summary and critique of the thermodynamic properties of rock-forming minerals. \emph{Am. J. Sci.} \bold{278-A}, 1--229. \url{http://www.worldcat.org/oclc/13594862}

  LaRowe, D. E. and Helgeson, H. C. (2007) Quantifying the energetics of metabolic reactions in diverse biogeochemical systems: electron flow and ATP synthesis. \emph{Geobiology} \bold{5}, 153--168. \url{http://dx.doi.org/10.1111/j.1472-4669.2007.00099.x}

  Majzlan, J., Navrotsky, A., McClesky, R. B. and Alpers, C. N. (2006) Thermodynamic properties and crystal structure refinement of ferricopiapite, coquimbite, rhomboclase, and Fe2(SO4)3(H2O)5. \emph{Eur. J. Mineral.} \bold{18}, 175--186. \url{http://dx.doi.org/10.1127/0935-1221/2006/0018-0175}

  Seewald, J. S. (1997) Mineral redox buffers and the stability of organic compounds under hydrothermal conditions. \emph{Mat. Res. Soc. Symp. Proc.} \bold{432}, 317--331. \url{http://lucy.mrs.org/meetings/spring96/Program/S.S96.html}

  Seewald, J. S. (2001) Aqueous geochemistry of low molecular weight hydrocarbons at elevated temperatures and pressures: Constraints from mineral buffered laboratory experiments \emph{Geochim. Cosmochim. Acta} \bold{65}, 1641--1664. \url{http://dx.doi.org/10.1016/S0016-7037(01)00544-0}

  Tagirov, B. and Schott, J. (2001) Aluminum speciation in crustal fluids revisited. \emph{Geochim. Cosmochim. Acta} \bold{65}, 3965--3992. \url{http://dx.doi.org/10.1016/S0016-7037(01)00705-0}

}

\keyword{primary}
