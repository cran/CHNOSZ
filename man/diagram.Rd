\name{diagram}
\alias{diagram}
\alias{strip}
\title{Calculate and Plot Equilibrium Chemical Activities of Species}
\description{
Plot chemical activity (speciation) or equal-activity (predominance) diagrams as a function of chemical activities of basis speecies, temperature and/or pressure. Or, plot values of chemical affinity, logK, logQ, Gibbs energies of basis species, species, and formation reactions, as a function of zero or one variables.
}

\usage{
  diagram(affinity, what = "logact", ispecies = NULL, balance = NULL, 
    names = NA, color = NA, add = FALSE, dotted = 0, cex = par("cex"),
    col = par("col"), pe = TRUE, pH = TRUE, ylog = TRUE, 
    main = NULL, cex.names = 1, legend.x = "topright", 
    lty = NULL, col.names = par("fg"), cex.axis=par("cex"),
    logact = NA, lwd = par("lwd"), alpha = FALSE,
    mar = NULL, residue = FALSE, yline = par("mgp")[1] + 1, 
    xrange = NULL, ylab = NULL, xlab = NULL, do.plot = TRUE,
    as.residue = FALSE, mam = TRUE, group = NULL, bg = par("bg"),
    side = 1:4, xlim = NULL, ylim = NULL)
  strip(affinity, ispecies = NULL, col = NULL, ns = NULL, 
    xticks = NULL, ymin = -0.2, xpad = 1)
}

\arguments{
  \item{affinity}{list, object returned by \code{\link{affinity}}.}
  \item{ispecies}{numeric, which species to consider (default of \code{NULL} is to consider all species).}
  \item{balance}{character, name of the balanced quantity in formation reactions, or numeric, coefficients to balance formation reactions.}
  \item{names}{character, names of species for activity lines or predominance fields.}
  \item{color}{character, colors of predominance fields, or colors of activity lines.}
  \item{add}{logical, add to current plot?}
  \item{dotted}{numeric, how often to skip plotting points on predominance field boundaries (to gain the effect of dotted or dashed boundary lines).}
  \item{cex}{numeric, character expansion (scaling relative to current).}
  \item{col}{character, color of predominance field boundaries and labels (\code{diagram}), or colors of bars in a strip diagram (\code{strip}).}
  \item{pe}{logical, convert an \samp{e-} axis to a \samp{pe} one? Default is \code{TRUE}; set this to \code{FALSE} to prevent this conversion.}
  \item{pH}{logical, convert an \samp{H+} axis to a \samp{pH} one?}
  \item{ylog}{logical, use a logarithmic y-axis (on 1D degree diagrams).}
  \item{main}{character, a main \code{\link{title}} for the plot. \code{NULL} means to plot no title.}
  \item{cex.names}{numeric, character expansion factor to be used for names of species on plots.}
  \item{legend.x}{character, description of legend placement passed to \code{\link{legend}}.}
  \item{lty}{numeric, line types to be used in plots.}
  \item{col.names}{character, colors for labels of species.}
  \item{cex.axis}{numeric, character expansion factor for names of axies.}
  \item{logact}{numeric, logarithm of total activity of balanced quantity (for speciation diagrams). If \code{NA}, the total activity of the balanced quantity is determined by the activities of the species.}
  \item{what}{character, what property to calculate and plot.}
  \item{lwd}{numeric, line width.}
  \item{alpha}{logical, for speciation diagrams, plot degree of formation instead of activities?}
  \item{mar}{numeric, margins of plot frame.}
  \item{residue}{logical, rewrite reactions to refer to formation of one mole of the balanced quantity (i.e., protein backbone group in proteins)?}
  \item{yline}{numeric, margin line on which to plot the y-axis name.}
  \item{xrange}{numeric, range of x-values between which predominance field boundaries are plotted.}
  \item{ylab}{character, label to use for y-axis.}
  \item{xlab}{character, label to use for x-axis.}
  \item{do.plot}{logical, make a plot?}
  \item{as.residue}{logical, make plot using activities of residues?}
  \item{mam}{logical, should maximum affinity method be used for 2-D diagrams?}
  \item{group}{list of numeric, groups of species to consider as a single effective species.}
  \item{bg}{character, background color for legend.}
  \item{side}{numeric, which sides of plot to draw axes.}
  \item{xlim}{numeric, limits of x-axis.}
  \item{ylim}{numeric, limits of y-axis.}
  \item{ns}{numeric, numbers of species, used to make inset plots for strip diagrams.}
  \item{xticks}{numeric, location of supplemental tick marks on x-axis.}
  \item{ymin}{numeric, lower limit of y-axis.}
  \item{xpad}{numeric, amount to extend x-axis on each side.}
}

\details{

  The \code{diagram} function takes as its primary input the results from \code{\link{affinity}} and displays diagrams representing either the thermodynamic properties of species, or the chemical activities of the species. The chemical activities are obtained by switching from an equal-activity reference state (values calculated using \code{\link{affinity}}) to an equal-affinity reference state (calculated using either a reaction matrix or the Maxwell-Boltzmann distribution). The activity diagrams that can be produced include chemical speciation diagrams as a function of one of temperature, pressure, or the chemical activities of the basis species, or equal-activity (predominance) diagrams as a function of two of these variables.

  To generate the stability relations from affinities of formation reactions, a reaction conservation rule is either automatically determined or specified by the user. For example, 3Fe2O3 = 2Fe3O4 + 1/2O2 is balanced on Fe (or a basis species containing Fe) and 4Fe2O3 + Fe = 3Fe3O4 is balanced on O (or a basis species containing O). The default action, if \code{balance} is \code{NULL}, is to balance on a basis species that appears in the formation reactions of all of the species of interest. The first such basis species (in order of their appearance in \code{thermo$basis}) will be used; this basis species is determined using \code{\link{which.balance}}. If a common basis species is not available, or if \code{balance} is \eqn{1}, the balance is set to unity. For proteins, an additional conservation rule is available; if \code{balance} is \code{PBB}, or if it is missing and all the species appear to be proteins (their names contain underscores) the metastability calculations will be balanced on protein length (number of protein backbone groups). 

  Predominance (2-D) diagrams are usually produced using the maximum affinity method, which is based on the notion that the predominant species at any point in the diagram is that one that has the greatest affinity of formation divided by the balanced quantity. This behavior can be altered by specifiying \code{mam} as \code{FALSE}, in which case the relative abundances of all species are calculated in the manner described below and the most abundant one at each grid point identified as the predominant species. However, this procedure can be very slow unless the reactions are cast in terms of residues (i.e., to use \code{\link{logact.mb}} instead of \code{\link{logact.react}}).

  When the coefficients of the balanced quantity are all equal, the location of the predominance field boundaries does not depend on the actual value of the activities of the species of interest, as long as they are all equal (these are "equal-activity" lines). Generally in the case in reactions among proteins, the coefficients of the balanced quantity are unequal in the different species, so the location of the predominance field boundaries does depend on the actual value chosen to represent equal activities of the species of interest. In this case the predominance field boundaries are "activity equal to X" lines. This is true for \code{mam} equal to TRUE; if \code{mam} is FALSE, the predominance fields really denote the most abundant species in a system, and the boundaries can represent different "equal-activity" values across the diagram.

  \code{residue}, if \code{TRUE}, instructs the function to rewrite the formation reactions so that each refers to formation of one mole of the balanced quantity (e.g., a residue of a protein for reactions conserving the protein backbone); the balance coefficients are then all unity. This is the recommended option for calculating relative abundances of proteins, because the large numbers that would otherwise be used as balance coefficients often create a situation of utter predominance of a single protein (see Dick, 2008). If \code{as.residue} is also \code{TRUE}, the calculated logarithms of activities of the residues are used in the plot and returned by the function, otherwise those of the species are shown. These options, however, are not restricted to systems of proteins. 

  If \code{group} is supplied, the activities of the species identified in each numeric vector of this list are summed together and subsequently treated as a single species. The names of the new species are taken from the names of this list. For 2-D diagrams, \code{group} only makes sense if \code{mam} is FALSE. An example of using the \code{group} argument is supplied in \code{\link{get.protein}}.

  For 1-D diagrams, the logarithmicity and limits of the y-axis can be controlled; the default (\code{ylog} is \code{TRUE}) is to plot logarithms of activities of species. If \code{alpha} is \code{TRUE}, the degrees of formation (ratios of activities of species to total activity) are plotted instead. The range of the y-axis on these diagrams can be controlled with \code{ylim} (which defaults to \code{c(0,1)} when \code{alpha} is \code{TRUE}). Also, a \code{\link{legend}} is placed at the location identified by \code{legend.x}, or omitted if \code{legend.x} is \code{FALSE}. 

  \code{what} indicates the type of calculation to perform. Usually, this means the metastable equilibrium logarithms of activities of the species of interest. Alternatively, \code{what} can be one of the (formulas of the) basis species, meaning to calculate the equilibrium logarithm of activity of that basis species in each of the formation reactions. Another option is to set \code{what} to \samp{A}, meaning to plot the affinities themselves (don't perform any equilibrium calculation). All of these actions assume that the property supplied in \code{a} is \samp{A} (i.e., affinities of formtation reactions). If it is anything else, \code{what} is ignored and the values themselves (as listed in \code{a}) are plotted.

  The \code{diagram} function by default starts a new plot, but if \code{add} is \code{TRUE} it adds to the current plot. If \code{do.plot} is FALSE, no plot will be generated but all the required computations will be performed and the results returned. The \code{name}s of the species, and the \code{color}s used to identify them (on chemical activity lines or predominance fields) can be changed; by default \code{\link{heat.colors}} are used to shade the predominance fields in on-screen diagrams. The \code{col} option controls the color of the predominance field boundaries, and \code{col.names} the color of the labels for the predominance fields. The line type (only on 1-D diagrams) and line width can be controlled with (\code{lty} and \code{lwd}, respectively). The line style on 2-D diagrams can be altered by supplying one or more non-zero integers in \code{dotted}, which indicates the fraction of line segments to omit. Finally, \code{cex}, \code{cex.names}, and \code{cex.axis} adjust the overall character expansion factors (see \code{\link{par}}) and those of the species names and axis labels.

  The x-axis label (1-D and 2-D diagrams) and y-axis label (2-D diagram) are automatically generated unless they are supplied in \code{xlab} and \code{ylab}. 
  A different incarnation of 1-D speciation diagrams is provided by \code{strip}. This function generates any number of strip diagrams in a single plot. The diagrams are made up of colors bars whose heights represent the relative abundances of species; the color bars are arranged in order of abundance and the total height of the stack of colors bars is constant. If \code{ispecies} is a list, the number of strip diagrams is equal to the number of elements of the list, and the elements of this list are numeric vectors that identify the species to consider for each diagram. The strips are labeled with the \code{\link{names}} of \code{ispecies}. If \code{col} is NULL, the colors of the bars are generated using \code{\link{rainbow}}. Supplemental ticks can be added to the x-axis at the locations specified in \code{xtick}; they are larger than the standard ticks and have colors corresponding to those of the color bars. \code{ymin} can be decreased in order to add more space at the bottom of the plot, and \code{xpad} can be changed in order to increase or decrease the size of the x-axis relative to the width of the strips. An inset dot-and-line plot is created below each strip if \code{ns} is given. This argument has the same format as \code{ispecies}, and can be used e.g. to display the relative numbers of species for comparison with the stability calculations.

}

\value{
  For speciation diagrams, an \code{\link{invisible}} list of the chemical activities of the species, or their degrees of formation (if \code{alpha} is \code{TRUE}), at each point. For predominance diagrams, an invisible list with elements \code{species}, the dataframe describing the species, \code{out}, which species predominates at each grid point, and \code{A}, a list of the calculated values of the chemical affinity (per balanced quantity) (log10 dimensionless) at each point.
}

\seealso{ 
  \code{\link{affinity}} for the source of the input for this function; \code{\link{par}} for how graphical parameters are set, \code{\link{protein}} and \code{\link{buffer}} for examples other than those shown below.
}

\examples{
\dontshow{data(thermo)}
## chemical affinities of formation (equal activity ref. state)
## and activities (equal affinity ref. state) of amino acids
opar <- par(mfrow=c(2,2))
mycol <- rev(rainbow(5))
basis("CHNOS+")
basis("H2S",-25)
aa <- c("glutamic acid","methionine","isoleucine",
  "leucine","tyrosine")
species(aa)
# affinities of reactions per CO2 at fixed conditions
a <- affinity()
diagram(a,what="A",col=mycol)
# affinities of reactions as a function of oxygen fugacity
a <- affinity(O2=c(-80,-66))
d <- diagram(a,what="A",balance=1,ylim=c(-80,30),
  cex.names=0.7,col=mycol,lwd=2)
title(main="Amino acid formation")
## logarithms of activities of amino acids
# balanced on amino acid backbone, using logact.react
d1 <- diagram(a,balance=1,legend.x=NULL,col=mycol,lwd=2)
# balanced on amino acid backbone, using logact.mb
d2 <- diagram(a,balance=1,residue=TRUE,do.plot=FALSE)
# the last two give the same results
stopifnot(all.equal(unlist(d2$logact),
  unlist(d1$logact)))
# in two dimensions
a <- affinity(H2O=c(-10,5),O2=c(-80,-66))
diagram(a,color=mycol)
# reset the plot device
par(opar)

### simple example of calculating the equilibrium logarithm of
### activity of a basis species
basis("CHNOS")
species(c("ethanol","acetic acid","lactic acid"))
a <- affinity(T=c(0,150))
diagram(a,what="O2")
# can also be done in 2 dimensions
a <- affinity(T=c(0,150),CO2=c(-10,0))
diagram(a,what="O2",ispecies=1,main=NULL)
diagram(a,what="O2",ispecies=2,col="red",add=TRUE)
diagram(a,what="O2",ispecies=3,col="blue",add=TRUE)
title(main="logfO2 ethanol/acetic/lactic (black/red/blue)")

### 1-D logarithm of activity plots
### (speciation diagrams)

## Aqueous sulfur species (after Seewald, 1997)
basis("CHNOS+")
basis("pH",5)
species(c("H2S","S2-2","S3-2","S2O3-2","S2O4-2","S3O6-2",
  "S5O6-2","S2O6-2","HSO3-","SO2","HSO4-"))
diagram(affinity(O2=c(-50,-15),T=325,P=350),
  ylim=c(-30,0),logact=-2,legend.x="topleft")
title(main=paste("Aqueous sulfur speciation, 325 degC, 350 bar\n",
  "After Seewald, 1997"))
# try it with and without the logact argument (total activity of 
# the balanced quantity, in this case H2S aka sulfur) 

\donttest{
## Degrees of formation of ionized forms of glycine
## as a function of pH
## After Fig. 1 of Aksu and Doyle, 2001 
basis("CHNOS+")
species(ispecies <- info(c("glycinium","glycine","glycinate")))
t <- affinity(pH=c(0,14))
diagram(t,alpha=TRUE,lwd=1)
title(main=paste("Degrees of formation of aqueous glycine species\n",
  "after Aksu and Doyle, 2001"))

## Degrees of formation of ATP species as a function of 
## temperature, after LaRowe and Helgeson, 2007
# Mg+2 here is the immobile component
# cf. LH07, where bulk composition of Mg+2 is specified
basis(c("CO2","NH3","H2O","H3PO4","O2","H+","Mg+2"),
  c(999,999,999,999,999,-5,-4))
species(c("HATP-3","H2ATP-2","MgATP-2","MgHATP-"))
diagram(affinity(T=c(0,120,25)),alpha=TRUE)  
title(main=paste("Degrees of formation of ATP species,\n",
  "pH=5, log(aMg+2)=-3. After LaRowe and Helgeson, 2007"),
  cex.main=0.9)

## using strip(): equilibrium abundances of 
## proteins from different mammals
organisms <- c("BOVIN","CANFA","HUMAN","MOUSE","PIG")
proteins <- c("AQP1","MYG","P53")
basis("CHNOS+")
species(rep(proteins,each=5),organisms)
a <- affinity(O2=c(-85,-65,128))
ispecies <- list(1:5,6:10,11:15)
desc <- c("(Aquaporin-1)","(Myoglobin)",
  "(Cellular tumor antigen p53)")
names(ispecies) <- paste(proteins,desc)
col <- rainbow(5)
strip(a,ispecies=ispecies,ymin=-0.7,col=col)
legend("bottomright",legend=organisms,col=col,
  lty=1,lwd=4,bty="n")
title(main=paste("Equilibrium abundances of proteins",
  "from different mammals",sep="\n"))

### predominance diagrams (equal activities of
### species as a function of two variables) 

## T-P stability diagram for SiO2
## after Ernst, 1976, Fig. 4.4
# the system is SiO2 (one component) but
# the basis species require all the elements:
# note that basis(c("SiO2","O2")) would identify SiO2(aq) 
# which is okay but brings in calculations of properties 
# of  water which are relatively slow.
basis(c("quartz","O2"))  # cr, gas
# browse the species of interest
info("SiO2 ")
# we'll take the crystalline ones
t <- info("SiO2","cr")
species(t)
# calculate chemical affinities
# the do.phases argument is necessary for 
# dealing with the phase transitions of minerals
t <- affinity(T=c(0,2000),P=c(0,30000),do.phases=TRUE)
diagram(t)
title(main="Phases of SiO2\nafter Ernst, 1976")

## Distribution of copper on Eh-pH diagram
## after Fig. 15 of Pourbaix, 1949
basis(c("Cu+2","Cl-","H2O","H+","e-"))
basis("Cl-",-2)
# two aqueous species, three solid ones
# (our CuCl is aq but it was cr in P's Fig. 15)
species(c("CuCl","Cu+2","copper","cuprite","tenorite"))
# (which is equivalent to ...)
# species(c("CuCl","Cu+2","Cu","Cu2O","CuO"),c("","","","","","cr"))
t <- affinity(pH=c(0,7),Eh=c(-0.1,1))
# this one corresponds to activity contours of 
# aqueous species at 10^-3 (the default aq activity in CHNOSZ)
diagram(t,color=NULL)
# here we set activities to unity; aq-cr boundaries change
species(c("CuCl","Cu+2"),c(0,0))
t <- affinity(pH=c(0,7),Eh=c(-0.1,1))
diagram(t,add=TRUE,names=NULL,col="blue",color=NULL)
water.lines()
title(main=paste("H2O-Cl-(Cu); activities of 10^-3 (black)\n",
  "or 0 (blue); after Pourbaix, 1949"))

## a pe-pH diagram for hydrated iron sulfides,
## goethite and pyrite, After Majzlan et al., 2006
# add some of these species to the database
danger()
basis(c("Fe+2","SO4-2","H2O","H+","e-"),c(0,log10(3),log10(0.75),999,999))
species(c("rhomboclase","ferricopiapite","hydronium jarosite",
  "goethite","melanterite","pyrite"))
a <- affinity(pH=c(-1,4),pe=c(-5,23))
diagram(a)
water.lines(yaxis="pe")
title(main=paste("Fe-S-O-H After Majzlan et al., 2006",
  describe(thermo$basis[2:3,],digits=2),sep="\n"))
# reset the database
data(thermo)

## cysteine-cysteinate-cystine Eh-pH
## at 25 and 100 deg C
basis("CHNOSe")
species(c("cysteine","cysteinate","cystine"))
a <- affinity(pH=c(5,10),Eh=c(-0.5,0))
diagram(a,color=NULL)
water.lines()
a <- affinity(pH=c(5,10),Eh=c(-0.5,0),T=100)
diagram(a,col="red",add=TRUE,names=NULL)
water.lines(T=convert(100,"K"),col="red")
title(main=paste("Cysteine Cysteinate Cystine",
  "25 and 100 deg C",sep="\n"))

## Soil Organic Matter CO2-H2O, O2-H2O (after Tardy et al., 1997)
# NH3 is conserved, and H2O is on an axis of the diagram
# formulas for aqueous species, names for phases ...
danger()
basis(c("NH3","water","CO2","O2"),c(999,999,-2.5,-28))
# switch to gaseous CO2 (aq is the default)
basis("CO2","gas")
# load the species of interest
species(c("microflore","plantes","acide fulvique",
  "acide humique","humine"))
# proceed with the diagrams
diagram(affinity(H2O=c(-0.6,0.1),CO2=c(-3,-1)))
title(main=paste("Relative stabilities of soil organic matter\n",
  "after Tardy et al., 1997"))
# this is the O2-H2O diagram
# diagram(affinity(H2O=c(-1,0.5),O2=c(-28.5,-27.5)))
data(thermo)

## Aqueous Aluminum Species F-/OH- (afer Tagirov and Schott, 2001)
# in order to reproduce this calculation, we have to 
# consider the properties of species given by these authors,
# which are not the default ones in the database
danger()
# The coefficients on H+ and O2 in all the formation reactions
# are zero, so the number of basis species here is three. Al+3 
# becomes the conservant, and F- and OH- are being plotted ...
# so their initial activities don't have to be set.
basis(c("Al+3","F-","OH-","H+","O2"),rep(999,5))
species(c("Al+3","Al(OH)4-","AlOH+2","Al(OH)2+","Al(OH)3",
  "AlF+2","AlF2+","AlF3","AlF4-","Al(OH)2F2-","Al(OH)2F","AlOHF2"))
# Increase the x- and y- resolution from the default and calculate
# and plot predominance limits. Names of charged basis species,
# such as "H+", "e-" and the ones shown here, should be quoted
# when given as arguments to affinity(). The OH- values shown here
# correspond to pH=c(0,14) (at unit activity of water).
a <- affinity("OH-"=c(-14,0),"F-"=c(-1,-8),T=200)
diagram(a)
title(main=paste("Aqueous aluminium species, T=200 C, P=Psat\n",
  "after Tagirov and Schott, 2001"))
# We could do this to overlay boundaries for a different pressure
#a.P <- affinity("OH-"=c(-14,0),"F-"=c(-1,-8),T=200,P=5000)
#diagram(a.P,names=NULL,color=NULL,col="blue",add=TRUE)
# restore thermodynamic database to default
data(thermo)

## Fe-S-O at 200 deg C, After Helgeson, 1970
basis(c("Fe","O2","S2"))
species(c("iron","ferrous-oxide","magnetite",
  "hematite","pyrite","pyrrhotite"))
a <- affinity(S2=c(-50,0),O2=c(-90,-10),T=200)
diagram(a,color="heat")
title(main=paste("Fe-S-O, 200 degrees C, 1 bar",
  "After Helgeson, 1970",sep="\n"))

## Temperature-Pressure: kayanite-sillimanite-andalusite
# this is a system of a single component (Al2SiO5)
# but we still have to use the same number of 
# basis species as the number of elements
basis(c("kyanite","quartz","oxygen"))
species(c("kyanite","sillimanite","andalusite"))
diagram(affinity(T=c(0,1000,200),P=c(500,5000,200)),color=NULL)
title(main="Al2SiO5")
} # end donttest
}

\references{

  Bowers, T. S., Jackson, K. J. and Helgeson, H. C., 1984. \emph{Equilibrium Activity Diagrams for Coexisting Minerals and Aqueous Solutions at Pressures and Temperatures to 5 kb and 600} \eqn{^{\circ}}{degrees }\emph{C}. Springer-Verlag, Heidelberg, 397 p. \url{http://www.worldcat.org/oclc/224591948}

  Dick, J. M., 2008. Calculation of the relative metastabilities of proteins using the CHNOSZ software package. \emph{Geochem. Trans.}, 9:10. \url{http://dx.doi.org/10.1186/1467-4866-9-10}

  Ernst, W. G., 1976. \emph{Petrologic Phase Equilibria}. W. H. Freeman, San Francisco, 333 p. \url{http://www.worldcat.org/oclc/2072721}
 
  Helgeson, H. C., 1970. A chemical and thermodynamic model of ore deposition in hydrothermal systems. \emph{Mineral. Soc. Amer. Spec. Pap.}, 3, 155-186. \url{http://www.worldcat.org/oclc/583263}

  Helgeson, H. C., Delany, J. M., Nesbitt, H. W. and Bird, D. K., 1978. Summary and critique of the thermodynamic properties of rock-forming minerals. \emph{Am. J. Sci.}, 278-A, 1-229. \url{http://www.worldcat.org/oclc/13594862}

  LaRowe, D. E. and Helgeson, H. C., 2007. Quantifying the energetics of metabolic reactions in diverse biogeochemical systems: electron flow and ATP synthesis. \emph{Geobiology}, 5, 153-168. \url{http://dx.doi.org/10.1111/j.1472-4669.2007.00099.x}

  Majzlan, J., Navrotsky, A., McClesky, R. B. and Alpers, C. N., 2006. Thermodynamic properties and crystal structure refinement of ferricopiapite, coquimbite, rhomboclase, and Fe2(SO4)3(H2O)5. \emph{Eur. J. Mineral.}, 18, 175-186. \url{http://dx.doi.org/10.1127/0935-1221/2006/0018-0175}

  Pourbaix, M. J. N., 1949. \emph{Thermodynamics of dilute aqueous solutions}. Edward Arnold & Co., London, 136 p. \url{http://www.worldcat.org/oclc/1356445}

  Seewald, J. S., 1997. Mineral redox buffers and the stability of organic compounds under hydrothermal conditions. \emph{Mat. Res. Soc. Symp. Proc.}, 432, 317-331. \url{http://lucy.mrs.org/meetings/spring96/Program/S.S96.html}

  Tagirov, B. and Schott, J., 2001. Aluminum speciation in crustal fluids revisited. \emph{Geochim. Cosmochim. Acta}, 65, 3965-3992. \url{http://dx.doi.org/10.1016/S0016-7037(01)00705-0}

  Tardy, Y., Schaul, R. and Duplay, J., 1997. Thermodynamic stability fields of humus, microflora and plants. \emph{C. R. Acad. Sci. Paris}, 324, 969-976. \url{http://dx.doi.org/10.1016/S1251-8050(97)83981-X}

}

\keyword{misc}
