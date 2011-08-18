\name{subcrt}
\alias{subcrt}
\alias{read.supcrt}
\alias{write.supcrt}
\title{Properties of Species and Reactions}
\description{
  Calculate the standard molal thermodynamic properties of one or more species or a reaction between species as a function of temperature and pressure. Import or export thermodynamic data in \acronym{SUPCRT} format.
}

\usage{
  subcrt(species, coeff = 1, state = NULL,
    property = c("logK","G","H","S","V","Cp"),
    T = seq(273.15,623.15,25), P = "Psat", grid = NULL, 
    convert = TRUE, do.phases = TRUE, logact = NULL,
    action.unbalanced = "warn", IS = 0)
  read.supcrt(file)
  write.supcrt(file = "supcrt.dat", obigt = thermo$obigt)
}

\arguments{
  \item{species}{character, name or formula of species, or numeric, rownumber of species in \code{thermo$obigt}.}
  \item{coeff}{numeric, reaction coefficients on species.}
  \item{state}{character, state(s) of species.}
  \item{property}{character, property(s) to calculate.}
  \item{T}{numeric, temperature(s) of the calculation.}
  \item{P}{numeric, pressure(s) of the calculation, or character, \samp{Psat}.}
  \item{grid}{character, type of \code{P}\eqn{\times}{x}\code{T} grid to produce (NULL, the default, means no gridding).}
  \item{do.phases}{logical, replace Gibbs energies of phases of minerals and of other species beyond their upper temperature limits of stability or extrapolation with 999999?}
  \item{logact}{numeric, logarithms of activities of species in reaction.}
  \item{file}{character, name of \acronym{SUPCRT92} data file.}
  \item{obigt}{dataframe, thermodynamic data.} 
  \item{convert}{logical, are input and output units of T and P those of the user (\code{TRUE}) (see \code{\link{nuts}}), or are they Kelvin and bar (\code{FALSE})?}
  \item{action.unbalanced}{character \samp{warn} or NULL, what action to take if unbalanced reaction is provided.}
  \item{IS}{numeric, ionic strength(s) at which to calculate apparent molal properties, mol kg\eqn{^{-1}}{^-1}.}
}

\details{

 
  \code{subcrt} calculates the standard molal thermodynamic properties of species and reactions as a function of temperature and pressure. For each of the \code{species} (a formula or name), optionally identified in a given \code{state}, the standard molal thermodynamic properties and equations-of-state parameters are retrieved via \code{\link{info}} (except for \eqn{\mathrm{H_2O}}{H2O} liquid). The standard molal properties of the species are computed using equations-of-state functions for aqueous species (\code{\link{hkf}}), crystalline, gas, and liquid species (\code{\link{cgl}}) and liquid or supercritical \eqn{\mathrm{H_2O}}{H2O} (\code{\link{water}}).

  \code{T} and \code{P} denote the temperature and pressure conditions for the calculations and should generally be of the same length, unless \code{P} is \samp{Psat} (the default; see \code{\link{water}}). Argument \code{\link{grid}} if present can be one of \code{T} or \code{P} to perform the computation of a \code{T}\eqn{\times}{x}\code{P} or \code{P}\eqn{\times}{x}\code{T} grid. The \code{property}s to be calculated can be one or more of those shown below:

 \tabular{lll}{
         \code{rho}     \tab Density of water                      \tab g cm\eqn{^{-3}}{^-3} \cr
         \code{logK}    \tab Logarithm of equilibrium constant     \tab dimensionless \cr
         \code{G}       \tab Gibbs energy                          \tab (cal | J) mol\eqn{^{-1}}{^-1} \cr
         \code{H}       \tab Enthalpy                              \tab (cal | J) mol\eqn{^{-1}}{^-1} \cr
         \code{S}       \tab Entropy                               \tab (cal | J) K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \cr
         \code{V}       \tab Volume                                \tab cm\eqn{^3} mol\eqn{^{-1}}{^-1} \cr
         \code{Cp}      \tab Heat capacity                         \tab (cal | J) K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} \cr
         \code{E}       \tab Exapansibility                        \tab cm\eqn{^3} K\eqn{^{-1}}{^-1} \cr
         \code{kT}      \tab Isothermal compressibility            \tab cm\eqn{^3} bar\eqn{^{-1}}{^-1} \cr
 }

  (Note that \code{E} and \code{kT} can only be calculated for aqueous species and only if the option (\code{\link{thermo}$opt$water}) for calculations of properties using \code{water} is set to \code{IAPWS}. On the other hand, if the \code{water} option is \samp{SUPCRT} (the default), \code{E} and \code{kT} can be calculated for water but not for aqueous species. (This is not an inherent limitation in either formulation, but it is just a matter of implementation.))

  Depending on the global units definition (see \code{\link{nuts}}) the values of \code{G}, \code{H}, \code{S} and \code{Cp} are returned in calories or Joule, but only if \code{convert} is \code{TRUE}. Likewise, the input values of \code{T} and \code{P} should be in units specified through \code{nuts}, but setting \code{convert} to \code{FALSE} forces \code{subcrt} to treat them as Kelvin and bar, respectively.

  A chemical reaction is defined if \code{coeff} is given. In this mode the standard molal properties of species are summed according to the stoichiometric \code{coeff}icients, where negative values denote reactants. Reactions that do not conserve elements are permitted; \code{subcrt} prints the missing composition needed to balance the reaction and produces a warning but computes anyway. Alternatively, if the \code{\link{basis}} species of a system were previously defined, and if the species being considered are within the compositional range of the basis species, an unbalanced reaction given in the arguments to \code{subcrt} will be balanced automatically, without altering the coefficients on the species identified in the arguments (unless perhaps they correspond to basis species), and without a warning. However, if a reaction is unbalanced and \code{action.unbalanced} is set to NULL, no warnings are generated and no attempt is made to balance the reaction.

  Minerals with phase transitions (denoted by having states \samp{cr1}, \samp{cr2} etc.) can be defined generically by \samp{cr} in the \code{state} argument. As of CHNOSZ-0.9-6, to consider phase transitions the \code{species} must be character, not numeric. \code{subcrt} in this case simultaneously calculates the requested properties of all the phases of each such mineral (phase species) and, using the values of the transition temperatures calculated from those at P = 1 bar given in the thermodynamic database together with functions of the entropies and volumes of transitions (see \code{\link{dPdTtr}}), determines the stable phase of the mineral at any grid point and substitutes the properties of this phase at that point for further calculations. If individual phase species of minerals are specified (by \samp{cr1}, \samp{cr2} etc. in \code{state}), the Gibbs energies of these minerals are assigned values of 999999 at conditions where another of the phases is stable (only if \code{do.phases} is \code{TRUE}). (NOTE: if you set \code{do.phases} to \code{FALSE} while investigating the properties of phases of minerals identified generically (by \samp{cr}) the function will identify the stable phase on the basis not of the transition temperatures but of the calculated Gibbs energies of the phase species. This is generally not advised, since the computed standard molal properties of a phase species lose their physical meaning if computed beyond the transition temperatures of the corresponding phase.)

  The chemical affinities of reactions are calculated if \code{logact} is provided. This argument indicates the logarithms of activities (fugacities for gases) of species in the reaction; if there are fewer values of \code{logact} than number of species those values are repeated as necessary. If the reaction was unbalanced to start, the logarithms of activities of any basis species added to the reaction are taken from the global basis definition. Columns appended to the output are \code{logQ} for the activity product of the reaction, and \code{A} for the chemical affinity (the latter in units corresponding to those of Gibbs energy). Note that \code{\link{affinity}} provides related functionality but is geared toward the properties of formation reactions of species from the basis species and can be performed in more dimensions. Calculations of chemical affinity in \code{subcrt} can be performed for any reaction of interest; however, they are currently limited to constant values of the logarithms of activities of species in the reactions, and hence of \code{logQ}, across the computational range.

  If \code{IS} is set to a single value other than zero, \code{\link{nonideal}} is used to calculate the apparent properties (\code{G}, \code{H}, \code{S} and \code{Cp}) of charged aqueous species at the given ionic strength. To perform calculations at a single \code{P} and \code{T} and for multiple values of ionic strength, supply these values in \code{IS}. Calculations can also be performed on a \code{P}-\code{IS}, \code{T}-\code{IS} or \code{P,T}-\code{IS} grid. Values of \code{logK} of reactions calculated for \code{IS} not equal to zero are consistent with the apparent Gibbs energies of the charged aqueous species.

  \code{subcrt} is modeled after the functionality of the \acronym{SUPCRT92} package (Johnson et al., 1992). Certain features of \acronym{SUPCRT92} are not available here, for example, calculations as a function of density of \eqn{\mathrm{H_2O}}{H2O} instead of pressure, or calculations of temperatures of univariant curves (i.e. for which \code{logK} is zero), or calculations of the molar volumes of quartz and coesite as a function of temperature (taken here to be constant). The informative messages produced by \code{SUPCRT92} when temperature or pressure limits of the equations of state are exceeded generally are not reproduced here. However, \code{NA}s may be produced in the output of \code{subcrt} if the requisite thermodynamic or electrostatic properties of water can not be calculated at given conditions. 

  \code{read.supcrt} and \code{write.supcrt} are used to read and write thermodynamic data \code{file}s in the format of \acronym{SUPCRT92}. \code{read.supcrt} parses the thermodynamic data contained in a \acronym{SUPCRT92}-format file into a dataframe that is compatible with the format used in CHNOSZ (see \code{\link{thermo}$obigt}). \code{write.supcrt} creates a \acronym{SUPCRT92} \code{file} using the dataframe given in the \code{obigt} argument (if missing the entire species database is processed). An edited version of the \samp{slop98.dat} data file (Shock et al., 1998) that is compatible with \code{read.supcrt} is available at \url{http://chnosz.net/download/} (small formatting quirks in the original version cause glitches with this function).

}

\section{Warning}{

  Comparison of the output of \code{subcrt} with that of \acronym{SUPCRT92} indicates the two give similar values of properties for neutral aqueous species up to 1000 \eqn{^{\circ}}{degrees }C and 5000 bar and for charged aqueous species over the temperature range 0 to 300 \eqn{^{\circ}}{degrees }C. At higher temperatures, there are significant divergences for charged species. For example, there is less than a 2 percent difference in the value of \eqn{{\Delta}G^{\circ}}{Delta G0} for K+(aq) at \eqn{^{\circ}}{degrees }C and 5000 bar, but this deviation increases with decreasing pressure at the same temperature. Further testing is necessary in CHNOSZ in connection with the \eqn{g}{g}- and \eqn{f}{f}-functions (Shock et al., 1992) for the partial derivatives of the omega parameter which are incorporated into the \code{\link{hkf}} function. (Note in particular the NaCl dissociation example under \code{\link{water}}, the results of which are noticeably different from the calculations of Shock et al., 1992.)
 
}

\value{
  For \code{subcrt}, a list of length two or three. If the properties of a reaction were calculated, the first element of the list (named \samp{reaction}) contains a dataframe with the reaction parameters; the second element, named \samp{out}, is a dataframe containing the calculated properties. Otherwise, the properties of species (not reactions) are returned: the first element, named \samp{species}, contains a dataframe with the species identification; the second element, named \samp{out}, is itself a list, each element of which is a dataframe of properties for a given species. If minerals with phase transitions are present, a third element (a dataframe) in the list indicates for all such minerals the stable phase at each grid point.

  \code{read.supcrt} returns a dataframe of properties and parameters of species with the same structure as \code{\link{thermo}$obigt}.
}

\seealso{
  \code{\link{info}} for querying the thermodynamic database; \code{\link{makeup}} for parsing formulas to check reaction balance; \code{\link{water}}, \code{\link{hkf}} and \code{\link{cgl}} for equations of state calculations.
}

\examples{
  \dontshow{data(thermo)}
  ## properties of species
  subcrt("water")
  # calculating at different temperatures
  subcrt("water",T=seq(0,100,10))
  # calculating at even increments  	
  subcrt("water",T=seq(500,1000,length.out=10),
    P=seq(5000,10000,length.out=10))
  # calculating on a temperature-pressure grid
  subcrt("water",T=c(500,1000),P=c(5000,10000),grid="P")
  # to calculate only selected properties
  subcrt("water",property=c("G","E"))	
  # the properties of multiple species
  subcrt(c("glucose","ethanol","CO2"))

  ## input/output in different units
  nuts(c("K","J"))
  subcrt("water")
  subcrt("water",T=298.15)
  nuts(c("cal","C"))   # back to defaults

  ## properties of reactions
  subcrt(c("H2O","H+","k-feldspar","kaolinite","K+","SiO2"),
    c(-1,-2,-2,1,2,4))
  subcrt(c("glucose","ethanol","CO2"),c(-1,2,2))
  # to specify the states
  subcrt(c("glucose","ethanol","CO2"),c(-1,2,2),c("aq","aq","gas")) 
  # this warns about an unbalanced reaction
  subcrt(c("glucose","ethanol"),c(-1,3))
  
  ## auto balancing reactions
  # the basis species must first be defined
  basis(c("CO2","H2O","NH3","H2S","O2"))
  subcrt(c("glucose","ethanol"),c(-1,3))
  # a bug in CHNOSZ <0.9 caused the following
  # to initiate an infinite loop
  basis(c("H2O","H2S","O2","H+"))
  subcrt(c("HS-","SO4-2"),c(-1,1))
  # note the next one is a non-reaction
  # (products same as reactants)
  subcrt("H2O",1)
  # but this one auto-balances into a reaction
  # (water to steam)
  subcrt("H2O",1,"gas")
  # properties of a species and a formation 
  # reaction for that species
  subcrt("C2H5OH")    # species
  basis("CHNOS")
  subcrt("C2H5OH",1)  # reaction
  
  ## properties of mineral phases
  # properties of phase species
  subcrt(c("pyrrhotite","pyrrhotite"),state=c("cr1","cr2"))
  # properties of the stable phases
  subcrt("pyrrhotite")
  # phase transitions in a reaction
  subcrt(c("pyrite","pyrrhotite","H2O","H2S","O2"),c(-1,1,-1,1,0.5))

  ## these produce NAs and warnings
  # Psat is not defined above the critical point
  subcrt("alanine",T=seq(0,5000,by=1000))
  # above the T, P limits for the H2O equations of state
  subcrt("alanine",T=c(2250,2251),P=c(30000,30001),grid="T")

  \donttest{
    ## heat capacity of Fe(cr)
    # compare calculated values of heat capacity with
    # values from Robie and Hemingway, 1995, (from which 
    # the parameters in the database were derived)
    nuts(c("J","K"))  # set the units
    # we set pressure here otherwise subcrt goes for 
    # Psat (saturation vapor pressure of H2O above 
    # 100 degrees C) which can not be calculated above 
    # the critical point of H2O (~647 K)
    s <- subcrt("Fe",T=seq(300,1800,20),P=1)
    plot(s$out[[1]]$T,s$out[[1]]$Cp,type="l",
      xlab=axis.label("T"),ylab=axis.label("Cp"))
    # add points from RH95
    RH95 <- read.csv(system.file("extdata/cpetc/RH95.csv",package="CHNOSZ"))
    points(RH95[,1],RH95[,2])
    title(main=paste("Heat capacity of Fe(cr)\n",
      "(points - Robie and Hemingway, 1995)"))
    # reset the units
    nuts(c("C","cal"))
    
    ## Skarn example after Johnson et al., 1992
    P <- seq(500,5000,500)
    # this is like the temperature specification used 
    # in the example by Johnson et al., 1992
    # T <- seq(0,1000,100)
    # we use this one to avoid warnings at 0 deg C, 5000 bar
    T <- c(2,seq(100,1000,100))
    subcrt(c("andradite","carbon dioxide","H2S","Cu+","quartz","calcite",
      "chalcopyrite","H+","H2O"),coeff=c(-1,-3,-4,-2,3,3,2,2,3),
      state=c("cr","g","aq","aq","cr","cr","cr","aq","liq"),
      P=P,T=T,grid="P")
    # the volumes differ from SUPCRT92, at least because CHNOSZ
    # doesn't have include volume changes for quartz
    # other differences are from updated properties of species
    # e.g. Cu+ from Shock et al., 1997

    ## Standard Gibbs energy of reactions with HCN and 
    ## formaldehyde, after Schulte and Shock, 1995 Fig. 1
    rxn1 <- subcrt(c("formaldehyde","HCN","H2O","glycolic acid","NH3"),
      c(-1,-1,-2,1,1),P=300)
    rxn2 <- subcrt(c("formaldehyde","HCN","H2O","glycine"),
      c(-1,-1,-1,1),P=300)
    plot(x=rxn1$out$T,rxn1$out$G/1000,type="l",ylim=c(-40,-10),
      xlab=axis.label("T"),ylab=axis.label("DG0r","k"))
    lines(rxn1$out$T,rxn2$out$G/1000)
    # write the reactions on the plot
    text(150,-14,describe(rxn1$reaction,
      use.name=c(TRUE,FALSE,FALSE,TRUE,FALSE)))
    text(200,-35,describe(rxn2$reaction,
      use.name=c(TRUE,FALSE,FALSE,TRUE)))
    title(main=paste("Standard Gibbs energy of reactions",
      "after Schulte and Shock, 1995",sep="\n"))

    ## Calculation of chemical affinities
    # after LaRowe and Helgeson, 2007
    # Fig. 3 (a): reduction of nicotinamide adenine 
    # dinucleotide (NAD) coupled to oxidation of glucose
    # list the available NAD species
    info("NAD ")
    T <- seq(0,120,10)
    # oxidation of glucose (C6H12O6)
    basis(c("glucose","H2O","NH3","CO2","H+"),c(-3,0,999,-3,-7))
    s <- subcrt(c("NAD(ox)-","NAD(red)-2"),c(-12,12),logact=c(0,0),T=T)
    # LH07's diagrams are shown per mole of electron (24 e- per 12 NAD)
    A <- s$out$A/24/1000
    plot(x=T,y=A,xlim=range(T),ylim=c(1.4,5.4),
      xlab=axis.label("T"),ylab=axis.label("A",opt="k"),type="l")
    text("NAD(ox)-/NAD(red)-2 = 1",x=median(T),y=median(A))
    # different activity ratio
    s <- subcrt(c("NAD(ox)-","NAD(red)-2"),c(-12,12),logact=c(1,0),T=T)
    A <- s$out$A/24/1000
    lines(x=T,y=A)
    text("NAD(ox)-/NAD(red)-2 = 10",x=median(T),y=median(A))
    # different activity ratio
    s <- subcrt(c("NAD(ox)-","NAD(red)-2"),c(-12,12),logact=c(0,1),T=T)
    A <- s$out$A/24/1000
    lines(x=T,y=s$out$A/24/1000)
    text("NAD(ox)-/NAD(red)-2 = 0.1",x=median(T),y=median(A))
    # this command prints the reaction on the plot
    text(40,4.5,c2s(s2c(describe(s$reaction,
      use.name=c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)),
      sep="="),sep="\n"))
    # label the plot
    title(main=paste("Chemical affinity of NAD reduction",
     "after LaRowe and Helgeson, 2007",sep="\n"),
     sub=describe(thermo$basis,T=NULL))

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
    lb <- length(basis())
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
    iplotfun(25,"black")
    iplotfun(40,"red")
    title(main=paste("activity coefficient (gamma) of -1,-2,-3,-4",
      "charged species at 0, 25, 40 deg C, after Alberty, 2003",
      sep="\n"),cex.main=0.95)
  }
}

\references{
  Alberty, R. A. (2003) \emph{Thermodynamics of Biochemical Reactions}, John Wiley & Sons, Hoboken, New Jersey, 397 p. \url{http://www.worldcat.org/oclc/51242181}
  
  Johnson, J. W., Oelkers, E. H. and Helgeson, H. C. (1992) SUPCRT92: A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000\eqn{^{\circ}}{degrees }C. \emph{Comp. Geosci.} \bold{18}, 899--947. \url{http://dx.doi.org/10.1016/0098-3004(92)90029-Q}

  LaRowe, D. E. and Helgeson, H. C. (2007) Quantifying the energetics of metabolic reactions in diverse biogeochemical systems: electron flow and ATP synthesis. \emph{Geobiology} \bold{5}, 153--168. \url{http://dx.doi.org/10.1111/j.1472-4669.2007.00099.x}

  Schulte, M. D. and Shock, E. L. (1995) Thermodynamics of Strecker synthesis in hydrothermal systems. \emph{Orig. Life Evol. Biosph.} \bold{25}, 161--173. \url{http://dx.doi.org/10.1007/BF01581580}

  Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \eqn{^{\circ}}{degrees }C and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \url{http://dx.doi.org/10.1039/FT9928800803}

  Shock, E. L. et al. (1998) \emph{SLOP98.dat} (computer data file). \url{http://geopig.asu.edu/supcrt92_data/slop98.dat}, accessed on 2005-11-05.

}

\keyword{primary}
