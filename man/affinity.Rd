\name{affinity}
\alias{affinity}
\alias{energy}
\alias{energy.args}
\title{Chemical Affinities of Formation Reactions}
\description{
  Calculate the chemical affinities of formation reactions of species. Do it for single values of temperature, pressure, ionic strength and chemical activities of the basis species, or as a function of one or two of these variables. Or, return other properties including standard molal Gibbs energies of basis species and species of interest, and standard molal Gibbs energies, equilibrium constants and reaction activity products of formation reactions. 
}

\usage{
  energy(what, vars, vals, lims, T=thermo$opt$Tr, P="Psat", IS=0, 
    sout=NULL, do.phases=FALSE)
  energy.args(args, quiet = FALSE)
  affinity(..., property=NULL, sout=NULL, do.phases=FALSE,
    return.buffer=FALSE, balance="PBB", quiet=FALSE, 
    iprotein=NULL, logact.protein=-3)
}

\arguments{
    \item{what}{character, name of property to calculate.}

    \item{vars}{character, names of variables over which to calculate a property.}

    \item{vals}{list of numeric, grid values for each variable.}

    \item{lims}{list of numeric, limits of the values for each variable.}

    \item{T}{numeric, temperature. Default is to take the temperature from \code{\link{thermo}$opt$Tr}, which corresponds to 25 degrees C.}

    \item{P}{numeric, pressure, or character "Psat" (default), which denotes 1 bar or the saturation vapor pressure of water above 100 degrees C (see \code{\link{water}}).}

    \item{IS}{numeric, ionic strength default is 0 mol kg\eqn{^{-1}}{^-1}.}

    \item{sout}{list, output from subcrt function.}

    \item{do.phases}{logical, allow subcrt to compute properties for stable phases?}

    \item{args}{list, defines the variables over which to calculate properties.}

    \item{quiet}{logical, print fewer messages to the screen?}

    \item{...}{numeric, zero or more named arguments, used to identify the variables of interest in the calculations.}

    \item{property}{character, denoting the property to be calculated. Default is \samp{A}, for chemical affinity of formation reactions of species of interest.}

    \item{return.buffer}{logical. If \code{TRUE}, and a \code{\link{buffer}} has been associated with one or more basis species in the system, return the values of the activities of the basis species calculated using the buffer (it is not necessary in this case to have defined any species of interest). Default is \code{FALSE}.}

    \item{balance}{character. This argument is passed to \code{\link{buffer}} to identify a conserved basis species or other quantity in a chemical activity buffer. Default is \samp{PBB}.}

    \item{iprotein}{numeric, indices of proteins in \code{\link{thermo}$protein} for which to calculate properties.}

    \item{logact.protein}{numeric, logarithms of activities of proteins identified in \code{protein}.}
  
}

\details{

  The chemical affinity relates the driving force for a reaction to proceed in a forward or reverse direction to the standard Gibbs energy of the reaction and to the physical and chemical environment. The calculation of chemical affinity \emph{\bold{A}} is expressed by \emph{\bold{A}}=\eqn{RT\ln (K/Q)}{RT*ln(K/Q)} (Kondepudi and Prigogine, 1998), where \eqn{K} denotes the equilibrium constant of the reaction and \eqn{Q} stands for the activity product of the species in the reaction. (The equilibrium constant is related to standard Gibbs energy of reaction, \eqn{\Delta G^{\circ}_r}{DeltaG0r}, by \eqn{\Delta G^{\circ}_r = -2.303RT\log K}{DeltaG0r = -2.303*RT*logK}, where \eqn{R} and \eqn{T} stand for, respectively, the gas constant and temperature). With the approach of a given reaction to a state of equilibrium, the chemical affinity tends toward zero, or \eqn{K = Q}. 

  In normal practice, the user calls \code{affinity} in order to calculate the chemical affinities of formation reactions of species of interest. This function in turn calls \code{energy.args} and \code{energy} to perform the calculations. These functions rely on the global definitions of the \code{\link{basis}} species and \code{\link{species}} of interest. 

  \code{energy} is the workhorse of the \code{affinity} calculations. Given \eqn{n} (which can be zero, one, or more) names of basis species and/or \samp{T}, \samp{P}, or \samp{IS} as the \code{var}s, it calculates the property given in \code{what} on an \eqn{n}-dimensional grid for each of the values (\code{vals}) of the corresponding variable. The limits for each variable given in \code{lims} indicate the minimum and maximum value and, if a third value is supplied, the resolution, or number of points in the given dimension. If \samp{T}, \samp{P}, and/or \samp{IS} are not among the \code{var}s, their constant values can be supplied in \code{T} (in Kelvin), \code{P} (in bar, or \samp{Psat}), and \code{IS} (in mol kg\eqn{^{-1}}{^-1}). \code{sout}, if provided, replaces the call to \code{\link{subcrt}} which can greatly speed up the calculations if this intermediate step is stored by other functions (e.g., \code{\link{transfer}}). \code{do.phases} is passed to \code{\link{subcrt}} so that the properties of stable mineral phases as a function of temperature can optionally be calculated.

  \code{energy.args} is used by \code{affinity} to generate the argument list for \code{energy}. \code{affinity} passes the \code{...} list as \code{args} to \code{energy.args}, which returns an argument list appropriate for \code{energy}. \code{energy.args} also has the job of converting \samp{Eh} to \samp{pe} as a temperature-dependent function (see \code{\link{convert}}), and converting \samp{pe} and \samp{pH} to logarithms of activities of the electron and protein, respectively (i.e., negating the values).

  The \code{property} argument to \code{affinity} is analogous to the \code{what} argument of \code{energy}. Valid properties are \samp{A} or NULL for chemical affinity, \samp{logK} or \samp{logQ} for logarithm of equilibrium constant and reaction activity product, or any of the properties available in \code{\link{subcrt}} except for \samp{rho}. The properties returned are those of the formation reactions of the species of interest from the basis species. It is also possible to calculate the properties of the species of interest themselves (not their formation reactions) by setting the \code{property} to \samp{G.species}, \samp{Cp.species}, etc.

  Zero, one, or more leading arguments to the function identify which of the chemical activities of basis species, temperature, pressure and/or ionic strength to vary. The names of each of these arguments may be the formula of any of the basis species of the system, or \samp{T}, \samp{P}, \samp{pe}, \samp{pH}, \samp{Eh}, or \samp{IS} (but names may not be repeated). The value of each argument is of the form \code{c(min,max)} or \code{c(min,max,res)} where \code{min} and \code{max} refer to the minimimum and maximum values of variable identified by the name of the argument, and \code{res} denotes the resolution, or number of points along which to do the calculations (this is assigned a default value of 128 if it is missing, or 72 if \samp{T} and \samp{P} are both variables). For any arguments that specify the names of basis species, the numerical values refer to the logarithms of the activities of that basis species, or logarithms of fugacities if it is a gas. Unlike the \code{energy} function, the units of \samp{T} and \samp{P} in \code{affinity} are those the user has set using \code{\link{nuts}} (on program start-up these are \eqn{^{\circ}}{degrees }C and bar, respectively). 

  If proteins are in the species definition (which are distinguished from other species by having an underscore character in the name), and the basis species are charged, and \code{\link{thermo}$opt$ionize} is TRUE, \code{affinity} calls \code{\link{ionize}} to calculate the properties of charged proteins. If one or more buffers are assigned to the definition of \code{\link{basis}} species, \code{affinity} calls \code{\link{buffer}} to calculate the logarithms of activities of these basis species from the buffer.

  The commutative properties of the protein group additivity algorithm can be utilized in order to calculate the properties of many proteins in a fraction of the time it would take to calculate them individually. The \code{iprotein} and \code{logact.protein} arguments can be used to compute the chemical affinities of formation of proteins that are not in the global \code{\link{species}} definition. The appropriate \code{\link{basis}} species (e.g., \samp{CHNOSZ} or \samp{CHNOSZ+}) still must be defined prior to calling \code{affinity}. The user can give the indices of desired proteins contained in \code{\link{thermo}$protein} and \code{affinity} will automatically add to the species list the amino acid \emph{residues} and protein backbone group (found as \dQuote{proteins} in \code{\link{thermo}$protein}) then calculate the properties of the reactions for the residues, and add them together to get those of the indicated proteins. The \code{iprotein} option is compatible with calculations for \code{\link{ionize}}d proteins.

  The values of chemical affinity of formation reactions of the species of interest are returned in dimensionless (base 10) units (i.e., \emph{\strong{A}}/\eqn{2.303RT}).

  To use the names of charged basis species such as \samp{K+} and \samp{SO4-2} as the arguments, they should be enclosed in quotes (see the  example for for aluminum speciation in \code{\link{diagram}}).

}

\value{

  For \code{energy}, a list the first element of which is \code{sout} (the intermediate results from \code{\link{subcrt}}) and the second element of which is \code{a}, which contains the calculated properties. The latter itself is a list, one element for each species of interest,the elements of which have dimensions that are the number of \code{var}iables passed to the function.

  For \code{energy.arg}, a list with elements \code{what}, \code{vars}, \code{vals}, \code{lims}, \code{T}, \code{P}, \code{IS} that are each of the type for the corresponding arguments in \code{energy}.

  For \code{affinity}, a list, elements of which are \code{sout}, \code{property} (name of the calculated property), \code{basis} and \code{species} (definition of basis species and species of interest in effect at runtime), \code{T} and \code{P} (temperature and pressure, in the system units of Kelvin and bar, of length two (corresponding to minimum/maximum values) if either one is a variable of interest or length one otherwise), \code{xname} (the name of the first variable of interest, "" if none is present), \code{xlim} (if a first variable of interest is present, numeric of length 3 specifying the (minimum, maximum, resolution) of this variable), \code{yname} (the name of the second variable of interest, "" if none present), \code{ylim} (analogous to \code{xlim} but for a second variable of interest), \code{values}. The latter is itself a list, each element of which corresponds to a species of interest (names of the elements in this list are the character versions of the index of the species in \code{\link{thermo}$obigt}) and which contains a numeric value if  the number of variables of interest is zero or matrix otherwise. 

  If \code{pe} or \code{Eh}, or \code{pH} is/are among the variables of interest, \code{xnames} and/or \code{ynames} become \samp{e-} or \samp{H+} (respective to the property) and the minimum and maximum values in \code{xlim} and/or \code{ylim} adjusted accordingly (using \code{\link{convert}}).

}

\seealso{
 For prerequisites to \code{affinity}, see \code{\link{basis}} and \code{\link{species}}. Some calculations may invoke \code{\link{ionize}} and \code{\link{buffer}}. To visualize the results, use \code{\link{diagram}}.
}

\examples{
  \dontshow{data(thermo)}

  ## set up a system and calculate
  ## chemical affinities of formation reactions.
  basis(c('SiO2','MgO','H2O','O2'),c(-5,-5,0,999))
  species(c('quartz','enstatite','forsterite'))
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
  affinity(property='logK')
  # Standard molal Gibbs energies of species,
  # dimensionless (-G/2.303RT)
  affinity(property='G.species')
  # Standard molal Gibbs energies of species,
  # user dimensions (default: calories)
  affinity(property='G')

  ## Activity of glycine as a function of those of
  ## HCN and formaldehyde (200 C, 300 bar)
  ## After Schulte and Shock, 1995, Fig. 5
  # we can define the basis as this:
  basis(c('formaldehyde','H2O','HCN','O2'))
  species('glycine')
  t <- affinity(HCHO=c(-10,-2,9),HCN=c(-18,-2,9),T=200,P=300)
  # that gave us *affinities* (dimensionless) for logact(glycine)=-3
  # (the default). we can now find the *activities* that
  # correspond to affinity=0
  logact.glycine <- species()$logact + t$values[[1]]
  # note transposition of the z-value matrix in the following command
  contour(x=-10:-2,y=seq(-18,-2,by=2),z=t(logact.glycine),
    xlab=axis.label('HCHO'),ylab=axis.label('HCN'),
    labcex=1,xaxs='i',yaxs='i')
  title(main=paste('log activity glycine, after Schulte and Shock, 1995\n',
    describe(x=NULL),', logaH2O = 1'))

  ### amino acid synthesis at low and high temperatures
  ### after Amend and Shock, 1998
  # select the basis species and species of interest
  # and set their activities (first for the 18 degree C case)
  basis(c('H2O','CO2','NH4+','H2','H+','H2S'),
    rep('aq',6),log10(c(1,1e-4,5e-8,2e-9,5e-9,1e-15)))
  species(c('alanine','argininate','asparagine','aspartate','cysteine',
    'glutamate','glutamine','glycine','histidine','isoleucine',
    'leucine','lysinium','methionine','phenylalanine','proline',
    'serine','threonine','tryptophan','tyrosine','valine'),
    log10(c(3.9,0.7,1.1,3.3,0.5,3.8,1.0,5.8,1.2,0.7,
    0.8,1.0,2.8,0.5,0.5,4.6,5.8,0.6,0.9,2.8)/1e9))
  Tc <- 18
  T <- convert(Tc,'K')
  # converting A (dimensionless) to G of reaction (cal/mol) 
  # is like converting log K to standard G of reaction 
  AS98.18 <- 
    convert(convert(as.numeric(affinity(T=Tc)$values),'G',T=T),'J')/1000
  # the 100 degree C case
  Tc <- 100
  T <- convert(Tc,'K')
  basis(c('H2O','CO2','NH4+','H2','H+','H2S'),
    log10(c(1,2.2e-3,2.9e-6,3.4e-4,1.9e-6,1.6e-3)))
  species(1:20,log10(c(2.8e-9,5.0e-10,7.9e-10,2.4e-9,3.6e-10,
                       2.7e-9,7.2e-10,4.2e-9,8.6e-10,5.0e-10,
                       5.7e-10,7.2e-10,2.0e-9,3.6e-10,3.6e-10,
                       3.3e-9,4.2e-9,4.3e-10,6.5e-10,2.0e-9)))
  AS98.100 <- 
    convert(convert(as.numeric(affinity(T=Tc)$values),'G',T=T),'J')/1000
  # the nominal carbon oxidation state
  Z.C <- ZC(as.character(thermo$obigt$formula[thermo$species$ispecies]))
  # put them together
  print(data.frame(T100=AS98.100,T18=AS98.18,Z.C=Z.C))
  # values not exactly reproducing AS98 - different amino acid parameters

  ## affinities of metabolic reactions
  ## after Amend and Shock, 2001, Fig. 7
  \dontshow{
     # forget species from last example
     species(delete=TRUE)
  }
  basis(c('CO2','H2','NH3','O2','H2S','H+'))
  basis(c('O2','H2'),'aq')   # O2 and H2 were gas
  species('H2O')
  doplot <- function(T) {
    res <- 20
    t <- affinity(H2=c(-10,0,res),O2=c(-10,0,res),T=T)
    T.K <- convert(T,'K')   # temperature in Kelvin
    t <- convert(t$values[[1]],'G',T.K)  # affinities (cal/mol)
    t <- convert(t,'J')  # affinities (Joule)
    contour(x=seq(-10,0,length.out=res),
      y=seq(-10,0,length.out=res),z=t(t/1000),
      labcex=1,xlab=axis.label('H2'),ylab=axis.label('O2'))
  }
  layout(matrix(c(1,1,2,3,4,5),ncol=2,byrow=TRUE),heights=c(1,4,4))
  T <- c(25,55,100,150)
  par(mar=c(0,0,0,0))
  plot.new()
  text(0.5,0.1,paste(c('H2(aq) + 0.5O2(aq) = H2O(liq)\n\n',
    'after Amend and Shock, 2001')),cex=2)
  par(mar=c(3,3,0.5,0.5),cex=1.3,mgp=c(2,1,0))
  for(i in 1:length(T)) doplot(T[i])
  # this is so the plots in the next examples show up OK
  layout(matrix(1))

  ## continuation of last example, affinity calculations 
  ## in three dimensions
  print(affinity(H2=c(-10,0,3),O2=c(-10,0,3),T=c(25,150,4))$values)

}

\references{

  Amend, J. P. and Shock, E. L., 1998. Energetics of amino acid synthesis in hydrothermal ecosystems. \emph{Science}, 281, 1659-1662.

  Amend, J. P. and Shock, E. L., 2001. Energetics of overall metabolic reactions of thermophilic and hyperthermophilic Archaea and Bacteria. \emph{FEMS Microbiol. Rev.}, 25, 175-243.

  Kondepudi, D. K. and Prigogine, I., 1998. \emph{Modern Thermodynamics: From Heat Engines to Dissipative Structures}, John Wiley & Sons, New York, 486 p.

  Schulte, M. D. and Shock, E. L.,  1995. Thermodynamics of Strecker synthesis in hydrothermal systems. \emph{Orig. Life Evol. Biosph.}, 25, 161-173. 


}


\keyword{misc}
