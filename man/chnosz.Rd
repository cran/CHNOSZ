\name{CHNOSZ}
\alias{CHNOSZ}
\title{Chemical Thermodynamics and Activity Diagrams}
\description{
  CHNOSZ is a software package for calculating the standard molal thermodynamic properties and chemical affinities of reactions in geobiochemical systems and for visualizing the stable or metastable equilibrium reference state on chemical speciation and predominance diagrams.
}

\usage{
  CHNOSZ()
}

\details{

  \pkg{CHNOSZ} is a package for thermodynamic calculations and graphical analysis.

  The functions in \pkg{CHNOSZ} perform three primary tasks: 1) calculation of the standard molal Gibbs energies and other thermodynamic properties of species and reactions. 2) calculation of the activity products of formation reactions of the species of interest from the basis species. The basis species stand for the descriptive compositional variables, i.e. perfectly mobile components (Korzhinskii, 1965). 3) combination of the results of (1) and (2) to generate diagrams showing the equilibrium distributions of species as a function of temperature, pressure, and chemical activities or fugacities of the basis species (which may include pH and Eh if the proton and electron are in the basis). 

  Note: The CHNOSZ function exists so that this overview page can be linked to the rest of the documentation. All this function does is run the \dQuote{getting started} examples below. Read further for the functions that really do the work.

  Here is some help on getting started with \pkg{CHNOSZ}:

  \itemize{

    \item Install the package from the command line with \samp{R CMD INSTALL CHNOSZ\_0.8.tar.gz}. You need the proper build environment on your machine because there is a FORTRAN source file in the package. Or, if you have gotten the binary distribution (e.g., for Windows or Mac), use the \samp{Install packages from local zip files} option from the GUI.
    \item To begin working with the package after installation, type \code{library(CHNOSZ)} at the command line (or select the name of the package from the GUI menu).
    \item Type the command \code{\link{examples}()} to run through all the examples provided in \code{CHNOSZ}. This takes about ten minutes depending on your system. Expect to get a few warning messages when the examples are finished. 
    \item The thermodynamic database, the system definition, and compositional data for proteins are stored in an object named \samp{thermo} (see \code{\link{thermo}}), which is added into the global namespace when you load the package. Running the examples or operating some of the functions may cause changes to this object; and correct operation of the examples generally demands that \code{thermo} is in its original state. Even trying to run the same example twice in succession can sometimes produce strange results. The command \code{data(thermo)} is your friend, as it causes the program to reload the data object, effectively forgetting any changes made to it.
    \item To learn how to update the thermodynamic database, look at its documentation in \code{\link{thermo}}.
    \item To view the documentation, use \code{help.start()} (or choose \R help in the GUI) to launch a browser, and select \samp{Packages} then \samp{CHNOSZ}. You can copy the text of individual examples from the browser window into the command line interface. The documentation topics are arranged alphabetically; to read them in the intended order, use this sequence (identical to the order of the \code{\link{examples}}): \code{\link{CHNOSZ}}, \code{\link{thermo}}, \code{\link{utilities}}, \code{\link{info}}, \code{eos}, \code{\link{water}}, \code{\link{subcrt}}, \code{\link{nuts}}, \code{\link{makeup}}, \code{\link{basis}}, \code{\link{species}}, \code{\link{affinity}}, \code{\link{diagram}}, \code{\link{buffer}}, \code{\link{protein}}, \code{\link{ionize}}, \code{\link{get.protein}}, \code{\link{diversity}}, \code{\link{transfer}}.
    \item Parts of some of the examples require internet or file access or user intervention, or are intentionally written to demonstrate conditions that lead to errors. This offensive code is hidden from \R's package checking mechanism using the \code{dontrun} tag. The user can experiment with \code{dontrun} examples by pasting the code to the command line interface.
    
  }

  Acknowledgement: This material is based upon work supported by the National Science Foundation under grant EAR-0847616. 

}

\references{

  Korzhinskii, D. S., 1965. The theory of systems with perfectly mobile components and processes of mineral formation. \emph{Am. J. Sci.}, 263, 193-205.


}

\examples{
  \dontshow{data(thermo)}

  ### Getting Started
  
  ## standard thermodynamic properties of species
  subcrt('H2O')
  subcrt('alanine')
  subcrt('LYSC_CHICK')  # a protein
  # custom temperature range
  T <- seq(0,500,100)
  subcrt('H2O',T=T,P=1000)
  # temperature - pressure grid
  P <- seq(1000,4000,1000)
  subcrt('H2O',T=T,P=P,grid='P')
  
  ## information about species
  # query the database using formulas
  info('C6H12O6')
  info('SiO2')
  # query using names
  info('quartz')
  t <- info(c('glucose','mannose'))
  # show the equations of state parameters
  info(t)
  # approximate matches for names or formulas
  info('acid ')
  info('SiO2 ')

  ## standard thermodynamic properties of reactions
  # fermentation example
  info(c('fructose','ethanol'))
  subcrt(c('fructose','C2H5OH','CO2'),c(-1,2,2))
  # weathering example -- also see transfer()
  subcrt(c('k-feldspar','H2O','H+','kaolinite','K+','SiO2'),
    c(-2,-1,-2,1,2,4))
  # partial reaction auto-completion is possible
  basis(c('SiO2','H2O','K+','H+','O2'))
  subcrt(c('k-feldspar','kaolinite'),c(-2,1))
 
  ## chemical affinities
  # basis species
  basis(c('CO2','H2O','O2'))
  # species of interest
  species(c('CH4','C2H4O2','CO2'))
  # chemical affinities of formation reactions
  # take off $values for complete output
  affinity()$values
  affinity(O2=c(-90,-60,5))$values
  

}


\keyword{package}
