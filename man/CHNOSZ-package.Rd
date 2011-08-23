\name{CHNOSZ-package}
\alias{CHNOSZ-package}
\docType{package}
\title{Chemical Thermodynamics and Activity Diagrams}
\description{
  CHNOSZ is a package for thermodynamic calculations, primarily with applications in geochemistry and biochemistry. It can be used to calculate the standard molal thermodynamic properties and chemical affinities of reactions relevant to geobiochemical processes, and to visualize the equilibrium activities of species on chemical speciation and predominance diagrams. The package can be used interactively and in batch mode, through the use of R source files containing a sequence of commands. The major features of the package are outlined below, with links to specific help topics in this document, which constitutes the primary technical description of the package. If you are a new user, the \sQuote{anintro} vignette (An introduction to CHNOSZ) may offer a more comfortable way to get started with using the package.
}

\details{

  Major features in CHNOSZ:

\itemize{

  \item Thermodynamic database - assembles literature values of the standard thermodynamic properties and equations of state parameters of minerals, aqueous organic and inorganic species, gases and liquids (\code{\link{thermo}}).

  \item Group additivity for proteins - estimate the standard thermodynamic properties and equations of state parameters for unfolded proteins from their amino acid composition; includes an additive calculation of ionization state of proteins as a function of temperature and pH (\code{\link{protein}}).

  \item File and internet access - read protein sequences from FASTA files, and download sequence information from UniProt (\code{\link{read.fasta}}, \code{\link{protein}}).

  \item Equations of state - calculate the standard thermodynamic properties of proteins or other species in the database, and reactions between them, as a function of temperature and pressure (\code{\link{hkf}}, \code{\link{cgl}}), \code{\link{subcrt}}.

  \item Stoichiometry - count elements in chemical formulas of species, check and optionally correct mass balance of chemical reactions (\code{\link{makeup}}).

  \item System of interest - define the basis species for a system together with one or more species of interest; compute the stoichiometries of the formation reactions of the species of interest (\code{\link{basis}}, \code{\link{species}}).

  \item Chemical affinity - calculate the chemical affinities of the formation reactions of the species of interest at a single point, or as a function of one or more of chemical activities of the basis species, temperature and/or pressure (\code{\link{affinity}}).

  \item Chemical activity - calculate the equilibrium activities of the species of interest as a function of the same variables used in the affinity calculation, using a reference state transformation (either the Boltzmann distribution or a reaction matrix approach). (\code{\link{diagram}}, \code{\link{equil.react}}, \code{\link{equil.boltz}}).

  \item Buffer calculations - compute activities of basis species that are determined by a buffer of one or more species (e.g., pyrite-pyrrhotite-magnetite; acetic acid-CO2) (\code{\link{buffer}}).

  \item Activity diagrams - plot the equilibrium activities at a single point (as barplots), or as a function of one (species activity diagrams) or two (predominance diagrams) variables (\code{\link{diagram}}).

  \item Activity statistics - calculate summary statistics for equilibrium activities of species (\code{\link{revisit}}).

  \item Multidimensional optimization (new in 0.9-3) - using an iterative gridded optimization, find a combination of chemical activities of basis species, temperature and/or pressure that maximize or minimize the value of a target statistic (\code{\link{findit}}).

  \item Mass transfer calculations (experimental) - calculate changes in solution composition and formation of secondary species as a function of incremental reaction of a mineral (or protein) (\code{\link{transfer}}).

}

  Here are some tips for new users:

  \itemize{
    \item Install the package from CRAN using \code{\link{install.packages}} or its GUI menu equivalent. 
    \item To begin working with the package after installation, type \code{\link{library}(CHNOSZ)} at the command line (or select the name of the package from the GUI menu).
    \item Running the examples shown in the various help topics is a great way to become more familiar with the usage of the functions. From \code{\link{help.start}}, select \sQuote{Packages} then \sQuote{CHNOSZ} and then select a function of interest. Individual examples can be run by pasting the example block directly into the R console.
    \item Type the command \code{\link{examples}()} to run all of the examples provided in CHNOSZ. This takes about five to ten minutes depending on your system. You should receive four warning messages when the examples are finished; don't worry, they are an expected result.
    \item Some of the examples require internet or file access or user intervention, or are intentionally written to demonstrate conditions that lead to errors. This offensive code is hidden from \R's package checking mechanism using the \code{dontrun} tag. You can experiment with \code{dontrun} examples by pasting the code to the R console.
    \item To learn how to update the thermodynamic database, look at its documentation in \code{\link{thermo}}.
  }

}


\section{Compatibility}{
The package depends on \R version 2.10.0 or greater (to read compressed data files). Before version 0.9-6 of the package, the dependency was given as \R version 2.7.0 or greater (major update of the X11 device in 2.7.0). However, without accessing the compressed data files in \code{extdata} it should be possible to run CHNOSZ under \R versions 2.4.0 or greater (availability of the \samp{stringsAsFactors} argument to \code{\link{data.frame}}). 
}

\section{Acknowledgements}{
  This package was made possible by the fearless leadership of the late Professor Harold C. Helgeson. He took no heed of disciplinary boundaries, and encouraged the author to pursue unconventional problems. Hal, and those in his research group, are in some way responsible for many parts of this package, such as the aqueous equations of state used here, the thermodynamic parameters for many species, and actual computer code (H2O92D.f, borrowed from the SUPCRT92 program).
  Work on this project at U.C. Berkeley (through 2008) was supported by research grants from the U.S. National Science Foundation and Department of Energy. In 2009--2011, research projects involving this package were supported by the National Science Foundation under grant EAR-0847616 (\url{http://www.nsf.gov/awardsearch/showAward.do?AwardNumber=0847616}). 
}

\examples{
  \dontshow{data(thermo)}
  ### Getting Started
  
  ## standard thermodynamic properties of species
  subcrt("H2O")
  subcrt("alanine")
  # names of proteins have an underscore
  subcrt("LYSC_CHICK")  
  # custom temperature range
  T <- seq(0,500,100)
  subcrt("H2O",T=T,P=1000)
  # temperature - pressure grid
  P <- seq(1000,4000,1000)
  subcrt("H2O",T=T,P=P,grid="P")
  
  ## information about species
  # query the database using formulas
  info("C6H12O6")
  info("SiO2")
  # query using names
  info("quartz")
  si <- info(c("glucose","mannose"))
  # show the equations of state parameters
  info(si)
  # approximate matches for names or formulas
  info("acid ")
  info("SiO2 ")

  ## standard thermodynamic properties of reactions
  # fermentation example
  info(c("fructose","ethanol"))
  subcrt(c("fructose","C2H5OH","CO2"),c(-1,2,2))
  # weathering example -- also see transfer()
  subcrt(c("k-feldspar","H2O","H+","kaolinite","K+","SiO2"),
    c(-2,-1,-2,1,2,4))
  # partial reaction auto-completion is possible
  basis(c("SiO2","H2O","K+","H+","O2"))
  subcrt(c("k-feldspar","kaolinite"),c(-2,1))
 
  ## chemical affinities
  # basis species
  basis(c("CO2","H2O","O2"))
  # species of interest
  species(c("CH4","C2H4O2","CO2"))
  # chemical affinities of formation reactions
  # take off $values for complete output
  affinity()$values
  affinity(O2=c(-90,-60,5))$values
}

\keyword{package}
