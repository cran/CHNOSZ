\name{CHNOSZ-package}
\alias{CHNOSZ-package}
\docType{package}
\title{Chemical Thermodynamics and Activity Diagrams}
\description{
  CHNOSZ is a package for thermodynamic calculations, primarily with applications in geochemistry and biochemistry. It can be used to calculate the standard molal thermodynamic properties and chemical affinities of reactions relevant to geobiochemical processes, and to visualize the equilibrium activities of species on chemical speciation and predominance diagrams. The package can be used interactively -- variations on a single command line can be used to calculate the standard thermodynamic properties of any number of reactions involving proteins and/or inorganic species such as minerals -- and in batch mode, through the use of R source files containing a sequence of commands. The major features of the package are outlined below, with links to specific help topics in this document, which constitutes the primary technical description of the package. If you are a new user, the \sQuote{anintro} vignette (An Introduction to CHNOSZ) may offer a more comfortable way to get started with using the package.
}

\details{

  Major features in CHNOSZ:

\itemize{

  \item Thermodynamic database - assembles literature values of the standard thermodynamic properties and equations of state parameters of minerals, aqueous organic and inorganic species, gases and liquids (\code{\link{thermo}}).

  \item Group additivity for proteins - estimate the standard thermodynamic properties and equations of state parameters for unfolded proteins from their amino acid composition; includes an additive calculation of ionization state of proteins as a function of temperature and pH (\code{\link{protein}}).

  \item File and internet access - read protein sequences from FASTA files, and download sequence information from UniProt (\code{\link{read.fasta}}, \code{\link{protein}}).

  \item Equations of state - calculate the standard thermodynamic properties of proteins or other species in the database, and reactions between them, as a function of temperature and pressure (\code{\link{hkf}}, \code{\link{cgl}}).

  \item Stoichiometry - count elements in chemical formulas of species, check and optionally correct mass balance of chemical reactions (\code{\link{makeup}}).

  \item System of interest - define the basis species for a system together with one or more species of interest; compute the stoichiometries of the formation reactions of the species of interest (\code{\link{basis}}, \code{\link{species}}).

  \item Chemical affinity - calculate the chemical affinities of the formation reactions of the species of interest at a single point, or as a function of one or more of chemical activities of the basis species, temperature and/or pressure (\code{\link{affinity}}).

  \item Chemical activity - calculate the equilibrium activities of the species of interest as a function of the same variables used in the affinity calculation, using a reference state transformation (either the Boltzmann distribution or a reaction matrix approach). These are often (in organic systems particularly) metastable equilibrium activities (i.e., the species of interest are less stable than the basis species). (\code{\link{diagram}}, \code{\link{equil.react}}, \code{\link{equil.boltz}}).

  \item Buffer calculations - compute the values of the activities of basis species that are determined by a buffer of one or more species (e.g., pyrite-pyrrhotite-magnetite; acetic acid-CO2) (\code{\link{buffer}}).

  \item Activity diagrams - plot the metastable equilibrium activities at a single point (as barcharts), or as a function of one (species activity diagrams) or two (predominance diagrams) variables (\code{\link{diagram}}).

  \item Activity statistics - calculate statistical quantities for metastable equilibrium activities of species. Examples of summary statistics: standard deviation, coefficient of variation, correlation with a (log)normal distribution, richness (number of species with activities above a defined minimum). Can also use a set of reference values in order to calculate root mean square deviation (RMSD), Pearson or Spearman correlation coefficients (\code{\link{revisit}}).

  \item Multidimensional optimization (new in 0.9-3) - using an iterative gridded optimization, find a combination of chemical activities of basis species, temperature and/or pressure that maximize or minimize the value of a target statistic (above) (\code{\link{findit}}).

  \item Mass transfer calculations (experimental) - calculate changes in solution composition and formation of secondary species as a function of incremental reaction of a mineral (or protein) (\code{\link{transfer}}).

}

  Here are some tips for new users:

  \itemize{

    \item Install the package from the command line with \samp{R CMD INSTALL CHNOSZ_0.9-6.tar.gz}. You need the proper build environment on your machine because there is a FORTRAN source file in the package. Or, if you have gotten the binary distribution (e.g., for Windows or Mac), use the \samp{Install packages from local zip files} option from the GUI.
    \item To begin working with the package after installation, type \code{\link{library}(CHNOSZ)} at the command line (or select the name of the package from the GUI menu).
    \item Running the examples shown in the various help topics is a great way to become more familiar with the usage of the functions. Type the command \code{\link{examples}()} to run all the examples provided in \code{CHNOSZ}. This takes about ten minutes depending on your system. You should receive some warning messages when the examples are finished; don't worry, they mean the functions are functioning as intended.
    \item To view the documentation, use \code{help.start()} (or choose \R help in the GUI) to launch a browser, and select \samp{Packages} then \samp{CHNOSZ}. You can copy the text of individual examples from the browser window into the command line interface. The documentation topics are arranged alphabetically; to read them in a more logical order, use this sequence (identical to the order of the \code{\link{examples}}): \code{\link{CHNOSZ-package}}, \code{\link{thermo}}, \code{util.*}, \code{\link{info}}, \code{eos}, \code{\link{water}}, \code{\link{subcrt}}, \code{\link{nuts}}, \code{\link{makeup}}, \code{\link{basis}}, \code{\link{species}}, \code{\link{affinity}}, \code{\link{diagram}}, \code{\link{buffer}}, \code{\link{protein}}, \code{\link{ionize}}, \code{\link{get.protein}}, \code{\link{revisit}}, \code{\link{transfer}}, \code{\link{eqdata}}.
    \item Parts of some of the examples require internet or file access or user intervention, or are intentionally written to demonstrate conditions that lead to errors. This offensive code is hidden from \R's package checking mechanism using the \code{dontrun} tag. The user can experiment with \code{dontrun} examples by pasting the code to the command line interface.
    \item To learn how to update the thermodynamic database, look at its documentation in \code{\link{thermo}}.
    \item The thermodynamic database, the system definition, and compositional data for proteins are stored in an object named \code{\link{thermo}}, which is added into the global search space when you load the package. The command \code{\link{data}(thermo)} is your friend, as it causes the program to reinitialize the object. This is useful for quickly "starting over" and in some cases it is the only effective solution to errors that may occur when configuring a calculation.
    
  }

}


\section{Compatibility}{
The package dependency is \R version 2.10.0 or greater (to read compressed data files). Before this, the package dependency was given as \R version 2.7.0 or greater (major update of the X11 device in 2.7.0). However, without accessing the compressed data files in \code{extdata} it should be possible to run CHNOSZ under \R versions 2.4.0 or greater (availability of the \samp{stringsAsFactors} argument to \code{\link{data.frame}}). 
}

\section{Acknowledgements}{
  This package was made possible by the fearless leadership of the late Professor Harold C. Helgeson. He took no heed of disciplinary boundaries, and encouraged the author to pursue unconventional problems. Hal, and those in his research group, are in some way responsible for many parts of this package, such as the aqueous equations of state used here, the thermodynamic parameters for many species, and actual computer code (H2O92D.f, borrowed from the SUPCRT92 program).
  Work on this project at U.C. Berkeley (through 2008) was supported by research grants from the U.S. National Science Foundation and Department of Energy. Effective starting in 2009: This material is based upon work supported by the National Science Foundation under grant EAR-0847616 (\url{http://www.nsf.gov/awardsearch/showAward.do?AwardNumber=0847616}). 
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
