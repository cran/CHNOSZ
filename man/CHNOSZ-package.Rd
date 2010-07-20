\name{CHNOSZ}
\alias{CHNOSZ-package}
\alias{CHNOSZ}
\docType{package}
\title{Chemical Thermodynamics and Activity Diagrams}
\description{
  CHNOSZ is a software package for calculating the standard molal thermodynamic properties and chemical affinities of reactions in geobiochemical systems and for visualizing the equilibrium activities of species on chemical speciation and predominance diagrams.
}

\details{

  The functions in \pkg{CHNOSZ} perform three primary tasks: 1) calculation of the standard molal Gibbs energies and other thermodynamic properties of species and reactions. 2) calculation of the activity products of formation reactions of the species of interest from the basis species. The basis species stand for the descriptive compositional variables, i.e. perfectly mobile components (Korzhinskii, 1965). 3) combination of the results of (1) and (2) to generate diagrams showing the equilibrium distributions of species as a function of temperature, pressure, and chemical activities or fugacities of the basis species (which may include pH and Eh if the proton and electron are in the basis). 

  Here is some help on getting started with \pkg{CHNOSZ}:

  \itemize{

    \item Install the package from the command line with \samp{R CMD INSTALL CHNOSZ_0.9.tar.gz}. You need the proper build environment on your machine because there is a FORTRAN source file in the package. Or, if you have gotten the binary distribution (e.g., for Windows or Mac), use the \samp{Install packages from local zip files} option from the GUI.
    \item To begin working with the package after installation, type \code{\link{library}(CHNOSZ)} at the command line (or select the name of the package from the GUI menu).
    \item Type the command \code{\link{examples}()} to run all the examples provided in \code{CHNOSZ}. This takes about ten minutes depending on your system. Expect to get four warning messages when the examples are finished (or more if you don't have some other package required by the examples).
    \item The thermodynamic database, the system definition, and compositional data for proteins are stored in an object named \code{\link{thermo}}, which is added into the global namespace when you load the package. Running the examples or operating some of the functions may cause changes to this object; and correct operation of some of the examples demands that \code{thermo} is in its original state. The command \code{\link{data}(thermo)} is your friend, as it causes the program to reload the data object, effectively forgetting any changes made to it.
    \item To learn how to update the thermodynamic database, look at its documentation in \code{\link{thermo}}.
    \item To view the documentation, use \code{help.start()} (or choose \R help in the GUI) to launch a browser, and select \samp{Packages} then \samp{CHNOSZ}. You can copy the text of individual examples from the browser window into the command line interface. The documentation topics are arranged alphabetically; to read them in a more logical order, use this sequence (identical to the order of the \code{\link{examples}}): \code{\link{CHNOSZ}}, \code{\link{thermo}}, \code{util.*}, \code{\link{info}}, \code{eos}, \code{\link{water}}, \code{\link{subcrt}}, \code{\link{nuts}}, \code{\link{makeup}}, \code{\link{basis}}, \code{\link{species}}, \code{\link{affinity}}, \code{\link{diagram}}, \code{\link{buffer}}, \code{\link{protein}}, \code{\link{ionize}}, \code{\link{get.protein}}, \code{\link{revisit}}, \code{\link{transfer}}, \code{\link{eqdata}}.
    \item Parts of some of the examples require internet or file access or user intervention, or are intentionally written to demonstrate conditions that lead to errors. This offensive code is hidden from \R's package checking mechanism using the \code{dontrun} tag. The user can experiment with \code{dontrun} examples by pasting the code to the command line interface.
    
  }

}


\section{Compatibility}{
CHNOSZ depends on \R version 2.7.0 or greater. The reason for this dependency is that some of the examples that create plots cause errors in earlier versions of \R on Linux (there was a major update to the X11 device in 2.7.0). Except for that issue, CHNOSZ can run on \R versions 2.4.0 or greater. The reason for the latter dependency is availability of the \samp{stringsAsFactors} argument to \code{\link{data.frame}}. In some functions in this package, \samp{stringsAsFactors} is set to FALSE to facilitate manipulating the \code{\link{thermo}$obigt} and \code{thermo$protein} objects.
}

\section{Errata}{
As this package develops around the papers that support and use it, some errors and omissions in these papers have become apparent. The values of \eqn{a_2}{a2} and \eqn{a_4}{a4} for [-CH2NH2] were inadvertently set to zero in Table 6 of Dick et al., 2006. The correct values for these parameters are contained in the thermodynamic database distributed with \pkg{CHNOSZ}. There is a plotting error in Dick, 2009 that caused the colors of the points in Figs. 5 and 6 of that paper to be incorrectly assigned. This error is addressed in one of the examples for \code{\link{get.protein}}. These errors are relatively minor and do not change the conclusions of the respective papers.
}

\section{Acknowledgements}{
Work on this project at U.C. Berkeley (through 2008) was supported by research grants from the U.S. National Science Foundation and Department of Energy. Professor Harold C. Helgeson (d. May 2007) was principal investigator on these grants and my research advisor during the initial stages of development of this software. Effective in 2009: This material is based upon work supported by the National Science Foundation under grant EAR-0847616 (\url{http://www.nsf.gov/awardsearch/showAward.do?AwardNumber=0847616}). Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
}


\references{

  Dick, J. M., LaRowe, D. E. and Helgeson, H. C., 2006. Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences}, 3, 311-336. \url{http://www.biogeosciences.net/3/311/2006/bg-3-311-2006.html}

  Dick, J. M., 2009. Calculation of the relative metastabilities of proteins in subcellular compartments of \emph{Saccharomyces cerevisiae}. \emph{BMC Syst. Biol.}, 3, 75. \url{http://dx.doi.org/10.1186/1752-0509-3-75}

  Korzhinskii, D. S., 1965. The theory of systems with perfectly mobile components and processes of mineral formation. \emph{Am. J. Sci.}, 263, 193-205. \url{http://www.ajsonline.org/cgi/content/abstract/263/3/193}

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
  t <- info(c("glucose","mannose"))
  # show the equations of state parameters
  info(t)
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
