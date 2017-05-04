\encoding{UTF-8}
\name{CHNOSZ-package}
\alias{CHNOSZ-package}
\docType{package}
\title{Chemical Thermodynamics and Activity Diagrams}
\description{
CHNOSZ is a package for thermodynamic calculations, primarily with applications in geochemistry and compositional biology.
It can be used to calculate the standard molal thermodynamic properties and chemical affinities of reactions relevant to geobiochemical processes, and to visualize the equilibrium activities of species on chemical speciation and predominance diagrams.
}

\section{Getting Help}{
The major features of the package are outlined in the Overview given below, with links to specific help topics.
See the vignette \emph{An Introduction to CHNOSZ} (\code{anintro.Rmd}) for more in-depth example of the usage the functions.
The help pages have been given either keywords or \dQuote{concept index entries} which are visible to \code{\link{help.search}} (aka \code{??}).
Use \code{??primary} to browse the most commonly used functions and \code{??secondary} to see other high-level, but less often-used functions.
\code{??protein} shows functions for working with proteins, and \code{??extra} lists functions with extra functionality (beyond the main workflow).
For help on the thermodynamic database, use \code{??"thermodynamic data"}.
Finally, \code{??utilities} (one of the standard \R keywords) can be used to locate utility functions in the package; these include useful tools for modifying the database, converting units, reading protein sequence files, parsing chemical formulas, etc.
}

\section{Warning}{
All thermodynamic data and examples are provided on an as-is basis.
It is up to you to check not only the accuracy of the data, but also the \emph{suitability of the data AND computational techniques} for your problem.
By combining data taken from different sources, it is possible to build an inconsistent and/or nonsensical calculation.
An attempt has been made to provide a primary database (OBIGT) that is internally consistent, but no guarantee can be made.
If there is any doubt about the accuracy or suitability of data for a particular problem, please consult the primary sources (see \code{\link{thermo.refs}}).
}

\section{Overview}{

  Major features in CHNOSZ:

\itemize{

  \item Thermodynamic database - assembles literature values of the standard thermodynamic properties and equations of state parameters of minerals, aqueous organic and inorganic species, gases and liquids (\code{\link{thermo}}).

  \item Group additivity for proteins - estimate the standard thermodynamic properties and equations of state parameters for unfolded proteins from their amino acid composition; includes an additive calculation of ionization state of proteins as a function of temperature and pH (\code{\link{protein}}).

  \item File and internet access - read protein sequences from FASTA files, and download sequence information from UniProt (\code{\link{read.fasta}}, \code{\link{protein}}).

  \item Equations of state - calculate the standard thermodynamic properties of proteins or other species in the database, and reactions between them, as a function of temperature and pressure (\code{\link{hkf}}, \code{\link{cgl}}, \code{\link{subcrt}}).

  \item Stoichiometry - count elements in chemical formulas of species, check and optionally correct mass balance of chemical reactions (\code{\link{makeup}}).

  \item System of interest - define the basis species for a system together with one or more species of interest; compute the stoichiometries of the formation reactions of the species of interest (\code{\link{basis}}, \code{\link{species}}).

  \item Chemical affinity - calculate the chemical affinities of the formation reactions of the species of interest at a single point, or as a function of one or more of chemical activities of the basis species, temperature and/or pressure (\code{\link{affinity}}).

  \item Chemical activity - calculate the equilibrium activities of the species of interest as a function of the same variables used in the affinity calculation, using a reference state transformation (either the Boltzmann distribution or a reaction matrix approach). (\code{\link{diagram}}, \code{\link{equil.reaction}}, \code{\link{equil.boltzmann}}).

  \item Activity diagrams - plot the equilibrium activities at a single point (as barplots), or as a function of one (species activity diagrams) or two (predominance diagrams) variables (\code{\link{diagram}}).

  \item Buffer calculations - compute activities of basis species that are determined by a buffer of one or more species (e.g., pyrite-pyrrhotite-magnetite; acetic acid-\CO2) (\code{\link{buffer}}).

  \item Activity coefficients (\bold{experimental}) - calculate activity coefficients of aqueous species using the extended Debye-Hückel equation (\code{\link{nonideal}}).

  \item Activity statistics (\bold{experimental}) - calculate summary statistics for equilibrium activities of species (\code{\link{revisit}}).

  \item Multidimensional optimization (\bold{experimental}) - using an iterative gridded optimization, find a combination of chemical activities of basis species, temperature and/or pressure that maximize or minimize the value of a target statistic (\code{\link{findit}}).

  \item Mass transfer calculations (\bold{experimental}) - calculate changes in solution composition and formation of secondary species as a function of incremental reaction of a mineral (or protein) (\code{\link{transfer}}).

}

  Here are some tips for new users:

  \itemize{
    \item Running the examples shown in the various help topics is a good way to become more familiar with the usage of the functions. From \code{\link{help.start}}, select \sQuote{Packages} then \sQuote{CHNOSZ} and then select a function of interest. Individual examples can be run by pasting the example block directly into the R console.
    \item Type the command \code{\link{examples}()} to run all of the examples provided in CHNOSZ. This takes a few minutes depending on your system. If things go as expected, the entire set will run without any warnings or errors.
    \item Also try out \code{\link{demos}()} to run the demos (longer examples that are not part of the main help pages).
    \item To learn how to update the thermodynamic database, look at \code{\link{add.obigt}}.
  }

}

\section{Acknowledgements}{
  This package would not exist without the leadership and encouragement of Professor Harold C. Helgeson.
  Hal and his students and collaborators are in some way responsible for many of the equations and data contained in this package.
  A direct contribution of code is the file \code{H2O92D.f}, taken from the SUPCRT92 distribution, with minor modifications (masking of WRITE and STOP statements) made for compatibility with an \R environment and for allowing sub-zero degrees C calculations (keep valTP flag TRUE).
  The revised Helgeson-Kirkham-Flowers equations of state are used in this package, together with the thermodynamic properties and parameters for many species taken from numerous papers coauthored by Helgeson.

  Work on this package at U.C. Berkeley from ~2003--2008 was supported by research grants to HCH from the U.S. National Science Foundation and Department of Energy.
  In 2009--2011, development of this package was based upon work supported by the National Science Foundation under grant EAR-0847616.
  The files in \code{extdata/bison} are derived from BLAST calculations made on the Saguaro high performance computer at Arizona State University. 
}

\section{Known Bugs}{
  \code{\link{subcrt}} does not correctly identify the stable polymorph of some minerals at high temperature.

  The values generated by \code{\link{buffer}} may not be applied correctly by \code{\link{affinity}} in calculating the affinities of the formation reactions. (The values returned by \code{affinity(..., return.buffer=TRUE)} do appear to be correct in the examples).

  There is an unidentified inconsistency in \code{\link{transfer}} causing the reaction boundaries in one of the examples (\code{apc("closed")}) to be offset from the stability diagram.
  On the other hand, \code{feldspar("closed")} appears to work correctly.
}

\keyword{package}
