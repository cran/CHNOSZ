\encoding{UTF-8}
\name{thermo}
\docType{data}
\alias{thermo}
\alias{OBIGT}
\alias{opt}
\alias{element}
\alias{refs}
\alias{groups}
\alias{buffers}
\title{Thermodynamic Database and System Definition}
\description{

  Please read the important Warning at \code{\link{CHNOSZ-package}}.
 
  The core data files provided with CHNOSZ are in the \code{data} directory of the package. These \code{*.csv} files are used to build the \code{thermo} data object on loading the package. Additional (extra) data files, supporting the examples and vignettes, are documented separately at \code{\link{extdata}}. 

  The \code{thermo} object holds the thermodynamic database of properties of species, some thermodynamic constants and operational parameters for functions in CHNOSZ, the properties of elements, references to literature sources of thermodynamic data, compositions of chemical activity buffers, and amino acid compositions of proteins. The \code{thermo} object also holds intermediate data used in calculations, in particular the definitions of basis species and species of interest input by the user, and the calculated properties of \code{\link{water}} so that subsequent calculations at the same temperature-pressure conditions can be accelerated.

The \code{thermo} object is a \code{\link{list}} composed of \code{\link{data.frame}}s or lists each representing a class of data.
The object is created in an environment named \code{thermo}; see \code{\link{sideeffects}} for details.
It is created upon loading the package, through a call to \code{data(thermo)} in \code{\link{.onAttach}}.
At any time, the user can restore the data object to its initial state by calling \code{data(thermo)}.
This is sometimes a useful command to use during an interactive session, when previous definitions of basis species and species of interest are longer desired.

  The function \code{\link{add.obigt}} is available to update the thermodynamic database in use in a running session. For example, one can run \code{add.obigt("mydata.csv")} after loading the package, and the data in that file will be added to the database. The format of this file must be the same as the \code{OBIGT.csv} file provided with CHNOSZ. Although changes made using \code{add.obigt} are lost when the current \R session is closed, the data can always be restored the next time as long as the user has the \code{mydata.csv} (or other) file available.

  The first example below shows how to find the installation locations of \code{OBIGT.csv} and other \code{*.csv} files. Making changes to these files is not recommended, because incompatible changes can leave the package unusable; also, the files will be overwritten whenver the package is installed (or updated). Instead, use these files as templates for creating your own database files.

}

\usage{data(thermo)}

\format{

  The items in the \code{thermo} object are documented below.

  \itemize{
     
    \item \code{thermo$opt} 
    List of operational parameters 
    \tabular{lll}{
      \code{Tr} \tab numeric \tab Reference temperature (K)\cr 
      \code{Pr} \tab numeric \tab Reference pressure (bar)\cr 
      \code{Theta} \tab numeric \tab \eqn{\Theta}{Theta} in the revised HKF equations of state (K)\cr 
      \code{Psi} \tab numeric \tab \eqn{\Psi}{Psi} in the revised HKF equations of state (bar)\cr
      \code{cutoff} \tab numeric \tab Cutoff below which values are taken to be zero (see \code{\link{makeup}})\cr
      \code{E.units} \tab character \tab The user's units of energy (\samp{cal} (default) or \samp{J})\cr
      \code{T.units} \tab character \tab The user's units of temperature (\samp{C} (default) or \samp{K})\cr
      \code{P.units} \tab character \tab The user's units of pressure (\samp{bar} (default) or \samp{MPa})\cr
      \code{state} \tab character \tab The default physical state for searching species (\samp{aq} by default)\cr
      \code{water} \tab character \tab Computational option for properties of water (\samp{SUPCRT} (default) or \samp{IAPWS})\cr
      \code{online} \tab logical \tab Allow online searches of protein composition? Default (\code{NA}) is to ask the user.\cr
      \code{G.tol} \tab numeric \tab Difference in value of G above which \code{\link{checkGHS}} produces a message (cal mol\eqn{^{-1}}{^-1}).\cr
      \code{Cp.tol} \tab numeric \tab Difference in value of Cp above which \code{\link{checkEOS}} produces a message (cal K\eqn{^{-1}}{-1} mol\eqn{^{-1}}{^-1}).\cr
      \code{V.tol} \tab numeric \tab Difference in value of V above which \code{\link{checkEOS}} produces a message (cm\eqn{^{3}}{3} mol\eqn{^{-1}}{^-1}).\cr
      \code{varP} \tab logical \tab Use variable-pressure standard state for gases (default \code{FALSE})? See \code{\link{subcrt}}.\cr
}

    \item \code{thermo$element}
    Dataframe containing the thermodynamic properties of elements taken from Cox et al., 1989 and Wagman et al., 1982. The standard molal entropy (\eqn{S}(\code{Z})) at 25 \eqn{^{\circ}}{°}C and 1 bar for the element of charge (\code{Z}) was calculated from \eqn{S}(H2,g) + 2\eqn{S}(\code{Z}) =  2\eqn{S}(H+), where the standard molal entropies of H2,g and H+ were taken from Cox et al., 1989. The mass of \code{Z} is taken to be zero. Accessing this data frame using \code{\link{mass}} or \code{\link{entropy}} will select the first entry found for a given element; i.e., values from Wagman et al., 1982 will only be retrieved if the properties of the element are not found from Cox et al., 1989.
      \tabular{lll}{
      \code{element}  \tab character  \tab Symbol of element\cr
      \code{state}  \tab character \tab Stable state of element at 25 \eqn{^{\circ}}{°}C and 1 bar\cr
      \code{source} \tab character \tab Source of data\cr
      \code{mass}  \tab numeric \tab Mass of element (in natural isotopic distribution;\cr
      \tab \tab referenced to a mass of 12 for \eqn{^{12}}{12}C)\cr
      \code{s}   \tab numeric \tab Entropy of the compound of the element in its stable\cr
      \tab \tab state at 25 \eqn{^{\circ}}{°}C and 1 bar (cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1})\cr
      \code{n}  \tab numeric \tab Number of atoms of the element in its stable\cr
      \tab \tab compound at 25 \eqn{^{\circ}}{°}C and 1 bar
    }

    \item \code{thermo$obigt}

  This dataframe is a thermodynamic database of standard molal thermodynamic properties and equations of state parameters of species. \acronym{OBIGT} is an acronym for \samp{OrganoBioGeoTherm}, which refers to a software package produced by Harold C. Helgeson and coworkers at the Laboratory of Theoretical Geochemistry and Biogeochemistry at the University of California, Berkeley. (There may be an additional meaning for the acronym: \dQuote{One BIG Table} of thermodynamic data.)

  As of CHNOSZ version 0.7, the data in \code{OBIGT.csv} represent 179 minerals, 16 gases, and 294 aqueous (largely inorganic) species taken from the data file included in the \acronym{SUPCRT92} distribution (Johnson et al., 1992), an additional 14 minerals, 6 gases, and 1049 aqueous organic and inorganic species from the \acronym{slop98.dat} file (Shock et al., 1998), and approximately 50 other minerals, 175 crystalline organic and biochemical species, 220 organic gases, 300 organic liquids, 650 aqueous inorganic, organic, and biochemical species, and 40 organic groups taken from the recent literature. Each entry is referenced to one or two literature sources listed in \code{thermo$refs}. Use \code{\link{browse.refs}} to display the references in a browser window.

  Note the following additional modifications:

  \itemize{
     \item Use corrected values of \eqn{a_2}{a2} and \eqn{a_4}{a4} for [-CH2NH2] (were incorrectly listed as zero in Table 6 of Dick et al., 2006).
     \item The standard molal thermodynamic properties and equations of state parameters of the aqueous electron are zero except for the standard molal entropy at 25 \eqn{^{\circ}}{°}C and 1 bar, which is the opposite of that for the element of charge (\code{Z}, see above).
     \item The properties and parameters of some reference unfolded proteins used by Dick et al., 2006 are included here. Their names have dashes, instead of underscores, so that they are not confused with proteins whose properties are generated at runtime.
     \item The standard molal Gibbs energies and enthalpies of formation of the elements and entropies at 25 \eqn{^{\circ}}{°}C and 1 bar of aqueous metal-amino acid (alanate or glycinate) complexes reported by Shock and Koretsky, 1995 were recalculated by adding to their values the differences in the corresponding properties between the values for aqueous alanate and glycinate used by Shock and Koretsky, 1995, and those used by Amend and Helgeson, 1997b and Dick et al., 2006.
     \item The standard molal properties and equations-of-state parameters of four phase species (see below) of Fe(cr) were generated from heat capacity data given by Robie and Hemingway, 1995. 
  } 

  These modifications are indicated in \code{OBIGT.csv} by having \samp{CHNOSZ} as one of the sources of data. Note also that some data appearing in the \acronym{slop98.dat} file (Shock et al., 1998) were corrected or modified as noted in that file, and are indicated in \code{OBIGT.csv} by having \samp{SLOP98} as one of the sources of data.

  In order to represent thermodynamic data for minerals with phase transitions, the different phases of these minerals are represented as phase species that have states denoted by \samp{cr1}, \samp{cr2}, etc. The standard molar thermodynamic properties at 25 \eqn{^{\circ}}{°}C and 1 bar (\eqn{T_r}{Pr} and \eqn{P_r}{Pr}) of the \samp{cr2} phase species of minerals were generated by first calculating those of the \samp{cr1} phase species at the transition temperature (\eqn{T_{tr}}{Ttr}) and 1 bar then taking account of the volume and entropy of transition (the latter can be retrieved by combining the former with the Clausius-Clapeyron equation and values of \eqn{(dP/dT)} of transitions taken from the \acronym{SUPCRT92} data file) to calculate the standard molar entropy of the \samp{cr2} phase species at \eqn{T_{tr}}{Ttr}, and taking account of the enthalpy of transition (\eqn{{\Delta}H^{\circ}}{DeltaH0}, taken from the \acronym{SUPCRT92} data file) to calculate the standard molar enthalpy of the \samp{cr2} phase species at \eqn{T_{tr}}{Ttr}. The standard molar properties of the \samp{cr2} phase species at \eqn{T_{tr}}{Ttr} and 1 bar calculated in this manner were combined with the equations-of-state parameters of the species to generate values of the standard molar properties at 25 \eqn{^{\circ}}{°}C and 1 bar. This process was repeated as necessary to generate the standard molar properties of phase species represented by \samp{cr3} and \samp{cr4}, referencing at each iteration the previously calculated values of the standard molar properties of the lower-temperature phase species (i.e., \samp{cr2} and \samp{cr3}). A consequence of tabulating the standard molar thermodynamic properties of the phase species is that the values of \eqn{(dP/dT)} and \eqn{{\Delta}H^{\circ}}{DeltaH0} of phase transitions can be calculated using the equations of state and therefore do not need to be stored in the thermodynamic database. However, the transition temperatures (\eqn{T_{tr}}{Ttr}) generally can not be assessed by comparing the Gibbs energies of phase species and are tabulated in the database.

  The identification of species and their standard molal thermodynamic properties at 25 \eqn{{^\circ}}{°}C and 1 bar are located in the first 12 columns of \code{thermo$obigt}:

    \tabular{lll}{
      \code{name}     \tab character \tab Species name\cr
      \code{abbrv}    \tab character \tab Species abbreviation\cr
      \code{formula}  \tab character \tab Species formula\cr
      \code{state}    \tab character \tab Physical state\cr
      \code{ref1}  \tab character \tab Primary source\cr
      \code{ref2}  \tab character \tab Secondary source\cr
      \code{date}     \tab character \tab Date of data entry (formatted as in SUPCRT92)\cr
      \code{G}        \tab numeric   \tab Standard molal Gibbs energy of formation\cr
      \tab \tab from the elements (cal mol\eqn{^{-1}}{^-1})\cr
      \code{H}        \tab numeric   \tab Standard molal enthalpy of formation\cr
      \tab \tab from the elements (cal mol\eqn{^{-1}}{^-1})\cr
      \code{S}        \tab numeric   \tab Standard molal entropy (cal mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1})\cr
      \code{Cp}       \tab numeric   \tab Standard molal isobaric heat capacity (cal mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1})\cr
      \code{V}	      \tab numeric   \tab Standard molal volume (cm\eqn{^3} mol\eqn{^{-1}}{^-1})
    }


   The meanings of the remaining columns depend on the physical state of a particular species. If it is aqueous, the values in these columns represent parameters in the revised HKF equations of state (see \code{\link{hkf}}), otherwise they denote parameters in a general equations for crystalline, gas and liquid species (see \code{\link{cgl}}). The names of these columns are compounded from those of the parameters in each of the equations of state (for example, column 13 is named \code{a1.a}). Scaling of the values by orders of magnitude is adopted for some of the parameters, following common usage in the literature.

    Columns 13-20 for aqueous species (parameters in the revised HKF equations of state):
    \tabular{lll}{
      \code{a1} \tab numeric \tab \eqn{a_1\times10}{a1 * 10} (cal mol\eqn{^{-1}}{^-1} bar\eqn{^{-1}}{^-1})\cr
      \code{a2} \tab numeric \tab \eqn{a_2\times10^{-2}}{a2 * 10^{-2}} (cal mol\eqn{^{-1}}{^-1})\cr
      \code{a3} \tab numeric \tab \eqn{a_3}{a3} (cal K mol\eqn{^{-1}}{^-1} bar\eqn{^{-1}}{^-1})\cr
      \code{a4} \tab numeric \tab \eqn{a_4\times10^{-4}}{a4 * 10^-4} (cal mol\eqn{^{-1}}{^-1} K)\cr
      \code{c1} \tab numeric \tab \eqn{c_1}{c1} (cal mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1})\cr
      \code{c2} \tab numeric \tab \eqn{c_2\times10^{-4}}{c2 * 10^-4} (cal mol\eqn{^{-1}}{^-1} K)\cr
      \code{omega} \tab numeric \tab \eqn{\omega\times10^{-5}}{omega * 10^-5} (cal mol\eqn{^{-1}}{^-1})\cr
      \code{Z}  \tab numeric \tab Charge
 
    }

    Columns 13-20 for crystalline, gas and liquid species (\eqn{Cp=a+bT+cT^{-2}+dT^{-0.5}+eT^2+fT^{\lambda}}{Cp = a + bT + cT^-2 + dT^-0.5 + eT^2 + fT^lambda}).
    \tabular{lll}{
      \code{a} \tab numeric \tab \eqn{a} (cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1})\cr
      \code{b} \tab numeric \tab \eqn{b\times10^3}{b * 10^3} (cal K\eqn{^{-2}}{^-2} mol\eqn{^{-1}}{^-1})\cr
      \code{c} \tab numeric \tab \eqn{c\times10^{-5}}{c * 10^-5} (cal K mol\eqn{^{-1}}{^-1})\cr
      \code{d} \tab numeric \tab \eqn{d} (cal K\eqn{^{-0.5}}{^-0.5} mol\eqn{^{-1}}{^-1})\cr
      \code{e} \tab numeric \tab \eqn{e\times10^5}{e * 10^5} (cal K\eqn{^{-3}}{^-3} mol\eqn{^{-1}}{^-1})\cr
      \code{f} \tab numeric \tab \eqn{f} (cal K\eqn{^{-\lambda-1}}{-lambda-1} mol\eqn{^{-1}}{^-1})\cr
      \code{lambda} \tab numeric \tab \eqn{\lambda}{lambda} (exponent on the \eqn{f} term)\cr
      \code{T} \tab numeric \tab Temperature of phase transition or upper\cr
      \tab \tab temperature limit of validity of extrapolation (K)
    }

    \item \code{thermo$source}
    Dataframe of references to sources of thermodynamic data.
    \tabular{lll}{
      \code{key} \tab character \tab Source key\cr
      \code{author} \tab character \tab Author(s)\cr
      \code{year} \tab character \tab Year\cr
      \code{citation} \tab character \tab Citation (journal title, volume, and article number or pages; or book or report title)\cr
      \code{URL} \tab character \tab URL\cr
    }

    \item \code{thermo$buffers}

    Dataframe which contains definitions of buffers of chemical activity. Each named buffer can be composed of one or more species, which may include any species in the thermodynamic database and/or any protein. The calculations provided by \code{\link{buffer}} do not take into account phase transitions of minerals, so individual phase species of such minerals must be specified in the buffers.
    \tabular{lll}{
      \code{name} \tab character \tab Name of buffer\cr
      \code{species} \tab character \tab Name of species\cr
      \code{state} \tab character \tab Physical state of species\cr
      \code{logact} \tab numeric \tab Logarithm of activity (fugacity for gases)
    }

    \item \code{thermo$protein}
    Data frame of amino acid compositions of selected proteins. Many of the compositions were taken from the SWISS-PROT/UniProt online database (Boeckmann et al., 2003) and the protein and organism names usually follow the conventions adopted there. In some cases different isoforms of proteins are identified using modifications of the protein names; for example, \samp{MOD5.M} and \code{MOD5.N} proteins of \samp{YEAST} denote the mitochondrial and nuclear isoforms of this protein. See \code{\link{iprotein}} to search this data frame by protein name, and other functions to work with the amino acid compositions.
    \tabular{lll}{
      \code{protein} \tab character \tab Identification of protein\cr
      \code{organism} \tab character \tab Identification of organism\cr
      \code{ref} \tab character \tab Reference key for source of compositional data\cr
      \code{abbrv} \tab character \tab Abbreviation or other ID for protein\cr
      \code{chains} \tab numeric \tab Number of polypeptide chains in the protein\cr
      \code{Ala}\dots\code{Tyr} \tab numeric \tab Number of each amino acid in the protein
    }

    \item \code{thermo$groups}
    This is a dataframe with 22 columns for the amino acid sidechain, backbone and protein backbone groups ([Ala]..[Tyr],[AABB],[UPBB]) whose rows correspond to the elements C, H, N, O, S. It is used to quickly calculate the chemical formulas of proteins that are selected using the \code{iprotein} argument in \code{\link{affinity}}.

    \item \code{thermo$basis}
    Initially \code{NULL}, reserved for a dataframe written by \code{\link{basis}} upon definition of the basis species. The number of rows of this dataframe is equal to the number of columns in \dQuote{...} (one for each element).
     \tabular{lll}{
        \code{...} \tab numeric \tab One or more columns of stoichiometric\cr
        \tab \tab coefficients of elements in the basis species\cr
        \code{ispecies} \tab numeric \tab Rownumber of basis species in \code{thermo$obigt}\cr
        \code{logact} \tab numeric \tab Logarithm of activity or fugacity of basis species\cr
        \code{state} \tab character \tab Physical state of basis species\cr
     }

    \item \code{thermo$species}
    Initially \code{NULL}, reserved for a dataframe generated by \code{\link{species}} to define the species of interest. The number of columns in \dQuote{...} is equal to the number of basis species (i.e., rows of \code{thermo$basis}).
    \tabular{lll}{
       \code{...} \tab numeric \tab One or more columns of stoichiometric\cr
       \tab \tab coefficients of basis species in the species of interest\cr
       \code{ispecies} \tab numeric \tab Rownumber of species in \code{thermo$obigt}\cr
       \code{logact} \tab numeric \tab Logarithm of activity or fugacity of species\cr
       \code{state} \tab character \tab Physical state of species\cr
       \code{name} \tab character \tab Name of species\cr
    }

    \item \code{thermo$water}
    The properties calculated with \code{\link{water}} at multiple T, P points (minimum of 26) are stored here so that repeated calculations at the same conditions can be done more quickly.

    \item \code{thermo$Psat}
    The values of Psat calculated with \code{water.SUPCRT} at multiple T points (minimum of 26) are stored here.

    \item \code{thermo$water2}
    The properties calculated with \code{water.SUPCRT} at multiple T, P points (minimum of 26) are stored here.


  }  % end of itemize with long descriptions
     

} % end of format

\seealso{ 
  \code{\link{add.obigt}} for thermodynamic data from local .csv files.
}

\examples{
\dontshow{data(thermo)}
## where are OBIGT.csv and the other data 
## files on your installation?
system.file("data",package="CHNOSZ")

## exploring thermo$obigt
# what physical states there are
unique(thermo$obigt$state)
# formulas of ten species at random
n <- nrow(thermo$obigt)
thermo$obigt$formula[runif(10)*n]
}

\references{

  Amend, J. P. and Helgeson, H. C. (1997b) Calculation of the standard molal thermodynamic properties of aqueous biomolecules at elevated temperatures and pressures. Part 1. L-\eqn{\alpha}{alpha}-amino acids. \emph{J. Chem. Soc., Faraday Trans.} \bold{93}, 1927--1941. \url{http://dx.doi.org/10.1039/a608126f}

  Cox, J. D., Wagman, D. D. and Medvedev, V. A., eds. (1989) \emph{CODATA Key Values for Thermodynamics}. Hemisphere Publishing Corporation, New York, 271 p. \url{http://www.worldcat.org/oclc/18559968}

  Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 3110--336. \url{http://www.biogeosciences.net/3/311/2006/bg-3-311-2006.html}

  Johnson, J. W., Oelkers, E. H. and Helgeson, H. C. (1992) SUPCRT92: A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000\eqn{^{\circ}}{°}C. \emph{Comp. Geosci.} \bold{18}, 899--947. \url{http://dx.doi.org/10.1016/0098-3004(92)90029-Q}

  Shock, E. L. and Koretsky, C. M. {1995} Metal-organic complexes in geochemical processes: Estimation of standard partial molal thermodynamic properties of aqueous complexes between metal cations and monovalent organic acid ligands at high pressures and temperatures. \emph{Geochim. Cosmochim. Acta} \bold{59}, 1497--1532. \url{http://dx.doi.org/10.1016/0016-7037(95)00058-8}

  Shock, E. L. et al. {1998} \emph{SLOP98.dat} (computer data file). http://geopig.asu.edu/supcrt92_data/slop98.dat, accessed on 2005-11-05. Current location: \url{http://geopig.asu.edu/sites/default/files/slop98.dat}.

  Wagman, D. D., Evans, W. H., Parker, V. B., Schumm, R. H., Halow, I., Bailey, S. M., Churney, K. L. and Nuttall, R. L. (1982) The NBS tables of chemical thermodynamic properties. Selected values for inorganic and C\eqn{_1}{1} and C\eqn{_2}{2} organic substances in SI units. \emph{J. Phys. Chem. Ref. Data} \bold{11} (supp. 2), 1--392. \url{http://www.nist.gov/data/PDFfiles/jpcrdS2Vol11.pdf}

}

\concept{Thermodynamic data}
