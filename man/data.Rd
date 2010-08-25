\name{thermo}
\docType{data}
\alias{thermo}
\alias{OBIGT}
\alias{proteins}
\alias{opt}
\alias{source}
\alias{stress}
\alias{PM90}
\alias{RH95}
\alias{RT71}
\alias{SGD}
\alias{ECO}
\alias{groups}
\alias{SOJSH}
\title{Thermodynamic Properties and Compositional Data}
\description{
 
  This data object holds the thermodynamic database of properties of species, along with operational parameters for \pkg{CHNOSZ}, the properties of elements, references to sources of thermodynamic and compositional data, compositions of chemical activity buffers, amino acid compositions of proteins, and miscellaneous other data taken from the literature. The \code{thermo} object also holds intermediate data used in calculations, in particular the definitions of basis species and species of interest input by the user, and the properties of \code{\link{water}} so that subsequent calculations at the same temperature-pressure conditions can be accelerated.

  The \code{thermo} object is a \code{\link{list}} composed of \code{\link{data.frame}}s or lists each representing a class of data. The object is created upon loading the package (by calling \code{\link{data}(thermo)} from within the \code{\link{.First.lib}} function) from \code{*.csv} files in the \code{data} directory of the package. \code{thermo} is globally accessible; i.e., it is present in the user's environment. After loading \pkg{CHNOSZ} you may run \code{ls()} to verify that \code{thermo} is present, or type \code{thermo} to print the entire contents of the object on the screen. The various elements of the \code{thermo} object can be accessed using \R's subsetting operators; for example, typing \code{thermo$opt} at the command line displays the current list of operational parameters (some of which can be altered using functions dedicated to this purpose; see e.g. \code{\link{nuts}}).

  To make persistent additions or changes to the thermodynamic database of your installation, including compositions of proteins, first locate the installation directory of the package. This will be different depending on your operating system and type of \R installation, but is something like /usr/lib/R/library/CHNOSZ, /Volumes/Macintosh HD/Library/Frameworks/R.framework/resources/ library/CHNOSZ, C:\\Program Files\\R\\R-2.10.0\\library\\CHNOSZ, or C:\\Users\\[User Name]\\Documents\\R \\win-library\\2.10\\CHNOSZ on Linux, Mac and Windows (XP and Vista) systems, respectively. To find the exact location of this directory on your system, use the command \code{\link{system.file}(package="CHNOSZ")}. Inside the \code{data} directory of the installation directory of the package are the \code{.csv} files that can be edited with a spreadsheet program. Edit and save the \code{OBIGT.csv} and/or \code{protein.csv} files as desired. The next time you start an R session, the new data will be available.

  Functions are available to interactively update the thermodynamic database or definitions of buffers (\code{\link{mod.obigt}} and \code{\link{mod.buffer}}, respectively; a function named \code{\link{change}} serves as a wrapper to both of these). Changes made using these functions, as well as any interactive definitions of basis species and species of interest, are lost when the current session is closed without saving or if the \code{thermo} object is reinitialized by the command \code{\link{data}(thermo)}.

}

\usage{data(thermo)}

\format{

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
      \code{ionize} \tab logical \tab Should \code{\link{affinity}} perform ionization calculations for proteins?\cr
      \code{water} \tab character \tab Computational option for properties of water (\samp{SUPCRT} (default) or \samp{IAPWS})\cr
      \code{online} \tab logical \tab Allow online searches of protein composition? Default (\code{NA}) is to ask the user.\cr
}

    \item \code{thermo$element}
    Dataframe containing the thermodynamic properties of elements taken from Cox et al., 1989 and Wagman et al., 1982. The standard molal entropy (\eqn{S}(\code{Z})) at 25 \eqn{^{\circ}}{degrees }C and 1 bar for the element of charge (\code{Z}) was calculated from \eqn{S}(H2,g) + 2\eqn{S}(\code{Z}) =  2\eqn{S}(H+), where the standard molal entropies of H2,g and H+ were taken from Cox et al., 1989. The mass of \code{Z} is taken to be zero. Accessing this dataframe using \code{\link{element}} will select the first entry found for a given element; i.e., values from Wagman et al., 1982 will only be retrieved if the properties of the element are not found from Cox et al., 1989.
      \tabular{lll}{
      \code{element}  \tab character  \tab Symbol of element\cr
      \code{state}  \tab character \tab Stable state of element at 25 \eqn{^{\circ}}{degrees }C and 1 bar\cr
      \code{source} \tab character \tab Source of data\cr
      \code{mass}  \tab numeric \tab Mass of element (in natural isotopic distribution;\cr
      \tab \tab referenced to a mass of 12 for \eqn{^{12}}{12}C)\cr
      \code{s}   \tab numeric \tab Entropy of the compound of the element in its stable\cr
      \tab \tab state at 25 \eqn{^{\circ}}{degrees }C and 1 bar (cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1})\cr
      \code{n}  \tab numeric \tab Number of atoms of the element in its stable\cr
      \tab \tab compound at 25 \eqn{^{\circ}}{degrees }C and 1 bar
    }

    \item \code{thermo$obigt}

  This dataframe is a thermodynamic database of standard molal thermodynamic properties and equations of state parameters of species. \acronym{OBIGT} is an acronym for \samp{OrganoBioGeoTherm}, which refers to a software package produced by Harold C. Helgeson and coworkers at the Laboratory of Theoretical Geochemistry and Biogeochemistry at the University of California, Berkeley. (There may be an additional meaning for the acronym: \dQuote{One BIG Table} of thermodynamic data.)

  As of \pkg{CHNOSZ} version 0.7, the data in \code{OBIGT.csv} represent 179 minerals, 16 gases, and 294 aqueous (largely inorganic) species taken from the data file included in the \acronym{SUPCRT92} distribution (Johnson et al., 1992), an additional 14 minerals, 6 gases, and 1049 aqueous organic and inorganic species from the \acronym{slop98.dat} file (Shock et al., 1998), and approximately 50 other minerals, 175 crystalline organic and biochemical species, 220 organic gases, 300 organic liquids, 650 aqueous inorganic, organic, and biochemical species, and 40 organic groups taken from the recent literature. Each entry is referenced to one or two literature sources listed in \code{thermo$source}. Some entries taken from the \acronym{SUPCRT92} or \acronym{slop98.dat} databases have been superseded, or duplicated, by later data that may be present in a separate file (OBIGT-2.csv; for more information see the help for \code{\link{danger}}). 

  Note the following modifications to some of the data:

  \itemize{
     \item The standard molal thermodynamic properties and equations of state parameters of the aqueous electron are zero except for the standard molal entropy at 25 \eqn{^{\circ}}{degrees }C and 1 bar, which is the opposite of that for the element of charge (\code{Z}, see above).
     \item The properties and parameters of some reference unfolded proteins used by Dick et al., 2006 are included here. Their names have dashes, instead of underscores, so that they are not confused with proteins whose properties are generated at runtime.
     \item The standard molal Gibbs energies and enthalpies of formation of the elements and entropies at 25 \eqn{^{\circ}}{degrees }C and 1 bar of aqueous metal-amino acid (alanate or glycinate) complexes reported by Shock and Koretsky, 1995 were recalculated by adding to their values the differences in the corresponding properties between the values for aqueous alanate and glycinate used by Shock and Koretsky, 1995, and those used by Amend and Helgeson, 1997b and Dick et al., 2006.
     \item The standard molal properties and equations-of-state parameters of four phase species (see below) of Fe(cr) were generated from heat capacity data given by Robie and Hemingway, 1995. 
  } 

  These modifications are indicated in \code{OBIGT.csv} by having \samp{CHNOSZ} as one of the sources of data. Note also that some data appearing in the \acronym{slop98.dat} file were corrected or modified as noted in that file, and are indicated in \code{OBIGT.csv} by having \samp{SLOP98} as one of the sources of data.

  In order to represent thermodynamic data for minerals with phase transitions, the different phases of these minerals are represented as phase species that have states denoted by \samp{cr1}, \samp{cr2}, etc. The standard molar thermodynamic properties at 25 \eqn{^{\circ}}{degrees }C and 1 bar (\eqn{T_r}{Pr} and \eqn{P_r}{Pr}) of the \samp{cr2} phase species of minerals were generated by first calculating those of the \samp{cr1} phase species at the transition temperature (\eqn{T_{tr}}{Ttr}) and 1 bar then taking account of the volume and entropy of transition (the latter can be retrieved by combining the former with the Clausius-Clapeyron equation and values of \eqn{(dP/dT)} of transitions taken from the \acronym{SUPCRT92} data file) to calculate the standard molar entropy of the \samp{cr2} phase species at \eqn{T_{tr}}{Ttr}, and taking account of the enthalpy of transition (\eqn{{\Delta}H^{\circ}}{DeltaH0}, taken from the \acronym{SUPCRT92} data file) to calculate the standard molar enthalpy of the \samp{cr2} phase species at \eqn{T_{tr}}{Ttr}. The standard molar properties of the \samp{cr2} phase species at \eqn{T_{tr}}{Ttr} and 1 bar calculated in this manner were combined with the equations-of-state parameters of the species to generate values of the standard molar properties at 25 \eqn{^{\circ}}{degrees }C and 1 bar. This process was repeated as necessary to generate the standard molar properties of phase species represented by \samp{cr3} and \samp{cr4}, referencing at each iteration the previously calculated values of the standard molar properties of the lower-temperature phase species (i.e., \samp{cr2} and \samp{cr3}). A consequence of tabulating the standard molar thermodynamic properties of the phase species is that the values of \eqn{(dP/dT)} and \eqn{{\Delta}H^{\circ}}{DeltaH0} of phase transitions can be calculated using the equations of state and therefore do not need to be stored in the thermodynamic database. However, the transition temperatures (\eqn{T_{tr}}{Ttr}) generally can not be assessed by comparing the Gibbs energies of phase species and are tabulated in the database.

  The identification of species and their standard molal thermodynamic properties at 25 \eqn{{^\circ}}{degrees }C and 1 bar are located in the first 12 columns of \code{thermo$obigt}:

    \tabular{lll}{
      \code{name}     \tab character \tab Species name\cr
      \code{abbrv}    \tab character \tab Species abbreviation\cr
      \code{formula}  \tab character \tab Species formula\cr
      \code{state}    \tab character \tab Physical state\cr
      \code{source1}  \tab character \tab Primary source\cr
      \code{source2}  \tab character \tab Secondary source\cr
      \code{date}     \tab character \tab Date of data entry\cr
      \code{G}        \tab numeric   \tab Standard molal Gibbs energy of formation\cr
      \tab \tab from the elements (cal mol\eqn{^{-1}}{^-1})\cr
      \code{H}        \tab numeric   \tab Standard molal enthalpy of formation\cr
      \tab \tab from the elements (cal mol\eqn{^{-1}}{^-1})\cr
      \code{S}        \tab numeric   \tab Standard molal entropy (cal mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1})\cr
      \code{Cp}       \tab numeric   \tab Standard molal isobaric heat capacity (cal mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1})\cr
      \code{V}	      \tab numeric   \tab Standard molal volume (cm\eqn{^3} mol\eqn{^{-1}}{^-1})
    }


   The meanings of the remaining columns depend on the physical state of a particular species. If it is aqueous, the values in these columns represent parameters in the revised HKF equations of state (see \code{\link{hkf}}), otherwise they denote parameters in a general equation of state for crystalline, gas and liquid species (see \code{\link{cgl}}). The names of these columns are compounded from those of the parameters in each of the equations of state (for example, column 13 is named \code{a1.a}). Scaling of the values by orders of magnitude is adopted for some of the parameters, following common usage in the literature.

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
    Dataframe of references to sources of thermodynamic data. Source keys with a leading underscore indicate abbreviations for journals.
    \tabular{lll}{
      \code{source} \tab character \tab Source key\cr
      \code{reference} \tab character \tab Reference\cr
    }

    \item \code{thermo$buffer}

    Dataframe which contains definitions of buffers of chemical activity. Each named buffer can be composed of one or more species, which may include any species in the thermodynamic database and/or any protein. The calculations provided by \code{\link{buffer}} do not take into account phase transitions of minerals, so individual phase species of such minerals must be specified in the buffers.
    \tabular{lll}{
      \code{name} \tab character \tab Name of buffer\cr
      \code{species} \tab character \tab Name of species\cr
      \code{state} \tab character \tab Physical state of species\cr
      \code{logact} \tab numeric \tab Logarithm of activity (fugacity for gases)
    }

    \item \code{thermo$protein}
    Dataframe of amino acid compositions of selected proteins. The majority of the compositions were taken from the SWISS-PROT online database (Boeckmann et al., 2003). N-terminal signal sequences were removed except for some cases where different isoforms of proteins have been identified (for example, \samp{MOD5.M} and \code{MOD5.N} proteins of \samp{YEAST} denote the mitochondrial and nuclear isoforms of this protein.) 
    \tabular{lll}{
      \code{protein} \tab character \tab Identification of protein\cr
      \code{organism} \tab character \tab Identification of organism\cr
      \code{source} \tab character \tab Source of compositional data\cr
      \code{abbrv} \tab character \tab Abbreviation or other ID for protein\cr
      \code{chains} \tab numeric \tab Number of polypeptide chains in the protein\cr
      \code{Ala}\dots\code{Tyr} \tab numeric \tab Number of each amino acid in the protein
    }

    \item \code{thermo$stress}
    Dataframe listing proteins identified in selected proteomic stress response experiments. The names of proteins begin at row 3, and columns are all the same length (padded as necessary at the bottom by \code{NA}s). Names correspond to ordered locus names (for \samp{SGD}) or gene names (for \samp{ECO}). The column names and first two rows give the following information:
    \tabular{lll}{
      \code{colname} \tab character \tab Name of the experiment\cr
      \code{organism} \tab character \tab Name of the organism (\samp{SGD} or \samp{ECO})\cr
      \code{source} \tab character \tab Source of the data\cr
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

    \item \code{thermo$expt} 
    List of experimental data that are used in the examples.

      \describe{
           \code{thermo$expt$PM90} Heat capacities of four unfolded aqueous proteins taken from Privalov and Makhatadze, 1990. Names of proteins are in the first column, temperature in \eqn{^{\circ}}{degrees }C in the second, and heat capacities in J mol\eqn{^{-1}}{^-1} K\eqn{^{-1}}{^-1} in the third.

           \code{thermo$expt$RH95} Heat capacity data for iron taken from Robie and Hemingway, 1995. Temperature in Kelvin is in the first column, heat capacity in J K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} in the second.

           \code{thermo$expt$RT71} pH titration measurements for unfolded lysozyme (\samp{LYSC_CHICK}) taken from Roxby and Tanford, 1971. pH is in the first column, net charge in the second.

           \code{thermo$expt$SOJSH92} Experimental equilibrium constants for the reaction NaCl(aq) = Na+ + Cl- as a function of temperature and pressure taken from Fig. 1 of Shock et al., 1992. Data were extracted from the figure using g3data (\url{http://www.frantz.fi/software/g3data.php}).

      }

    \item \code{thermo$SGD}
    Dataframe of amino acid composition of proteins from the \emph{Saccharomyces} Genome Database.

    \describe{
      Contains twenty-two columns. Values in the first column are the rownumbers, the second column (\code{OLN}) has the ordered locus names of proteins, and the remaining twenty columns (\code{Ala}..\code{Val}) contain the numbers of the respective amino acids in each protein; the columns are arranged in alphabetical order based on the three-letter abbreviations for the amino acids. The source of data for \samp{SGD.csv} is the file \samp{protein_properties.tab} found on the FTP site of the SGD project on 2008-08-04. Blank entries were replaced with "NA" and column headings were added.

    }

    \item \code{thermo$ECO}

    \describe{
      Contains 24 columns. Values in the first column correspond to rownumbers, the second column {\code{AC}} holds the accession numbers of the proteins, the third column (\code{Name}) has the names of the corresponding genes, and the fourth column {\code{OLN}} lists the ordered locus names of the proteins. The remaining twenty columns (\code{A}..\code{Y}) give the numbers of the respective amino acids in each protein and are ordered alphabetically by the one-letter abbreviations of the amino acids. The sources of data for \samp{ECO.csv} are the files \samp{ECOLI.dat} \url{ftp://ftp.expasy.org/databases/hamap/complete_proteomes/entries/bacteria} and \samp{ECOLI.fas} \url{ftp://ftp.expasy.org/databases/hamap/complete_proteomes/fasta/bacteria} downloaded from the HAMAP (High-quality Automated and Manual Annotation of microbial Proteomes system) FTP site (Gattiker et al., 2003) on 2007-12-20.
    }


    \item \code{thermo$yeastgfp}

    \describe{
      Has 28 columns; the names of the first five are \code{yORF}, \code{gene name}, \code{GFP tagged?}, \code{GFP visualized?}, and \code{abundance}. The remaining columns correspond to the 23 subcellular localizations considered in the YeastGFP project (Huh et al., 2003 and Ghaemmaghami et al., 2003) and hold values of either \code{T} or \code{F} for each protein. \samp{yeastgfp.csv} was downloaded on 2007-02-01 from http://yeastgfp.ucsf.edu using the Advanced Search, setting options to download the entire dataset and to include localization table and abundance, sorted by orf number.
    }
     
  }  % end of itemize

} % end of format

\seealso{ \code{\link{add.protein}} and \code{\link{add.obigt}} for adding data from local .csv files.
}

\examples{
  \dontshow{data(thermo)}
  ## exploring thermo$obigt
  # what physical states there are
  unique(thermo$obigt$state)
  # formulas of ten random species
  n <- nrow(thermo$obigt)
  thermo$obigt$formula[runif(10)*n]

  ## make a table of duplicated species
  name <- thermo$obigt$name
  state <- thermo$obigt$state
  source <- thermo$obigt$source1
  species <- paste(name,state)
  dups <- species[which(duplicated(species))]
  id <- numeric()
  for(i in 1:length(dups)) id <- c(id,which(species \%in\% dups[i]))
  data.frame(name=name[id],state=state[id],source=source[id])
}

\references{

  Amend, J. P. and Helgeson, H. C., 1997b. Calculation of the standard molal thermodynamic properties of aqueous biomolecules at elevated temperatures and pressures. Part 1. L-\eqn{\alpha}{alpha}-amino acids. \emph{J. Chem. Soc., Faraday Trans.}, 93, 1927-1941. \url{http://dx.doi.org/10.1039/a608126f}

  Cox, J. D., Wagman, D. D. and Medvedev, V. A., eds., 1989. \emph{CODATA Key Values for Thermodynamics}. Hemisphere Publishing Corporation, New York, 271 p. \url{http://www.worldcat.org/oclc/18559968}

  Dick, J. M., LaRowe, D. E. and Helgeson, H. C., 2006. Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences}, 3, 311-336. \url{http://www.biogeosciences.net/3/311/2006/bg-3-311-2006.html}

  Gattiker, A., Michoud, K., Rivoire, C., Auchincloss, A. H., Coudert, E., Lima, T., Kersey, P., Pagni, M., Sigrist, C. J. A., Lachaize, C., Veuthey, A.-L., Gasteiger, E. and Bairoch, A., 2003. Automatic annotation of microbial proteomes in Swiss-Prot. \emph{Comput. Biol. Chem.}, 27, 49-58. \url{http://dx.doi.org/10.1016/S1476-9271(02)00094-4}

  Ghaemmaghami, S., Huh, W., Bower, K., Howson, R. W., Belle, A., Dephoure, N., O'Shea, E. K. and Weissman, J. S., 2003. Global analysis of protein expression in yeast. \emph{Nature}, 425, 737-741. \url{http://dx.doi.org/10.1038/nature02046}

  HAMAP system. HAMAP FTP directory, \url{ftp://ftp.expasy.org/databases/hamap/}, accessed on 2007-12-20.

  Huh, W. K., Falvo, J. V., Gerke, L. C., Carroll, A. S., Howson, R. W., Weissman, J. S. and O'Shea, E. K., 2003. Global analysis of protein localization in budding yeast. \emph{Nature}, 425, 686-691. \url{http://dx.doi.org/10.1038/nature02026}

  Johnson, J. W., Oelkers, E. H. and Helgeson, H. C., 1992. SUPCRT92: A software package for calculating the standard molal thermodynamic properties of minerals, gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000\eqn{^{\circ}}{degrees }C. \emph{Comp. Geosci.}, 18, 899-947. \url{http://dx.doi.org/10.1016/0098-3004(92)90029-Q}

  Privalov, P. L. and Makhatadze, G. I., 1990. Heat capacity of proteins. II. Partial molar heat capacity of the unfolded polypeptide chain of proteins: Protein unfolding effects. \emph{J. Mol. Biol.}, 213, 385-391. \url{http://dx.doi.org/10.1016/S0022-2836(05)80198-6}

  Robie, R. A. and Hemingway, B. S., 1995. \emph{Thermodynamic Properties of Minerals and Related Substances at 298.15 K and 1 Bar (\eqn{10^5} Pascals) Pressure and at Higher Temperatures}. U. S. Geol. Surv., Bull. 2131, 461 p. \url{http://www.worldcat.org/oclc/32590140}

  Roxby, R. and Tanford, C., 1971. Hydrogen ion titration curve of lysozyme in 6 M guanidine hydrochloride. \emph{Biochemistry}, 10, 3348-3352. \url{http://dx.doi.org/10.1021/bi00794a005}

  SGD project. \emph{Saccharomyces} Genome Database, \url{http://www.yeastgenome.org}, accessed on 2008-08-04.

  Shock, E. L. and Koretsky, C. M., 1995. Metal-organic complexes in geochemical processes: Estimation of standard partial molal thermodynamic properties of aqueous complexes between metal cations and monovalent organic acid ligands at high pressures and temperatures. \emph{Geochim. Cosmochim. Acta}, 59, 1497-1532. \url{http://dx.doi.org/10.1016/0016-7037(95)00058-8}

  Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C., 1992. Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \eqn{^{\circ}}{degrees }C and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.}, 88, 803-826. \url{http://dx.doi.org/10.1039/FT9928800803}

  Shock, E. L. et al., 1998. slop98.dat (computer data file). \url{http://geopig.asu.edu/supcrt92_data/slop98.dat}, accessed on 2005-11-05.

  Wagman, D. D., Evans, W. H., Parker, V. B., Schumm, R. H., Halow, I., Bailey, S. M., Churney, K. L. and Nuttall, R. L., 1982. The NBS tables of chemical thermodynamic properties. Selected values for inorganic and C\eqn{_1}{1} and C\eqn{_2}{2} organic substances in SI units. \emph{J. Phys. Chem. Ref. Data}, 11 (supp. 2), 1-392. \url{http://www.nist.gov/srd/PDFfiles/jpcrdS2Vol11.pdf}

  YeastGFP project. Yeast GFP Fusion Localization Database, http://yeastgfp.ucsf.edu, accessed on 2007-02-01. Current location: \url{http://yeastgfp.yeastgenome.org}

}

\keyword{datasets}
