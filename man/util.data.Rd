\encoding{UTF-8}
\name{util.data}
\alias{util.data}
\alias{add.obigt}
\alias{mod.obigt}
\alias{browse.refs}
\alias{checkEOS}
\alias{checkGHS}
\alias{check.obigt}
\alias{obigt2eos}
\alias{today}
\alias{RH2obigt}
\title{Functions to Work with the Thermodynamic Database}
\description{
  Add species to or alter properties of species in the thermodynamic database or in the buffer definition table. Show references for sources of thermodynamic data in a web browser. Check internal consistency of individual entries in database.
}

\usage{
  add.obigt(file = system.file("extdata/thermo/OBIGT-2.csv", 
    package = "CHNOSZ"), force = FALSE, E.units = "cal")
  mod.obigt(...)
  today()
  browse.refs(key=NULL)
  checkEOS(eos, state, prop, ret.diff = FALSE)
  checkGHS(ghs, ret.diff = FALSE)
  check.obigt()
  obigt2eos(obigt, state, fixGHS = FALSE)
  RH2obigt(compound = NULL, state = "cr", 
    file = system.file("extdata/thermo/RH98_Table15.csv", package = "CHNOSZ"))
}

\arguments{
  \item{file}{character, path to a file}
  \item{force}{logical, force replacement of already existing species?}
  \item{E.units}{character, units of energy, \samp{cal} or \samp{J}}
  \item{...}{character or numeric, properties of species to modify in the thermodynamic database}
  \item{key}{character, numeric, or list, containing reference key(s) for which to show URL(s) in browser}
  \item{eos}{dataframe, containing equations-of-state parameters in the format of \code{thermo$obigt}}
  \item{state}{character, physical state of species}
  \item{prop}{character, property of interest (\samp{Cp} or \samp{V})}
  \item{ret.diff}{logical, return the difference between calculated and tabulated values?}
  \item{ghs}{dataframe, containing G, H and S, in the format of \code{thermo$obigt}}
  \item{obigt}{dataframe, in the format of \code{thermo$obigt}}
  \item{fixGHS}{logical, calculate one of missing G, H, or S?}
  \item{compound}{character, name of compound(s) in group additivity calculation}
}

\details{

  \code{add.obigt} reads data from the specified \code{file} and adds it to \code{\link{thermo}$obigt}. Set \code{force} to TRUE to replace species that exist in the thermodynamic database (each unique combination of a name and a state in the database is considered to be a unique species). \code{force}, if not specified, reverts to TRUE if the \code{file} is left at its default. Given the default setting of \code{E.units}, the function does not perform any unit conversions. If \code{E.units} is set to \samp{J}, then the thermodynamic parameters are converted from units of Joules to calories, as used in the CHNOSZ database. 

\code{mod.obigt} changes one or more of the properties of species or adds species to the thermodynamic database.
These changes are lost if you reload the database by calling \code{\link{data}(thermo)} or if you quit the \R session without saving it.
The name of the species to add or change must be supplied as the first argument of \code{...} or as a named argument (named \samp{name}).
When adding new species, a \samp{formula} should be included along with the values of any of the thermodynamic properties.
Additional arguments refer to the name of the property(s) to be updated and are matched to any part of compound column names in \code{\link{thermo}$obigt}, such as \samp{z} or \samp{T} in \samp{z.T}.
Unless \samp{state} is specified as one of the properties, its value is taken from \code{thermo$opt$state}.
When adding species, properties that are not specified become NA (except for \samp{state}).
The values provided should be in the units specifed in the documentation for the \code{thermo} data object, including any order-of-magnitude scaling factors.

  \code{today} returns the current date in the format adopted for \code{thermo$obigt} (inherited from SUPCRT-format data files) e.g. \samp{13.May.12} for May 13, 2012.

  \code{change} is a wrapper function to \code{mod.obigt} and \code{\link{mod.buffer}}. The \code{name} provided in the argument refers to the name or numeric index of the species to update or add using \code{mod.obigt}, unless the \code{name} begins with an underscore character, in which case the remaining part of the name (after the underscore) is passed to \code{mod.buffer}. The arguments in \code{...} are sent without change to the subordinate function.

  \code{browse.refs} with default arguments uses \code{\link{browseURL}} to display the sources of thermodynamic data in \code{thermo$refs}, with the URLs in that table showing as hyperlinks in the browser. Otherwise, if \code{key} is character, the URLs of those reference \code{key}s are opened in the browser. If \code{key} is numeric, the values refer to the species in those rows of \code{thermo$obigt}, and the URLs for each listed reference (\code{thermo$obigt$ref1}, \code{thermo$obigt$ref2}) are opened. If \code{key} is a list, it is interpreted as the result of a call to \code{subcrt}, and the reference URLs for each species involved in the calculation are opened.

  \code{checkEOS} compares heat capacity and volume calculated from equation-of-state parameters with reference (tabulated) values at 25 \eqn{^{\circ}}{°}C and 1 bar and prints a message and returns the calculated value if tolerance is exceeded. The Helgeson-Kirkham-Flowers equations of state parameters are in \code{eos}, which is a data frame with columns (and column names) in the same format as \code{\link{thermo}$obigt}. The \code{prop}erty can be one of \samp{Cp} or \code{V}. The code only distinguishes between \code{state}s of \samp{aq} and all others. The default tolerances, given in \code{thermo$opt$Cp.tol} and \code{thermo$opt$V.tol}, are 1 cal/K.mol for Cp and 1 cm3/mol for V. If \code{ret.diff} is TRUE, the differences are returned irrespective of their values, and no messages are printed.

  \code{checkGHS} compares G (standard molal Gibbs energy of formation from the elements) calculated from H (standard molal enthalpy of formation) and S (standard molal entropy) with reference (tabulated) values of G at 25 \eqn{^{\circ}}{°}C and 1 bar. A message is printed and the calculated difference is returned if it exceeds the value given in \code{thermo$opt$G.tol}, which has a default value of 100 cal/mol. The calculation requires that G, H and S, and the chemical formula of the species all be present. \code{checkEOS} and \code{checkGHS} are used by \code{\link{info}} when retrieving species parameters from the database.

  \code{check.obigt} is a function to check self-consistency of each entry (with itself) in the thermodynamic database, using \code{checkEOS} and \code{checkGHS}. The function checks data in both \code{\link{thermo}$obigt} and \code{\link{extdata}/thermo/OBIGT-2.csv}. The output is a table listing only species that exceed at least one of the tolerance limits, giving the name of the database (OBIGT or OBIGT-2), the species index (rownumber in the database), species name and state, and DCp, DV and DG, for the calculated differences (only those above the tolerances are given). This function is used to generate the file found at \code{extdata/thermo/obigt_check.csv}.

  \code{obigt2eos} takes a data frame in the format of \code{thermo$obigt} of one or more rows, removes scaling factors from equations-of-state parameters, and applies new column names depending on the \code{state}. This function is used by both \code{\link{info}} and \code{\link{subcrt}} when retrieving entries from the thermodynamic database. If \code{fixGHS} is TRUE a missing one of G, H or S for any species is calculated from the other two and the chemical formula of the species (used by \code{subcrt} when retrieving database entries).

  \code{RH2obigt} implements a group additivity algorithm for standard molal thermodynamic properties and equations of state parameters of crystalline and liquid organic molecules from Richard and Helgeson, 1998. The names of the \code{compound}s and their physical \code{state} are searched for in the indicated \code{file}, that also contains chemical formulas and group stoichiometries; the names of the groups are stored in the column names of this file, and must be present in \code{\link{thermo}$obigt}. The default \code{file} (\code{\link{extdata}/thermo/RH98_Table15.csv}) includes data taken from Table 15 of Richard and Helgeson, 1998 for high molecular weight compounds in \samp{cr}ystalline and \samp{liq}uid states. An error is produced if any of the \code{compound}-\code{state} combinations is not found in the \code{file}, if any of the group names for a given \code{compound}-\code{state} combination is not found in \code{thermo$obigt}, or if the chemical formula calculated from group additivity (with the aid of \code{\link{i2A}} and \code{\link{as.chemical.formula}}) is not identical to that listed in the \code{file}.

}

\section{Warning}{
  \code{add.obigt} affects the order of entries in \code{thermo}$obigt; therefore, it should be called before any \code{basis} or \code{species} definition. Also, there is no attempt made to keep the order of physical states in the database (aq-cr-liq-gas); the function simply adds new rows to the end of \code{thermo}$obigt. As a result, retrieving the properties of an added aqueous species using \code{\link{info}} requires an explicit \code{state="aq"} argument to that function if a species with the same name is in one of the cr, liq or gas states.
}

\value{
  The values returned (\code{\link{invisible}}-y) by \code{mod.obigt} are the rownumbers of the affected species.
}

\seealso{
  The default supplementary thermodynamic database for \code{add.obigt} (\code{\link{extdata}/thermo/OBIGT-2.csv}) is used in some of the example calculations in the help page for \code{\link{diagram}} and also in \code{\link{anim.TCA}}.

  \code{\link{mod.buffer}} for updating the list of available buffers.
}

\examples{
\dontshow{data(thermo)}
## modify an existing species with fake properties
ialanine <- mod.obigt("alanine", state="cr", G=0, H=0, S=0)
# we have made the values of G, H, and S inconsistent
# with the elemental composition of alanine so the following 
# produces a message about that
info(ialanine)
## add a species
inew <- mod.obigt("newspecies", formula="CHNOSZ", G=0, H=0, date=today())
info(inew)
\dontrun{
## marked dontrun because they open pages in a browser
# show the contents of thermo$refs
browse.refs()
# Internet needed for the following examples:
# open URL for Helgeson et al., 1998
browse.refs("HOK+98")
# open two URLs for alanine
browse.refs(info("alanine"))
# open three URLs for species in the reaction
s <- subcrt(c("O2","O2"),c("gas","aq"),c(-1,1))
browse.refs(s) }

## calculate thermodynamic properties of organic compounds
## using group additivity, after Richard and Helgeson, 1998
RH2obigt()
}

\references{

  Richard, L. and Helgeson, H. C. (1998) Calculation of the thermodynamic properties at elevated temperatures and pressures of saturated and aromatic high molecular weight solid and liquid hydrocarbons in kerogen, bitumen, petroleum, and other organic matter of biogeochemical interest. \emph{Geochim. Cosmochim. Acta} \bold{62}, 3591--3636. \url{http://dx.doi.org/10.1016/S0016-7037(97)00345-1}

}


\keyword{util}
