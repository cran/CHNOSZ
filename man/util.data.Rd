\name{util.data}
\alias{util.data}
\alias{add.obigt}
\alias{mod.obigt}
\alias{add.protein}
\alias{change}
\alias{browse.refs}
\alias{checkEOS}
\alias{checkGHS}
\alias{check.obigt}
\alias{obigt2eos}
\title{Functions to Work with the Thermodynamic Database}
\description{
  Add species to or alter properties of species in the thermodynamic database or in the buffer definition table. Show references for sources of thermodynamic data in a web browser. Check internal consistency of individual entries in database.
}

\usage{
  add.obigt(file = system.file("extdata/thermo/OBIGT-2.csv", 
    package = "CHNOSZ"), force = FALSE, E.units = "cal")
  mod.obigt(species, ..., missingvalues = NA)
  add.protein(file="protein.csv")
  change(name, ...)
  browse.refs(key=NULL)
  checkEOS(eos, state, prop, ret.diff = FALSE)
  checkGHS(ghs, ret.diff = FALSE)
  check.obigt()
  obigt2eos(obigt, state, fixGHS = FALSE)
}

\arguments{
  \item{file}{character, path to a file.}
  \item{force}{logical, force replacement of already existing species?}
  \item{E.units}{character, units of energy, \samp{cal} or \samp{J}.}
  \item{species}{character, formulas or names of species to modify or add to the thermodynamic database.}
  \item{...}{character or numeric, properties of species to modify in the thermodynamic database.}
  \item{missingvalues}{numeric, values to assign to undefined properties.}
  \item{name}{character or numeric, name (or numeric index) of species or name of buffer to be modified.}
  \item{key}{character, numeric, or list, containing reference key(s) for which to show URL(s) in browser.}
  \item{eos}{dataframe, containing equations-of-state parameters in the format of \code{thermo$obigt}.}
  \item{state}{character, physical state of species (\samp{aq} or other).}
  \item{prop}{character, property of interest (\samp{Cp} or \samp{V}).}
  \item{ret.diff}{logical, return the difference between calculated and tabulated values?}
  \item{ghs}{dataframe, containing G, H and S, in the format of \code{thermo$obigt}.}
  \item{obigt}{dataframe, in the format of \code{thermo$obigt}.}
  \item{fixGHS}{logical, calculate one of missing G, H, or S?}
}

\details{

  \code{add.obigt} and \code{add.protein} read data from the specified \code{file} and add it to either \code{\link{thermo}$obigt} or \code{thermo$protein}, as appropriate. Set \code{force} to TRUE to replace species that exist in the thermodynamic database (each unique combination of a name and a state in the database is considered to be a unique species). \code{force}, if not specified, reverts to TRUE if the \code{file} is left at its default. Given the default setting of \code{E.units}, the function does not perform any unit conversions. If \code{E.units} is set to \samp{J}, then the thermodynamic parameters are converted from units of Joules to calories, as used in the CHNOSZ database. 


  \code{mod.obigt} changes one or more of the properties of one or more species or adds species to the thermodynamic database. These changes are lost if you reload the database by calling \code{\link{data}(thermo)} or if you quit the \R session without saving it. To modify the properties of species, give the names in the \code{species} argument and supply other arguments: if one of these arguments is \code{state}, species in those states will be updated. Additional arguments refer to the name of the property(s) to be updated and correspond to the column names of \code{\link{thermo}$obigt} (the names of the properties are matched to any part of compound column names, such as \samp{z.T}). The values provided should be in the units specifed in the documentation for the \code{thermo} data object. To add species, supply the new names in \code{species} and provide an argument named \code{formula} with the corresponding chemical formulas. Additional arguments refer to any of the properties you wish to specify. Properties that are not specified are assigned the value of \code{missingvalues} which is \code{NA} by default (however if \code{state} is missing it is set to the value of \code{thermo$opt$state}). The values returned (\code{\link{invisible}}-y) by \code{mod.obigt} are the rownumbers of the affected species.

  \code{change} is a wrapper function to \code{mod.obigt} and \code{\link{mod.buffer}}. The \code{name} provided in the argument refers to the name or numeric index of the species to update or add using \code{mod.obigt}, unless the \code{name} begins with an underscore character, in which case the remaining part of the name (after the underscore) is passed to \code{mod.buffer}. The arguments in \code{...} are sent without change to the subordinate function.

  \code{browse.refs} with default arguments uses \code{\link{browseURL}} to display the sources of thermodynamic data in \code{thermo$refs}, with the URLs in that table showing as hyperlinks in the browser. Otherwise, if \code{key} is character, the URLs of those reference \code{key}s are opened in the browser. If \code{key} is numeric, the values refer to the species in those rows of \code{thermo$obigt}, and the URLs for each listed reference (\code{thermo$obigt$ref1}, \code{thermo$obigt$ref2}) are opened. If \code{key} is a list, it is interpreted as the result of a call to \code{subcrt}, and the reference URLs for each species involved in the calculation are opened.

  \code{checkEOS} compares heat capacity and volume calculated from equation-of-state parameters with reference (tabulated) values at 25 degrees C and 1 bar and prints a message and returns the calculated value if tolerance is exceeded. The Helgeson-Kirkham-Flowers equations of state parameters are in \code{eos}, which is a data frame with columns (and column names) in the same format as \code{\link{thermo}$obigt}. The \code{prop}erty can be one of \samp{Cp} or \code{V}. The code only distinguishes between \code{state}s of \samp{aq} and all others. The tolerances are 1 cal/K.mol for Cp and 1 cm3/mol for V. If \code{ret.diff} is TRUE, the differences are returned irrespective of their values, and no messages are printed.

  \code{checkGHS} compares G (standard molal Gibbs energy of formation from the elements) calculated from H (standard molal enthalpy of formation) and S (standard molal entropy) with reference (tabulated) values of G at 25 degrees C and 1 bar. A message is printed and the calculated difference is returned if it exceeds a tolerance of 500 cal/mol. The calculation requires that G, H and S, and the chemical formula of the species all be present. \code{checkEOS} and \code{checkGHS} are used by \code{\link{info}} when retrieving species parameters from the database.

  \code{check.obigt} is a function to check self-consistency of each entry (with itself) in the thermodynamic database, using \code{checkEOS} and \code{checkGHS}. The function checks data in both \code{\link{thermo}$obigt} and \code{\link{extdata}/thermo/OBIGT-2.csv}. The output is a table listing only species that exceed at least one of the tolerance limits, giving the name of the database (OBIGT or OBIGT-2), the species index (rownumber in the database), species name and state, and DCp, DV and DG, for the calculated differences (only those above the tolerances are given). This function is used to generate the file found at \code{extdata/thermo/obigt_check.csv}.

  \code{obigt2eos} takes a data frame in the format of \code{thermo$obigt} of one or more rows, removes scaling factors from equations-of-state parameters, and applies new column names depending on the \code{state}. This function is used by both \code{\link{info}} and \code{\link{subcrt}} when retrieving entries from the thermodynamic database. If \code{fixGHS} is TRUE a missing one of G, H or S for any species is calculated from the other two and the chemical formula of the species (used by \code{subcrt} when retrieving database entries).


}

\seealso{
  The default supplementary thermodynamic database for \code{add.obigt} (\code{\link{extdata}/thermo/OBIGT-2.csv}) is used in some of the example calculations in the help page for \code{\link{diagram}} and also in \code{\link{anim.TCA}}.

  \code{\link{mod.buffer}} for updating the list of available buffers.
}

\section{Side Effects}{
  These functions can modify the \code{obigt}, \code{protein} and \code{buffers} dataframes in the global \code{\link{thermo}} list.
}

\examples{
  \dontshow{data(thermo)}
  ## modify/add species
  info(ia <- info("alanine","cr"))
  mod.obigt("alanine",state="cr",G=0,H=0,S=0)
  # now the values of G, H, and S are inconsistent
  # with the elemental composition of alanine
  info(ia)
  # add a species
  mod.obigt("myspecies",formula="CHNOSZ",G=0,H=0)
  info(im <- info("myspecies"))

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
}

\keyword{util}
