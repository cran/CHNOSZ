\name{util.data}
\alias{util.data}
\alias{add.obigt}
\alias{mod.obigt}
\alias{add.protein}
\alias{change}
\title{Modify the Thermodynamic Database}
\description{
  These functions add species to or alter properties of species in the thermodynamic database (\code{thermo$obigt}) or in the buffer definition table (\code{thermo$buffer}).
}

\usage{
  add.obigt(file = system.file("extdata/OBIGT-2.csv", 
    package = "CHNOSZ"), force = FALSE)
  mod.obigt(species, ..., missingvalues = NA)
  add.protein(file="protein.csv")
  change(name, ...)
}

\arguments{
  \item{file}{character, path to a file.}
  \item{force}{logical, force replacement of already existing species?}
  \item{...}{character or numeric, properties of species to modify in the thermodynamic database.}
  \item{species}{character, formulas or names of species to modify or add to the thermodynamic database.}
  \item{missingvalues}{numeric, values to assign to undefined properties.}
  \item{name}{character or numeric, name (or numeric index) of species or name of buffer to be modified.}
}

\details{

  \code{add.obigt} and \code{add.protein} read data from the specified \code{file} and add it to either \code{\link{thermo}$protein} or \code{\link{thermo}$obigt}, as appropriate. Set \code{force} to TRUE to replace species that exist in the thermodynamic database (each unique combination of a name and a state in the database is considered to be a unique species). \code{force}, if not specified, reverts to TRUE if the \code{file} is left at its default. The default supplementary thermodynamic database (\code{extdata/OBIGT-2.csv}) is used in the example calculations in the help page for \code{\link{diagram}} and also in \code{\link{anim.TCA}}.

  \code{mod.obigt} changes one or more of the properties of one or more species or adds species to the thermodynamic database. These changes are lost if you reload the database by calling \code{\link{data}(thermo)} or if you quit the \R session without saving it. To modify the properties of species, give the names in the \code{species} argument and supply other arguments: if one of these arguments is \code{state}, species in those states will be updated. Additional arguments refer to the name of the property(s) to be updated and correspond to the column names of \code{\link{thermo}$obigt} (the names of the properties are matched to any part of compound column names, such as \samp{z.T}). The values provided should be in the units specifed in the documentation for the \code{thermo} data object. To add species, supply the new names in \code{species} and provide an argument named \code{formula} with the corresponding chemical formulas. Additional arguments refer to any of the properties you wish to specify. Properties that are not specified are assigned the value of \code{missingvalues} which is \code{NA} by default (however if \code{state} is missing it is set to the value of \code{thermo$opt$state}). The values returned (\code{\link{invisible}}-y) by \code{mod.obigt} are the rownumbers of the affected species.

  \code{change} is a wrapper function to \code{\link{mod.obigt}} and \code{\link{mod.buffer}}. The \code{name} provided in the argument refers to the name or numeric index of the species to update or add using \code{\link{mod.obigt}}, unless the \code{name} begins with an underscore character, in which case the remaining part of the name (after the underscore) is passed to \code{\link{mod.buffer}}. The arguments in \code{...} are sent without change to the subordinate function.


}

\seealso{
  \code{\link{mod.buffer}} for updated the list of available buffers.
}

\section{Side Effects}{
  These functions can modify the \code{obigt}, \code{protein} and \code{buffer} dataframes in the global \code{\link{thermo}} list.
}

\examples{
  \dontshow{data(thermo)}
  ## modify/add species
  info(t <- info("alanine","cr"))
  mod.obigt("alanine",state="cr",G=0,H=0,S=0)
  # now the values of G, H, and S are inconsistent
  # with the elemental composition of alanine
  info(t)
  # add a species
  mod.obigt("myspecies",formula="CHNOSZ",G=0,H=0)
  info(t <- info("myspecies"))
}

\keyword{misc}
