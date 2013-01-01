\name{sideeffects}
\alias{sideeffects}
\title{Side effects of functions in CHNOSZ}

\description{
  Some functions in the package access thermodynamic data and system definitions contained in the \code{thermo} data object, as well as modify the object. This help topic should help users understand the major side effects, but does not contain a comprehensive description of these interactions (the code is the ultimate reference).
}

\details{
  When the package is loaded, \code{\link{sys.source}} is used to create a \code{\link{list}} object named \code{\link{thermo}} that is placed in an \code{\link{environment}} named \samp{CHNOSZ:thermo}. Because that environment is \code{\link{attach}ed}, the \code{thermo} object, including all of its components, is accessible in the \code{\link{search}} path. Some functions in CHNOSZ have side effects that modify the contents of \code{thermo}; all such changes can be reverted, and the object restored to its original state, by calling \code{\link{data}(thermo)}.

  In the functions in the package, the greatest number of accessions are to the thermodynamic database (\code{thermo$obigt}), followed by the basis and species definitions (\code{thermo$basis} and \code{thermo$species}). For example, \code{\link{info}} can be used to look up thermodynamic data in \code{thermo$obigt} by the name or chemical formula of a species. As another example, \code{\link{subcrt}} attempts to balance unbalanced chemical reactions with the user-defined basis species in \code{thermo$basis}.

  Some functions modify the thermodynamic database or system definition in \code{thermo}. These are examples of \dQuote{side effects}, since the functions have an effect on the state of the program that persists beyond the lifetime of the objects returned by the functions. In the code, side effects can be recognized by the presence of the \dQuote{superassignment} operator \code{\link{<<-}}. 

  Side effects are not highly desirable in functional programming languages such as \R. The reason this design is adopted in CHNOSZ is that interactive use of \code{\link{basis}} and \code{\link{species}} appeared to the author, in the early stages of developing the package and of learning \R, to be facilitated by not requiring users to assign the results of these functions to objects. Instead, using side effects, the program \dQuote{remembers} the results of these function calls. Experience has shown that this design is usable (especially for new users), and is adaptable to many usage scenarios, but the dependence on side effects probably should be eliminated in the future.

 The two \emph{major} side effects, that most users will encounter, are the basis and species definitions. These functions and a few other modifications (writing) and accessions (reading) of data objects are listed below. The names of objects in this table refer to the components of the \code{thermo} object; for example, one can type \code{thermo$opt} at the command line to access all of the contents of the \code{opt} component, including those not listed in the table.

  \tabular{llll}{
    object \tab writer \tab reader \tab notes \cr
    \code{obigt} \tab \code{\link{mod.obigt}} \tab \code{\link{info}} \tab thermodynamic database \cr
    \code{basis} \tab \code{\link{basis}} \tab \code{\link{species}}, \code{\link{subcrt}} \tab basis definition \cr
    \code{species} \tab \code{\link{species}} \tab \code{\link{affinity}} \tab species definition \cr
    \code{opt$T.units} \tab \code{\link{T.units}} \tab \code{\link{convert}} \tab units \cr
    \code{opt$water} \tab -- \tab \code{\link{water}} \tab formulation for properties of water \cr
    \code{opt$Tr}, \code{Pr} \tab -- \tab \code{\link{GHS}} \tab reference temperature and pressure \cr
    \code{opt$state} \tab -- \tab \code{\link{info}} \tab physical state \cr
    \code{opar} \tab \code{\link{thermo.plot.new}} \tab -- \tab graphical parameters \cr
  }

  Beginning with CHNOSZ version 0.9-9, instances of \code{\link{<<-}} in the code have the effect of \dQuote{superassignment} to an enclosing environment (\samp{CHNOSZ:thermo} on the search path), instead of the \emph{very} highly discouraged assignment to the global environment (user's workspace) used in previous versions of the package. This convention means that if, from the command line, \emph{you} wish to alter something in \code{thermo}, you also should use the \code{\link{<<-}} operator; examples of changing \code{thermo$opt$water} in this manner can found in the help page for \code{\link{water}}.
  
}

\seealso{
  A discussion on \dQuote{What does the "<<-" operator mean?} appeared on the R-help mailing list, \url{https://stat.ethz.ch/pipermail/r-help/2011-April/275905.html}.
}

\examples{
data(thermo)    # side effect: reset the system definition
basis()         # NULL
basis("CHNOS")  # side effect: define the basis species
basis()         # not NULL
data(thermo)    # side effect: reset the system definition
basis()         # NULL
}

\keyword{secondary}
