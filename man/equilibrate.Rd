\name{equilibrate}
\alias{equilibrate}
\alias{balance}
\alias{equil.boltzmann}
\alias{equil.reaction}
\title{Equilibrium Chemical Activities of Species}
\description{
Calculate equilibrium chemical activities of species from the affinities of formation of the species at unit activity. 
}

\usage{
  equilibrate(aout, balance=NULL, loga.balance=NULL,
    ispecies=1:length(aout$values), normalize=FALSE, stay.normal=FALSE)
  balance(aout, balance=NULL)
  equil.boltzmann(Astar, n.balance, loga.balance)
  equil.reaction(Astar, n.balance, loga.balance)
}

\arguments{
  \item{aout}{list, output from \code{\link{affinity}}}
  \item{balance}{character or numeric, how to balance the transformations}
  \item{ispecies}{numeric, which species to include}
  \item{normalize}{logical, normalize the molar formulas of species by the balancing coefficients?}
  \item{stay.normal}{logical, report results for the normalized formulas?}
  \item{Astar}{numeric, affinities of formation reactions excluding species contribution}
  \item{n.balance}{numeric, number of moles of conserved component in the formation reactions of the species of interest}
  \item{loga.balance}{numeric, logarithm of total activity of balanced quantity}
}

\details{
\code{equilibrate} provides an interface to calculate the chemical activities of species in metastable equilibrium, in an open system at constant temperature and pressure and chemical activities of basis species, and with linear balancing constraints on transformations.

It takes as input \code{aout}, the output from \code{\link{affinity}}, which may be calculated from a multidimensional grid of conditions.
The equilibrium chemical activities of species are calculated using either the \code{equil.reaction} or \code{equil.boltzmann} functions, the latter only if the balance is on one mole of species.

As \code{aout} contains the chemical affinities of formation reactions of each species of interest, \code{equilibrate} in order to function needs to be provided constraints on how to balance the reactions representing transformations between the species.
\code{balance} returns the balancing coefficients, where \code{balance} indicates the balancing constraints, according to the following scheme:
\tabular{ll}{
      \samp{NULL} \tab autoselect using \code{\link{which.balance}}\cr
      \emph{name of basis species} \tab balance on this basis species\cr
      \samp{length} \tab balance on length of proteins\cr
      \samp{1} \tab balance on one mole of species\cr
      \emph{numeric vector} \tab user-defined constraints\cr
}

The default value of NULL for \code{balance} indicates to select the first shared basis species in all formation reactions identified using \code{\link{which.balance}}, or if that fails, to set the balance to \samp{1}.
However, if all the species (as listed in code \code{aout$species}) are proteins (have an underscore character in their names), the default value of NULL for \code{balance} indicates to use \samp{length} as the balance.

NOTE: the summation of activities assumes an ideal system, so \sQuote{molality} is implied by \sQuote{activity} in the following.
\code{loga.balance} gives the logarithm of the total activity of \code{balance} (which is total activity of species for \samp{1} or total activity of amino acid residue-equivalents for \samp{length}).
If \code{loga.balance} is missing, its value is taken from the activities of species listed in \code{aout}; this default is usually the desired operation.

\code{normalize} if TRUE indicates to normalize the molar formulas of species by the balance coefficients.
This operation is intended for systems of polymers, such as proteins, whose conventional formulas are much larger than the basis speices.
The normalization also applies to the balancing coefficients, which as a result consist of \samp{1}s.
\code{normalize} has the same effect as did \code{diagram(..., residue=TRUE)} in versions of CHNOSZ before 0.9-9.
After normalization and equilibration, the equilibrium activities are then un-normalized (for the original formulas of the species), unless \code{stay.normal} is TRUE.

\code{ispecies} can be supplied to identify a subset of the species to include in the calculation.

\code{equil.boltzmann} is used to calculation the equilibrium activities if \code{balance} is \samp{1} (including the normalized result when \code{normalize} is TRUE), otherwise \code{equil.reaction} is called.

}

\section{Algorithms}{

The input values to \code{equil.reaction} and \code{equil.boltzmann} are in a list, \code{Astar}, all elements of the list having the same dimensions; they can be vectors, matrices, or higher-dimensionsal arrays.
Each list element contains the chemical affinities of the formation reactions of one of the species of interest (in dimensionless base-10 units, i.e. A/2.303RT), calculated at unit activity of the species of interest.
The equilibrium activities (in base-10 log units) of the species of interest returned by either function satisfy the constraints that 1) the final chemical affinities of the formation reactions of the species are all equal and 2) the total activity of the conserved component is equal to (\code{loga.balance}). 
The first constraint does \emph{not} impose a complete equilibrium, where the affinities of the formation reactions are all equal to zero, but allows for a metastable equilibrium, where the affinities of the formation reactions are equal to each other.

In \code{equil.reaction} (the algorithm described in Dick, 2008 and the only one available prior to CHNOSZ-0.8), the calculations of relative abundances of species are based on a solving a system of equations representing the two constraints stated above.
A close approximation of the solution is found using \code{\link{uniroot}}.
Prior to CHNOSZ_0.9-9, the values in the \code{Astar} were used to generate initial guesses of the logarithms of activities of species; values of \code{loga.balance} that were hugely different from these guesses could generate errors from \code{uniroot} such as "f() values at end points not of opposite sign".
Now (from version 0.9-9), a more flexible method for generating guesses is in place. 

In \code{equil.boltzmann} (algorithm available beginning with CHNOSZ-0.8), the chemical activities of species are calculated using the Boltzmann distribution.
This calculation is faster than the algorithm of \code{equil.reaction}, but is limited to systems where the transformations are all balanced on one mole of species.

}

\seealso{
\code{\link{diagram}} has examples of using \code{equilibrate} to make equilibrium activity diagrams. \code{\link{revisit}} can be used to perform further analysis of the equilibrium activities.
\code{\link{palply}} is used by both \code{equil.reaction} and \code{equil.boltzmann} to parallelize intensive parts of the calculations if \pkg{parallel} is loaded.
}

\value{
\code{equil.reaction} and \code{equil.boltzman} each return a list with dimensions and length equal to those of \code{Astar}, giving the \code{\link{log10}} of the equilibrium activities of the species of interest.
\code{equilibrate} returns a list, containing first the values in \code{aout}, to which are appended \code{m.balance} (the balancing coefficients if \code{normalize} is TRUE, a vector of \samp{1}s otherwise), \code{n.balance} (the balancing coefficients if \code{normalize} is FALSE, a vector of \samp{1}s otherwise) and \code{loga.equil} (the calculated equilibrium activities of the species).
\code{balance} returns a list containing the balancing coefficients (\code{n}) and a textual description (\code{description}).
}

\examples{
\dontshow{data(thermo)}
## equilibrium in a simple system:
## ionization of carbonic acid
basis("CHNOS+")
species(c("CO2", "HCO3-", "CO3-2"))
# set unit activity of the species (0 = log10(1))
species(1:3, 0)
# calculate Astar (for unit activity)
res <- 100
Astar <- affinity(pH=c(0, 14, res))$values
# the logarithms of activity for a total activity
# of the balanced quantity (C or CO2) equal to 0.001
loga.boltz <- equil.boltzmann(Astar, c(1, 1, 1), 0.001)
# calculated another way
loga.react <- equil.reaction(Astar, c(1, 1, 1), 0.001)
# probably close enough for most purposes
stopifnot(all.equal(loga.boltz, loga.react))
# the first ionization constant (pKa)
ipKa <- which.min(abs(loga.boltz[[1]] - loga.boltz[[2]]))
pKa.equil <- seq(0, 14, length.out=res)[ipKa]
# calculate logK directly
logK <- subcrt(c("CO2","H2O","HCO3-","H+"), c(-1, -1, 1, 1), T=25)$out$logK
# we could decrease tolerance here by increasing res
stopifnot(all.equal(pKa.equil, -logK, tol=1e-2))
}

\references{
  Dick, J. M. (2008) Calculation of the relative metastabilities of proteins using the CHNOSZ software package. \emph{Geochem. Trans.} \bold{9}:10. \url{http://dx.doi.org/10.1186/1467-4866-9-10}
}

\keyword{secondary}
