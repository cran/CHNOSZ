\name{logact}
\alias{logact}
\alias{logact.mb}
\alias{logact.react}
\title{Calculate Logarithms of Activity of Species}
\description{
Calculate equilibrium chemical activities of species (i.e., equal-affinity frame of reference) taking as input the affinities of formation of the species (in an equal-activity frame of reference). 
}

\usage{
  logact.mb(Astar, nbalance, thisloga)
  logact.react(Astar, av, nbalance, thisloga)
}

\arguments{
  \item{Astar}{numeric, affinities of formation reactions excluding species contribution.}
  \item{nbalance}{numeric, vector of balance coefficients.}
  \item{thisloga}{numeric, logarithm of total activity of balanced thing.}
  \item{av}{numeric, values of affinities of formation reactions.}
}

\details{

  The chemical activities of species shown by \code{\link{diagram}} are calculated using either the \code{logact.mb} or \code{logact.react} functions. The former is used if \code{residue} in \code{diagram} is set to \code{TRUE}, and the latter if \code{residue} is set to \code{FALSE}. 

  The input values are in a list, \code{Astar}, all element of the list having the same dimensions, and each element contains the chemical affinities of the formation reactions of one of the species of interest at unit activity. The metastable equilibrium activities calculated using either function satisfy the constraints that 1) the resulting chemical affinities of the formation reactions of the species are all equal and 2) the total activity of the conserved component (\code{thisloga}) is unchanged.

  In \code{logact.react} (the algorithm described in Dick, 2008 and the only one available prior to CHNOSZ-0.8), the calculations of relative abundances of species use the activities in the \code{\link{affinity}} output as initial guesses, and attempt to solve a system of equations that represent the two constraints stated above. So, if you supply a value for \code{thisloga} that is much different from that of the initial guess you may end up with errors from \code{\link{uniroot}} such as "f() values at end points not of opposite sign". 

  In \code{logact.mb} (available beginning with CHNOSZ-0.8), the chemical activities of species are calculated using the Maxwell-Boltzmann distribution, which is much faster than the equation-solving approach used above. However, this algorithm is limited to systems where the balance coefficients are all unity. Therefore, this function is only called by \code{diagram} for systems of proteins when \code{residue} in that function is set to TRUE.

}

\value{
  A numeric value (possibly an array) with dimensions equal to those of one of the list elements of \code{Astar}.
}

\keyword{misc}
