\name{revisit}
\alias{revisit}
\alias{richness}
\alias{extremes}
\alias{where.extreme}
\title{Diversity Calculations for Chemical Species}
\description{
  Calculate species richness, or standard deviation, coefficient of variation or Shannon diversity index of activities or logarithms of activities of chemical species, and plot the results. 
}

\usage{
  revisit(d, target = "cv", loga.ref = NULL,
    do.plot = NULL, col = par("fg"), yline = 2, ylim = NULL, 
    ispecies = NULL, add = FALSE, cex = par("cex"), lwd = par("lwd"), 
    mar = NULL, side = 1:4, xlim = NULL, labcex = 0.6, pch = 1, 
    legend = "", legend.x = NULL, lpch = NULL, main = NULL, 
    lograt.ref = NULL, plot.ext = TRUE)
  extremes(z, target)
  where.extreme(z, target, do.sat = FALSE)
}

\arguments{
  \item{d}{list, output from \code{\link{diagram}}, or list of logarithms of activities of species.}
  \item{target}{character, what statistic to calculate.}
  \item{loga.ref}{numeric, logarithm of activities for comparison statistics}
  \item{do.plot}{logical, make a plot?}
  \item{col}{character, color to use for points or liness.}
  \item{yline}{numeric, margin line for y-axis label.}
  \item{ylim}{numeric, limits of y axis.}
  \item{ispecies}{numeric, which species to consider.}
  \item{add}{logical, add to an existing plot?}
  \item{cex}{numeric, character expansion factor.}
  \item{lwd}{numeric, line width.}
  \item{mar}{numeric, plot margin specifications.}
  \item{side}{numeric, which sides of plot to draw axes.}
  \item{xlim}{numeric, limits of x axis.}
  \item{labcex}{numeric, character expansion factor for species labels.}
  \item{pch}{numeric, plotting symbol(s) to use for points.}
  \item{legend}{character, text to use for legend.}
  \item{legend.x}{character, placement of legend.}
  \item{lpch}{numeric, plotting symbol(s) to use in legend.}
  \item{main}{character, main title for plot.}
  \item{lograt.ref}{numeric, log10 of reference abundance ratios.}
  \item{plot.ext}{logical, show the location of the extreme value(s)?}
  \item{z}{numeric, matrix of values.}
  \item{do.sat}{logical, identify multiple extreme values.}
}

\details{

  The purpose of \code{richness} is to calculate and visualize summary statistics for logarithms of activities of chemical species. For most uses, supply the output of \code{\link{diagram}} as the value for \code{d}. Alternatively, \code{d} can be a list of logarithms of activities; the list elements each correspond to a different species and can be vectors, matrices, or higher-dimensional arrays, but they must all have the same dimensions. (This is always the case for \code{d$logact} if \code{d} is the output from \code{diagram}; the dimensionality is determined by the number of variables used in the calculations of \code{\link{affinity}}.) The type of statistic to be calculations is indicated by \code{target}, as summarized in the following table.

  \tabular{llll}{
    target \tab description \tab extremum \tab additional arguments \cr
    \code{sd} \tab standard deviation \tab min \tab none \cr
    \code{cv} \tab coefficient of variation \tab min \tab none \cr
    \code{shannon} \tab Shannon diversity index \tab max \tab none \cr
    \code{qqr} \tab correlation coefficient on q-q plot (normal distribution) \tab max \tab none \cr
    \code{richness} \tab species richness \tab max \tab loga.ref \cr
    \code{cvrmsd} \tab coefficient of variation of RMSD \tab min \tab loga.ref \cr
    \code{spearman} \tab Spearman correlation coefficient \tab max \tab loga.ref \cr
    \code{pearson} \tab Pearson correlation coefficient \tab max \tab loga.ref \cr
  }

  \code{sd}, \code{cv}, \code{shannon} and \code{qqrr} all operate on just the sample values. \code{richness} counts the numbers of species whose logarithms of activities are above \code{log.min}. \code{cvrmsd}, \code{spearman} and \code{pearson} are comparison statistics where \code{loga.target} represents the observed values. \code{ratio} determines the correlation coefficient of a predicted change in loga ratios (\code{d$logact} vs. loga.ref) plotted agains observed changed in loga ratios (e.g., from changes in protein expression deduced from microarray experiments; given in loga.target)

  If \code{do.plot} is TRUE, \code{d} is the output from \code{diagram}, and the number of variables is 1 or 2, the results are plotted -- a line diagram in 1 dimension or a contour plot in 2 dimensions. 

  The value of \samp{extremum} in the table shows whether the extreme value that optimizes the system is the minimum (\code{sd}, \code{cv}, \code{cvrmsd}) or the maximum (all the others). On plots the location of the extreme value is indicated (by a dashed vertical line on a 1-D plot or a point marked by an asterisk on a 2-D plot). On 2-D plots the valleys (or ridges) leading to the location of the extremum are plotted. The ridges or valleys are plotted as dashed lines and are colored green for the \eqn{x}{x} values returned by \code{extremes} and blue for the \eqn{y}{y} values returned by \code{extremes}.

  The location of the extreme value in a matrix or vector \code{z} is calculated using \code{where.extreme}. Whether the extreme is the minimum or the maximum value depends on the value of \code{target}. For matrices, if \code{do.sat} is TRUE, if the extreme value is repeated, the row and columns numbers for all instances are returned. Given a matrix of numeric values in \code{z}, \code{extremes} locates the maximum or minimum values in both dimensions. That is, the \eqn{x}{x} values that are returned are the column numbers where the extreme is found for each row, and the \eqn{y}{y} values that are returned are the row numbers where the extreme is found for each column. 

  If \code{lograt.ref} is provided, these values are the reference values for logarithm of abundance ratio.
 
  The function name was changed from \code{diversity} to \code{revisit} in CHNOSZ-0.9 because there is a function named \code{\link[vegan]{diversity}} in the \pkg{vegan} package. Note that while \code{diversity} takes a matrix with species on the columns, \code{revisit} takes a list with species as the elements of the list. 

}

\value{

  \code{revisit} returns a list containing at least an element named \samp{H} giving the calculated values for the target statistic. This has the same dimensions as a single element of \code{d} (or \code{d$logact}, if \code{d} was the output from \code{diagram}). For calculations as a function of one or two variables, the output also contains the elements \code{ix} (location of the extremum in the first direction), \code{x} (x-value at the extremum), and \code{extval} (extreme value). For calculations as a function of two variables, the output also contains the elements \code{iy} (location of the extremum in the second direction) and \code{y} (y-value at the extremum).

}

\examples{
  \dontshow{data(thermo)}
  \donttest{
    ### using grep.file, read.fasta, add.protein
    # calculations for Pelagibacter ubique
    f <- system.file("extdata/fasta/HTCC1062.faa.xz",package="CHNOSZ")
    # what proteins to select (set to "" for all proteins)
    w <- "ribosomal"
    # locate entries whose names contain w
    j <- grep.file(f,w)
    # get the amino acid compositions of these protein
    p <- read.fasta(f,j)
    # add these proteins to CHNOSZ's inventory
    i <- add.protein(p)
    # set up a the chemical system
    basis("CHNOS+")
    # calculate affinities of formation in logfO2 space
    a <- affinity(O2=c(-90,-60),iprotein=i)
    # show the equilibrium activities
    d <- diagram(a,cex=1.5,logact=0)
    # make a title
    expr <- as.expression(substitute(x~y~"proteins in"~
      italic("P. ubique"),list(x=length(j),y=w)))
    mtitle(c("Equilibrium activities of",expr),cex=1.5)
    # show the coefficient of variation
    revisit(d,"CV",cex=1.5)
    mtitle(c("CV of equilibrium activities of",expr),cex=1.5)
    # calculate affinities in logfO2-logaH2O space
    a <- affinity(O2=c(-90,-60),H2O=c(-20,0),iprotein=i)
    # calculate the equilibrium activities
    d <- diagram(a,do.plot=FALSE,mam=FALSE,logact=0)
    # show the coefficient of variation
    revisit(d,"CV",cex=1.5)
    mtitle(c("CV of equilibrium activities of",expr),cex=1.5)
  }
}

\keyword{primary}
