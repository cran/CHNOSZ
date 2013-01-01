\name{revisit}
\alias{revisit}
\alias{richness}
\alias{extremes}
\alias{optimal.index}
\title{Plots and Optima of Objective Functions}
\description{
Calculate values of an objective function from logarithms of activities of chemical species and (for some objectives) reference logarithms of activity.
Make line or contour plots showing the values of the objective function and the positions of the optima (minimum or maximum).
}

\usage{
  revisit(eout, objective = "CV", loga2 = NULL, ispecies = NULL,
    col = par("fg"), yline = 2, ylim = NULL, cex = par("cex"),
    lwd = par("lwd"), mar = NULL, side = 1:4, xlim = NULL, labcex = 0.6,
    pch = 1, main = NULL, plot.it = NULL, add = FALSE, plot.optval = TRUE,
    style.2D = "contour")
  optimal.index(z, objective)
  extremes(z, objective)
}

\arguments{
  \item{eout}{list, output from \code{\link{equilibrate}}, containing logarithms of activities of species}
  \item{objective}{character, name of \code{\link{objective}} function}
  \item{loga2}{numeric array, reference values of logarithm of activities}
  \item{ispecies}{numeric, which species to consider}
  \item{col}{character, color to use for points or lines}
  \item{yline}{numeric, margin line for y-axis label}
  \item{ylim}{numeric, limits of y axis}
  \item{cex}{numeric, character expansion factor}
  \item{lwd}{numeric, line width}
  \item{mar}{numeric, plot margin specifications}
  \item{side}{numeric, which sides of plot to draw axes}
  \item{xlim}{numeric, limits of x axis}
  \item{labcex}{numeric, character expansion factor for species labels}
  \item{pch}{numeric, plotting symbol(s) to use for points}
  \item{main}{character, main title for plot}
  \item{plot.it}{logical, make a plot?}
  \item{add}{logical, add to an existing plot?}
  \item{plot.optval}{logical, show the location of the optimal value(s)?}
  \item{style.2D}{character, type of 2-D plot}
  \item{z}{numeric, matrix of values (\samp{H} element of output from \code{revisit})}
}

\details{

Most commonly, the output of \code{\link{equilibrate}} is used as the value for \code{eout}.
The type of calculation is indicated by \code{objective}, giving the name of an \code{\link{objective}} function.
The list of logarithms of chemical activities in \code{eout$loga.equil} is passed as \code{loga1} to the objective function.
If the objective function has an argument \code{a1} instead of \code{loga1}, the activities instead of their logarithms are passed.
Reference logarithms of activities, for those objective functions that require them, are supplied in the \code{loga2} argument.
Generally, \code{loga2} must be a numeric vector with length equal to that of \code{loga1} (i.e., number of species).
If a single numeric value is supplied for \code{loga2}, it is recycled to the length of \code{loga1}.

If \code{plot.it} is TRUE, and \code{eout} is the output from \code{equilibrate}, and the number of variables is 1 or 2, the results are plotted --- a line diagram in 1 dimension or a contour plot in 2 dimensions.
No plot is made for calculations at a single condition (0-D, no variation), except for the \code{qqr} objective, showing a quantile-quantile plot (\code{\link{qqnorm}}).
\code{style.2D} can be set to \code{\link{image}} to fill the plot with colors instead of the \code{\link{contour}} plot that is the default.

If \code{plot.optval} is TRUE, the location of the optimum (or optima) is indicated by a dashed vertical line(s) on a 1-D plot or a point(s) marked by an asterisk on a 2-D plot.
Also, on 2-D plots, the locations of the optima at each grid line perpendicular to the \eqn{x}{x} and \eqn{y}{y} axes are plotted.
These points follow major ridges or valleys, and are plotted as dashed lines colored green for the \eqn{x}{x} and blue for the \eqn{y}{y} values returned by \code{extremes}.

An alternative source for the \code{eout} argument is any list of numeric values, each element of which corresponds to a different observation (such as a single species), all having the same dimensions (as vectors, matrices or higher-dimensional arrays)
In this case, plotting is disabled, since the names of the variables are not in the input.

\samp{revisit} is a loose anagram of \samp{diversity}, which was the provisional name of the function but was changed in CHNOSZ-0.9.
While the \code{\link[vegan]{diversity}} function (in \pkg{vegan}) operates on a matrix with (biological) species on the columns, \code{revisit} operates on a list with (chemical) species as the elements of the list.
The name of the \samp{H} output value is the conventional symbol for the Shannon diversity index, which was the first target statistic to be implemented in \code{revisit}.

}

\value{

\code{revisit} returns a list containing at least an element named \samp{H}, giving the calculated values of the objective function.
For 1 or 2 dimensions of varibility of chemical conditions, the output also contains the elements \code{ixopt} and \code{iyopt} (1-D and 2-D) and \code{iyopt} and \code{yopt} (2-D) indicating the positions and values of the optimum.
For calculations in more than two dimensions, the output contains \code{iopt} which is a matrix 

\code{optimum.index} returns the index (or indices) of the optimal values in \code{z} using \code{\link{array}(..., arr.ind=TRUE)}.
\code{extremes} returns \eqn{x}{x} values that are the column numbers where the optimum is found for each row, and \eqn{y}{y} values that are the row numbers where the optimum is found for each column.
The \samp{optimum} attribute of the \code{objective} function indicates whether minimal or maximal values are used.

}

\seealso{

\code{\link{findit}} for gridded search of chemical activities, temperature and/or pressure that optimize the objective function.

}

\examples{
## example of defining a new objective function
# count the species with logarithms of activity greater than loga2
count <- function(loga1, loga2) rowSums(loga1 > loga2)
# set the attribute indicating the type of optimum
attr(count, "optimum") <- "maximal"
# equilibrate a system of amino acids
basis("CHNOS")
species(aminoacids(""))
a <- affinity(O2=c(-80, -60))
e <- equilibrate(a)
# make a plot
r <- revisit(e, "count", -5)
title(main="amino acids with metastable log activities > -5")

\donttest{# can also make a 2-D plot
a <- affinity(O2=c(-74, -60), H2O=c(-3, 3))
e <- equilibrate(a)
r <- revisit(e, "count", -5, style.2D="image", plot.optval=FALSE)
title(main="amino acids with metastable log activities > -5")

## 'revisit' calculations for amino acids
opar <- par(mfrow=c(2, 2))
basis("CHNOS+")
species(aminoacids(""))
# chemical affinities as a function of logarithm of oxygen fugacity
a <- affinity(O2=c(-85, -60))
# shows the equilibrium abundances of the amino acids
e <- equilibrate(a)
diagram(e)
mtitle(c("20 amino acids", "balanced on CO2"))
# show a legend with input constraints
db <- describe.basis(ibasis=3)
dp <- describe.property("T", 25)
legend("bottomright", c(dp, db))
# default is to plot coefficient of variation
r <- revisit(e)
# show a title with the optimal conditions
mincv <- format(r$optimum, digits=3)
t1 <- paste("minimum coeff of variation,", mincv, "at:")
# the logfO2 that minimized the C.V.
basis("O2", r$x)
t2 <- describe.basis(ibasis=5)
mtitle(c(t1, t2))
# chemical affinities as a function of two other variables
a <- affinity(NH3=c(-10, 10, 40), T=c(0, 80, 40))
diagram(a, fill="heat")
# show a legend with input constraints
db <- describe.basis(ibasis=5)
legend("bottomright", db)
# contour plot of the CV
e <- equilibrate(a)
r <- revisit(e)
# show a title with the optimal conditions
mincv <- format(r$optimum, digits=3)
t1 <- paste("minimum coeff of variation,", mincv, "at:")
# the logaNH3 and T that minimized the C.V.
basis("NH3", r$x)
db <- describe.basis(ibasis=3)
dp <- describe.property("T", r$y)
t2 <- substitute(list(dp, db), list(dp=dp[[1]], db=db[[1]]))
mtitle(c(t1, as.expression(t2)))
par(opar)

## calculations for proteins in Pelagibacter ubique
## using grep.file, read.fasta, add.protein
f <- system.file("extdata/fasta/HTCC1062.faa.xz", package="CHNOSZ")
# what proteins to select (set to "" for all proteins)
w <- "ribosomal"
# locate entries whose names contain w
j <- grep.file(f, w)
# get the amino acid compositions of these protein
aa <- read.fasta(f, j)
# add these proteins to CHNOSZ's inventory
ip <- add.protein(aa)
# set up a the chemical system
basis("CHNOS+")
# calculate affinities of formation in logfO2 space
a <- affinity(O2=c(-85, -60), iprotein=ip)
# show the equilibrium activities
opar <- par(mfrow=c(2, 2))
e <- equilibrate(a, loga.balance=0, normalize=TRUE)
diagram(e, names=NULL)
# make a title
expr <- as.expression(substitute(x~y~"proteins in"~
  italic("P. ubique"), list(x=length(j), y=w)))
mtitle(c("Equilibrium activities of", expr))
# show the coefficient of variation
revisit(e, "CV")
mtitle(c("CV of equilibrium activities of", expr))
# calculate affinities in logfO2-logaH2O space
a <- affinity(O2=c(-85,-65), H2O=c(-10,0), iprotein=ip)
# show the predominances
diagram(a, normalize=TRUE, fill="heat")
# calculate the equilibrium activities
e <- equilibrate(a, loga.balance=0, normalize=TRUE)
# show the coefficient of variation
r <- revisit(e, "CV")
stopifnot(r$ix==37)
stopifnot(r$iy==53)
mtitle(c("CV of equilibrium activities of", expr))
par(opar)
}
}

\keyword{primary}
