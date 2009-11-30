\name{revisit}
\alias{revisit}
\alias{richness}
\alias{extremes}
\alias{where.extreme}
\alias{draw.diversity}
\title{Diversity Calculations for Chemical Species}
\description{
  Calculate species richness, or standard deviation, coefficient of variation or Shannon diversity index of activities or logarithms of activities of chemical species, and plot the results. 
}

\usage{
  revisit(logact, which = "cv", as.is = FALSE)
  richness(logact, logactmin = -4, mult = 1, as.is = FALSE)
  extremes(z, what = "minimum")
  where.extreme(z, what= "minimum")
  draw.diversity(d, which = "cv", logactmin = -4, col = par("fg"),
    as.is = FALSE, yline = 2, ylim = NULL, ispecies = NULL,
    add = FALSE, cex = par("cex"), lwd = par("lwd"), mar = par("mar"),
    mult = 1, side = 1:4)
}

\arguments{
  \item{logact}{list, logarithms of activities of species.}
  \item{d}{list, output from \code{\link{diagram}}.}
  \item{which}{character, type of diversity index to calculate.}
  \item{as.is}{logical, do not convert logarithms of activities to activities?}
  \item{logactmin}{numeric, cutoff value for a species to be counted in richness calculations.}
  \item{z}{numeric, matrix of values.}
  \item{what}{character, \samp{minimum} or \samp{maximum}.}
  \item{col}{character, color to use for line plots.}
  \item{yline}{numeric, margin line for y-axis label.}
  \item{ylim}{numeric, limits of y axis.}
  \item{ispecies}{numeric, which species to consider.}
  \item{add}{logical, add to an existing plot?}
  \item{cex}{numeric, character expansion factor.}
  \item{lwd}{numeric, line width.}
  \item{mar}{numeric, plot margin specifications.}
  \item{mult}{numeric, multiply richness counts by this factor.}
  \item{side}{numeric, which sides of plot to draw axes.}
}

\details{


  \code{revisit} and \code{richness} both take as input a list of logarithms of activities of chemical species calculated using \code{\link{diagram}}. The calculations available in \code{revisit} are the Shannon index, standard deviation, and coefficient of variation, which are specified in \code{which} by \samp{shannon}, \samp{sd} and \samp{cv}, respectively. By default the input is assumed to contain logarithms of activities which are exponentiated so that activities of species are used in the calculations, but if \code{as.is} is set to \code{TRUE}, the values are not changed. 


  \code{richness} simply counts the numbers of species whose logarithms of activities are above \code{logactmin}. The species can be weighted differently by supplying values of \code{mult} other than unity. If \code{logactmin} is any character value, its values are calculated as the mean of activities at each point (or mean of logarithms of activities if as.is is TRUE).

  Given a matrix of numeric values in \code{z}, \code{extremes} locates the maximum or minimum values in both dimensions. That is, the \eqn{x}{x} values that are returned are the column numbers where the extreme is found for each row, and the \eqn{y}{y} values that are returned are the row numbers where the extreme is found for each column. \code{where.extreme} takes a matrix and returns the row and column numbers where the overall maximum or minimum value is found. If the extreme value is repeated, the row and column numbers for all instances are returned.

  \code{draw.diversity} takes output from \code{\link{diagram}}, calls either \code{diversity} or \code{richness}, and plots the results. For systems of two variables \code{where.extreme} (and \code{extremes}) are used to indicate the location of the maximum Shannon index or richness (and the ridges around it), or minimum standard deviation or coefficient of variation (and the valleys around it). The ridges and valleys are plotted as dashed lines and are colored green for the \eqn{x}{x} values returned by \code{extremes} and blue for the \eqn{y}{y} values returned by \code{extremes}.

  The name was changed from \code{diversity} to \code{revisit} in \pkg{CHNOSZ}-0.9 because there is a function named \code{\link[vegan]{diversity}} in the \pkg{vegan} package. There is some overlap in the related functions from the different packages. The major difference in the interfaces is that \code{diversity} takes a matrix with species on the columns, while \code{revisit} takes a list with species as the elements of the list. See the examples for a calculation that uses \code{revisit}, as well as \code{\link[vegan]{radfit}} from \pkg{vegan}.

}

\value{

  \code{diversity} and \code{richness} return a numeric (vector or matrix) that has the same dimensions as a single element of \code{logact}. \code{extremes} and \code{where.extreme} both return a list of \eqn{x}{x} and \eqn{y}{y} values. \code{draw.diversity} returnes a list with elements \code{H} (the result of either \code{diversity} or \code{richness}), \code{ixmin} and \code{iymin} (the column and row numbers where the maximum or minimum value is found), \code{xmin} and \code{ymin} (the values of the x and y variables where the maximum or minimum is found) and \code{extval} (the value of \code{H} at the maximum or minimum).

}

\examples{
  \dontshow{data(thermo)}
  ### these examples can be run using longex("revisit")
  \dontrun{
    ## carboxylases from different organisms
    # ribulose phosphate carboxylase/oxygenase
    rubisco <- c("Q6JAI0_9RHOD","A5CKC7_9CHRO","RBL_SYNJA","A1E8R4_9CHLO",
      "A8C9T6_9MYCO","A3EQE1_9BACT","A6YF84_9PROT", "RBL_BRAJA",
      "A1RZJ5_THEPD","RBL_METJA","A3DND9_STAMF","RBL_PYRHO") 
    # acetyl-coenzyme A carboxylase
    accoaco <- c("Q05KD0_HYDTH","Q9F7M8_PRB01","A6CDM2_9PLAN","A0GZU2_9CHLR",
      "ACCA_DEIRA","A1VC70_DESVV","A7WGI1_9AQUI","Q2JSS7_SYNJA",
      "A4AGS7_9ACTN","ACCA_AQUAE","ACCA_CAUCR","A6VIX9_METM7")  
    # calculate affinities as a function of T with
    # buffered logfH2 and CO2.. adjust the glutathione buffer 
    # (more oxidizing than default)
    mod.buffer("GSH-GSSG",c("GSH","GSSG"),logact=c(-3,-7))
    # add a CO2 gas saturation buffer
    mod.buffer("CO2","CO2","gas",-1)
    basis(c("CO2","H2O","NH4+","SO4-2","H2","H+"),
      c("CO2",0,-4,-4,"GSH-GSSG",-7))
    species(c(rubisco,accoaco))
    a <- affinity(T=c(0,160))
    # create speciation diagram with colors
    col <- c(rep("red",12),rep("blue",12))
    d <- diagram(a,residue=TRUE,color=col,ylim=c(-6,-1),legend.x=NULL)
    legend("topleft",col=c("red","blue"),lty=1,
      legend=c("ribulose phosphate carboxylase",
      "acetyl-coenzyme A carboxylase"))
    title(main=paste("Calculated relative abundances of",
      "carboxylases from different organisms",sep="\n"))
    # ... now on to a species richness diagram
    # all the proteins, then rubisco and accoaco
    draw.diversity(d,"richness",logactmin=-3.6)
    draw.diversity(d,"richness",logactmin=-3.6,
      ispecies=1:12,col="red",add=TRUE)
    draw.diversity(d,"richness",logactmin=-3.6,
      ispecies=13:24,col="blue",add=TRUE)
    legend("bottomleft",col=c("red","blue","black"),lty=1,
      legend=c("ribulose phosphate carboxylase",
      "acetyl-coenzyme A carboxylase","all"))
    title(main=paste("Carboxylases with activities",
      "greater than 10^(-3.6)",sep="\n"))

    ## continuing from above ... make a rank-abundance
    ## diagram and fit with a lognormal distribution
    if(require(vegan)) {
      basis("H2",-4)
      a <- affinity()
      logact <- diagram(a,residue=TRUE,do.plot=FALSE)$logact
      act <- 10^as.numeric(logact)
      # we use family=Gamma because our species have activities
      # (i.e., proportional abundances) and not integer counts
      mod <- rad.lognormal(act,family=Gamma)
      plot(mod,main=paste("Relative abundances of carboxylases",
        "fit with lognormal distribution",sep="\n"))
      # calculate Shannon diversity index
      # using revisit (CHNOSZ)
      H1 <- revisit(act,"shannon",as.is=TRUE)
      legend("topright",legend=paste("H =",round(H1,2)),pch=1)
      # using diversity (vegan)
      H2 <- diversity(matrix(act,nrow=1))
      stopifnot(isTRUE(all.equal(H1,H2)))
    }

    ### using grep.file, read.fasta, add.protein
    # calculations for Pelagibacter ubique
    f <- system.file("HTCC1062.faa",package="CHNOSZ")
    # line numbers of all entries in the file
    j <- grep.file(f)  # length = 1354
    # locate entries whose names contain DNA
    j <- grep.file(f,"hypothetical")
    # get the amino acid compositions of these protein
    p <- read.fasta(f,j)
    # add these proteins to CHNOSZ's inventory
    i <- add.protein(p)
    # set up a thermodynamic system
    basis("CHNOS+")
    # calculate affinities in logfO2-pH space
    a <- affinity(H2O=c(-10,-2),O2=c(-82,-76),iprotein=i)
    # calculate the logarithms of activities
    d <- diagram(a,do.plot=FALSE,mam=FALSE)
    # show the protein richness
    draw.diversity(d,"richness",logactmin="")
    mtitle(c("Richness of hypothetical proteins in",
      expression(italic("Pelagibacter ubique"))))
    # show the coefficient of variation of activities
    draw.diversity(d,"cv")
    mtitle(c("Coefficient of variation of hypothetical",
      expression("protein activities in"~italic("P. ubique"))))
  }
}

\keyword{misc}
