\name{findit}
\alias{findit}
\alias{plot.findit}
\title{Multidimensionally Optimize Chemical Activities}
\description{
  Use a gridded search to find a combination of one or more of chemical activities of basis species, temperature and/or pressure that maximize or minimize a target function of the metastable equilibrium chemical activities of the species of interest.
}

\usage{
  findit (lims = list(), target = "cv", n = NULL, iprotein = NULL, 
    do.plot = TRUE, T = 25, P = "Psat", res = NULL, labcex = 0.6, 
    loga.ref = NULL, as.residue = FALSE, loga.tot = 0)
  \method{plot}{findit}(x, which=NULL, mar=c(3.5,5,2,2), xlab="iteration", \dots)
}

\arguments{
  \item{lims}{list, specification of search limits.}
  \item{target}{character, target statistic to optimize.}
  \item{n}{numeric, number of iterations.}
  \item{res}{numeric, grid resolution (number of points on one edge).}
  \item{iprotein}{numeric, indices of proteins.}
  \item{do.plot}{logical, make a plot?}
  \item{T}{numeric, temperature.}
  \item{P}{numeric, pressure; or character, "Psat".}
  \item{labcex}{numeric, character expansion for plot labels.}
  \item{loga.ref}{numeric, reference logarithms of activity of species.}
  \item{as.residue}{logical, use the activities of residues instead of proteins?}
  \item{loga.tot}{numeric, logarithm of total activity of residues (or other immobile component).}
  \item{x}{list, object of class \code{findit}.}
  \item{which}{numeric, which of the parameters to plot.}
  \item{mar}{numeric, plot margin specification.}
  \item{xlab}{character, x-axis label.}
  \item{\dots}{additional arguments passed to \code{plot}.}
}

\details{
  \code{findit} implements a gridded optimization to find the minimum or maximum value of \code{target} as a function of one or more of the chemical activities, temperature and/or pressure whose ranges are listed in \code{lims}. Any target available in \code{\link{revisit}} can be optimized. Generally, the system (\code{\link{basis}} species and \code{\link{species}} of interest) must be set up before calling this function. If \code{iprotein} is supplied, indicating a set of proteins to use in the calculation, the definition of the \code{species} is not required. 

  \code{lims} is a list, each element of which contains two values indicating the range of the parameter indicated by the \code{\link{names}} of \code{lims}. The \code{names} should be formula(s) of one or more of the basis speices, \samp{T} and/or \samp{P}. If the latter two are missing, the calculations are performed at isothermal and/or isobaric conditions indicated by \code{T} and \code{P}. 

  It can be used for one, two, or more parameters. If \eqn{nd} is the number of parameters (dimensions), default values of \code{n} and \code{res} come from the following table. These settings are selected mostly for quickness of testing the function in high-dimensional space. Detailed studies of a system might have to use more iterations and/or higher resolutions.

  \tabular{rrrrr}{
    \code{nd} \tab \code{n} \tab \code{res} \tab \code{res^nd} \tab \code{rat} \cr
    1 \tab 4 \tab 128 \tab 128 \tab 2/(1+sqrt(5)) \cr
    2 \tab 6 \tab 64 \tab 4096 \tab 2/(1+sqrt(5)) \cr
    3 \tab 6 \tab 16 \tab 4096 \tab 0.75 \cr
    4 \tab 8 \tab 8 \tab  4096 \tab 0.8 \cr
    5 \tab 12 \tab 6 \tab 7776 \tab 0.8 \cr
    6 \tab 12 \tab 4 \tab 4096 \tab 0.85 \cr
    7 \tab 12 \tab 4 \tab 16384 \tab 0.85 \cr
  }

  The function performs \code{n} iterations. At first, the limits of the parameters given in \code{lims} define the extent of a \eqn{nd}-dimensional box around the space of interest. The value of \code{target} is calculated at each of the \eqn{res^{nd}}{res^nd} grid points and and optimum value located (see \code{\link{revisit}} and \code{\link{where.extreme}}). In the next iteration the new search box is centered on the location of the optimum value, and the edges are shrunk so their length is \code{rat} * the length in the previous step. If the limits of any of the parameters extend beyond those in \code{lims}, they are pushed in to fit (preserving the difference between them).

  \code{plot.findit} plots the values of the parameters and the target statistic as a function of the number of iterations.

}

\value{
  \code{findit} returns a list having class \code{findit} with elements \code{value} (values of the parameters, and value of the target statistic, at each iteration), \code{lolim} (lower limits of the parameters) and \code{hilim} (upper limits of the parameters).
}

\examples{
  \dontshow{data(thermo)}
  \donttest{
    # an inorganic example: sulfur species
    basis("CHNOS+")
    basis("pH",5)
    species(c("H2S","S2-2","S3-2","S2O3-2","S2O4-2","S3O6-2",
      "S5O6-2","S2O6-2","HSO3-","SO2","HSO4-"))
    # to minimize the standard deviations of the loga of the species
    target <- "sd"
    # this one gives logfO2=-27.8
    f1 <- findit(list(O2=c(-50,-15)),target,T=325,P=350,n=3)
    # this one gives logfO2=-30.0 pH=9.38
    f2 <- findit(list(O2=c(-50,-15),pH=c(0,14)),target,T=325,P=350,res=16,n=4)

    # an organic example: 
    # find chemical activities where metastable activities of
    # selected proteins in P. ubique have high correlation
    # with a lognormal distribution (i.e., maximize r of q-q plot)
    f <- system.file("extdata/HTCC1062.faa",package="CHNOSZ")
    # search for three groups of proteins
    myg <- c("ribosomal","nucle","membrane")
    g <- lapply(myg,function(x) grep.file(f,x))
    # note that some proteins match more than one search term
    uug <- unique(unlist(g))
    # read their amino acid compositions from the file
    p <- read.fasta(f,uug)
    # add these proteins to thermo$protein
    ip <- add.protein(p)
    # load a predefined set of uncharged basis species
    # (speeds things up as we won't model protein ionization)
    basis("CHNOS")
    # make colors for the diagram
    rgbargs <- lapply(1:3,function(x) as.numeric(uug \%in\% g[[x]]))
    col <- do.call(rgb,c(rgbargs,list(alpha=0.5)))
    # get point symbols (use 1,2,4 and their sums)
    pch <- colSums(t(list2array(rgbargs)) * c(1,2,4))
    # plot 1: calculated logarithms of chemical activity
    # as a function of logfO2 ... a bundle of curves near logfO2 = -77
    a <- affinity(O2=c(-90,-50),iprotein=ip)
    d <- diagram(a,logact=0,col=col)
    # plot 2: q-q correlation coefficient
    # it shows lognormal distribution favored near logfO2 = -73.6
    r <- revisit(d,"qqr")
    # plot 3: q-q at a single value of logfO2
    basis("O2",-73.6)
    a <- affinity(iprotein=ip)
    d <- diagram(a,logact=0,do.plot=FALSE)
    qqr3 <- revisit(d,"qqr",pch=pch)$H
    legend("topleft",pch=c(1,2,4),legend=myg)
    # plot 4: findit... maximize qqr as a function of O2-H2O-NH3-CO2
    # it shows an optimum at low logaH2O, logaNH3
    f1 <- findit(list(O2=c(-90,-70),H2O=c(-30,0),CO2=c(-20,10),NH3=c(-15,0)),
      "qqr",iprotein=ip,n=8)
    # plot 5: q-q plot at the final loga O2, H2O, CO2, NH3
    # higher correlation coefficient than plot 3
    a <- affinity(iprotein=ip)
    d <- diagram(a,logact=0,do.plot=FALSE)
    qqr5 <- revisit(d,"qqr",pch=pch)$H
    legend("topleft",pch=c(1,2,4),legend=myg)
    # plot 6: trajectory of O2, H2O, CO2, NH3, and the
    # q-q correlation coefficient in the search
    plot(f1,mar=c(2,5,1,1),mgp=c(4,1,0))
  }
}


\keyword{misc}
