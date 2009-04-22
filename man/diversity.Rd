\name{diversity}
\alias{grep.file}
\alias{read.fasta}
\alias{diversity}
\alias{richness}
\alias{extremes}
\alias{where.extreme}
\alias{draw.diversity}
\title{Calculate Diversity Indices for Chemical Activities}
\description{
  Search a FASTA file for protein sequences, read protein sequences from a file and count numbers of amino acids in each sequence, calculate species richness or standard deviations or coefficient of variation of activities or logarithms of activities of chemical species. Plot the diversity values. 
}

\usage{
  grep.file(file, x = "", y = NULL, ignore.case = TRUE, startswith = ">")
  read.fasta(file, i = NULL, aa = TRUE)
  diversity(logact, which = "cv", as.log = FALSE)
  richness(logact, logactmin = -4)
  extremes(z, what = "minimum")
  where.extreme(z, what= "minimum")
  draw.diversity(d, which = "cv", logactmin = -4, col = par("fg"),
    as.log = FALSE, yline = NULL, ylim = NULL, ispecies = NULL,
    add = FALSE)
}

\arguments{
  \item{file}{character, path to a FASTA amino acid sequence file.}
  \item{x}{character, term to search in sequence headers.}
  \item{y}{character, term to exclude in searching sequence headers.}
  \item{ignore.case}{logical, ignore differences between upper- and lower-case?}
  \item{startswith}{character, only lines starting with this expression are matched.}
  \item{i}{numeric, line numbers of sequence headers for protein sequences to read.}
  \item{aa}{logical, count amino acids in sequences and format as dataframe?}
  \item{logact}{list, logarithms of activities of species.}
  \item{d}{list, output from \code{\link{diagram}}.}
  \item{which}{character, type of diversity statistic to calculate.}
  \item{as.log}{logical, calculate standard deviations of logarithms of activities?}
  \item{logactmin}{numeric, cutoff value for a species to be counted in richness calculations.}
  \item{z}{numeric, matrix of values.}
  \item{what}{character, \samp{minimum} or \samp{maximum}.}
  \item{col}{character, color to use for line plots.}
  \item{yline}{numeric, margin line for y-axis label.}
  \item{ylim}{numeric, limits of y axis.}
  \item{ispecies}{numeric, which species to consider.}
  \item{add}{logical, add to an existing plot?}
}

\details{

  \code{grep.fasta} and \code{read.fasta} are used to search for and retrieve entries in a FASTA file. \code{grep.fast} takes a search term in \code{x} and optionally a term to exclude in \code{y}. The \code{ignore.case} option is passed to \code{\link{grep}}, which does the work of finding lines that match. Only lines that start with the expression in \code{startswith} are searched; the default setting reflects the format of the header line for each sequence in a FASTA file.

  \code{diversity} and \code{richness} both take as input a list of logarithms of activities of chemical species calculated using \code{\link{diagram}}. The calculations available in \code{diversity} are the Shannon index, standard deviation, and coefficient of variation, which are specified in \code{which} by \samp{shannon}, \samp{sd} and \samp{cv}, respectively. By default the \emph{activities} of species are used in the calculations, but if \code{as.log} is set to \code{TRUE}, the \emph{logarithms} of activities are used instead.

  Given a matrix of numeric values in \code{z}, \code{extremes} locates the maximum or minimum values in both dimensions. That is, the x values that are returned are the column numbers where the extreme is found for each row, and the y values that are returned are the row numbers where the extreme is found for each column. \code{where.extreme} takes a matrix and returns the row and column numbers where the maximum or minimum value is found. If the extreme value is repeated, the row and column numbers for all instances are returned.

  \code{draw.diversity} takes output from \code{\link{diagram}} that includes logarithms of activities of species, calls either \code{diversity} or \code{richness}, plots the results. For systems of two variables \code{where.extreme} (and \code{extremes}) are used to indicate the location of the maximum Shannon index or richness (and the ridges around it), or minimum standard deviation or coefficient of variation (and the valleys around it). The ridges and valleys are plotted as dashed lines and are colored green for the x values returned by \code{extremes} and blue for the y values returned by \code{extremes}.

}

\value{

  \code{grep.fasta} returns numeric, indicating the line numbers of the file that match the search constraints. \code{read.fasta} returns a list of sequences (if \code{aa} is \code{TRUE}) or a data frame with amino acid compositions of proteins, with columns corresponding to those in \code{\link{thermo}$protein}. \code{diversity} and \code{richness} return a numeric (vector or matrix) that has the same dimensions as a single element of \code{logact}. \code{extremes} and \code{where.extreme} both return a list of x and y values. \code{draw.diversity} returnes a list with elements \code{H} (the result of either \code{diversity} or \code{richness}), \code{ixmin} and \code{iymin} (the column and row numbers where the maximum or minimum value is found), \code{xmin} and \code{ymin} (the values of the x and y variables where the maximum or minimum is found) and \code{extval} (the value of \code{H} at the maximum or minimum).

}

\seealso{
  Use \code{\link{add.protein}} to add the amino acid compositions returned by \code{read.fasta} to the data frame in \code{\link{thermo}$protein}. When computing relative abundances of many proteins that might be found with \code{grep.fasta}, consider using the \code{iprotein} arugment of \code{\link{affinity}} to speed things up.
}

\examples{
  \dontshow{data(thermo)}

  ### diversity plots
  
  ## carboxylases from different organisms
  # ribulose phosphate carboxylase/oxygenase
  rubisco <- c('Q6JAI0_9RHOD','A5CKC7_9CHRO','RBL_SYNJA','A1E8R4_9CHLO',
    'A8C9T6_9MYCO','A3EQE1_9BACT','A6YF84_9PROT', 'RBL_BRAJA',
    'A1RZJ5_THEPD','RBL_METJA','A3DND9_STAMF','RBL_PYRHO') 
  # acetyl-coenzyme A carboxylase
  accoaco <- c('Q05KD0_HYDTH','Q9F7M8_PRB01','A6CDM2_9PLAN','A0GZU2_9CHLR',
    'ACCA_DEIRA','A1VC70_DESVV','A7WGI1_9AQUI','Q2JSS7_SYNJA',
    'A4AGS7_9ACTN','ACCA_AQUAE','ACCA_CAUCR','A6VIX9_METM7')  
  # calculate affinities as a function of T with
  # buffered logfH2 and CO2
  # adjust the glutathione buffer 
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

  ### using grep.file, read.fasta, add.protein
  # calculations for Pelagibacter ubique
  f <- system.file("HTCC1062.faa",package="CHNOSZ")
  # line numbers of all entries in the file
  j <- grep.file(f)  # length = 1354
  # locate entries whose names contain DNA
  j <- grep.file(f,"DNA")
  # get the amino acid compositions of these protein
  p <- read.fasta(f,j)
  # add these proteins to CHNOSZ's inventory
  i <- add.protein(p)
  # set up a thermodynamic system
  basis('CHNOS+')
  # calculate affinities in logfO2-pH space
  a <- affinity(pH=c(0,12),O2=c(-80,-65),iprotein=i)
  # calculate the logarithms of activities
  d <- diagram(a,residue=TRUE,do.plot=FALSE,mam=FALSE)
  # show the protein richness
  draw.diversity(d,"richness",logactmin=-4)
  title(main=paste("Predicted protein richness in",
   "Pelagibacter ubique",sep="\n"))
  # show the coefficient of variation of activities
  draw.diversity(d,"cv")
  title(main=paste("Coefficient of variation of protein",
    "activities in P. ubique",sep="\n"))


  ### the rest is a rather long example
  ### comparing calculated and observed abundances
  ### of proteins and organisms in two different environments
  # (the calculations use a buffer of proteins from the 
  # four most abundant organisms to predict the occurrence 
  # of other species in a metastable population)
  ## bacterial species and abundances from Spear et al., 2005 (Fig. 4)
  ORGANISMS <- c('Aquificales','Thermotogales','Bacillus','Thermodesulfobacteria',
    'Thermus/Deinococcus','Hydrogenobacter', 'Hydrogenobaculum','Bacteroidetes',
    'a-proteobacterium','d-proteobacterium','g-proteobacterium','OD-1')
  # which species are most abundant in low and high sulfur environments
  # the first three of each are aquificales, thermotogales, thermodesulfobacteria
  lowS.which <- c(1,2,4,3,5)
  lowS.values <- c(75,12,4,7,2)
  highS.which <- c(1,2,4,11,6,9,7,8,10,13)
  highS.values <- c(63,2,10,2,5,2,9,2,3,1)
  # what chemical species we use to represent each organism
  PROTEINS <- c('ACCA_AQUAE','A8F7V4_THELT','RBL1_THIDA',
    'Q93EV7_9BACT','ACCA_DEIRA','Q05KD0_HYDTH','A7WGI1_9AQUI','Q64VW6_BACFR',
    'ACCA_CAUCR','A1VC70_DESVV','ACCA_PSEAE')
  ## to make pie charts for calculated protein abundances
  plot.pie.calc <- function(which="high",T=25,main='') {
    # first clean up the buffer definition in case we have
    # been run already
    thermo$buffer <<- thermo$buffer[thermo$buffer$name!='PROTEINS',]
    # we take four predominant proteins from SWM+05
    myprot <- PROTEINS[get(paste(which,"S.which",sep=""))][1:4]
    mypercent <- get(paste(which,"S.values",sep=""))[1:4]
    # use these four proteins to create a buffer
    mybufprot <- paste(myprot,'RESIDUE',sep='.')
    mod.buffer('PROTEINS',mybufprot,'aq',log10(mypercent/100))
    # our species are the residues of all proteins
    species(delete=TRUE)
    add.protein(protein.residue(PROTEINS))
    species(paste(PROTEINS,"RESIDUE",sep="."))
    # assign the buffer to three basis species
    basis(c('CO2','H2O','NH3'),'PROTEINS')
    # calculate the buffered activities
    a <- affinity(return.buffer=TRUE,balance=1,T=T)
    # make the titles
    sub <- c2s(paste(names(a)[1:3],round(as.numeric(a)[1:3])))
    main <- paste('\n',main,'\n',sub)
    # set the total species activities to those in the buffer
    species(1:nrow(species()),-99)
    species(mybufprot,log10(mypercent/100))
    # get the activities back to numeric
    basis(names(a)[1:3],as.numeric(a)[1:3])
    thermo$basis$logact <<- as.numeric(thermo$basis$logact)
    # colors
    col <- rep('white',99)
    col[match(myprot,PROTEINS)] <- heat.colors(4)
    # calculate the distribution of species
    mylogaH2O <- thermo$basis$logact[rownames(thermo$basis)=='H2O']
    a <- affinity(H2O=c(mylogaH2O,mylogaH2O-1,2),T=T)
    logacts <- diagram(a,do.plot=FALSE,residue=TRUE,balance=1)$logact
    # assemble the names and logarithms of activities
    # of species that are above a minimum value
    names <- character()
    values <- numeric()
    cols <- character()
    logactmin <- -2
    for(i in 1:length(logacts)) {
      myvalue <- logacts[[i]][1]
      if(myvalue > logactmin) {
        names <- c(names,ORGANISMS[i])
        values <- c(values,myvalue)
        cols <- c(cols,col[i])
      }
    }
    # remove the logarithms
    values <- 10^values
    # sort them by abundance
    isort <- sort(values,index.return=TRUE,decreasing=TRUE)$ix
    names <- names[isort]
    values <- values[isort]
    cols <- cols[isort]
    # make a pie chart 
    pie(values,names,clockwise=FALSE,main=main,col=cols,radius=0.7)
  }
  ## to plot pie charts for observed abundances of organisms
  plot.pie.obs <- function(which="low") {
    # the values from SWM+05
    names <- ORGANISMS[get(paste(which,"S.which",sep=""))]
    values <- get(paste(which,"S.values",sep=""))
    main <- paste("observed at",which,"H2S")
    # colors for the four dominant species
    mycol <- heat.colors(4)
    # colors for the rest
    mycol <- c(mycol,rep('white',length(names)-length(mycol)))
    # sort the values
    isort <- sort(values,index.return=TRUE,decreasing=TRUE)$ix
    values <- values[isort]
    names <- names[isort]
    mycol <- mycol[isort]
    pie(values,names,clockwise=FALSE,main=main,col=mycol,radius=0.7)
  }
  ## to plot pie diagrams showing calculated protein activities
  ## and observed abundances of organisms
  plot.pie <- function(T=80) {
    # first deal with the layout
    layout(matrix(1:4,byrow=TRUE,nrow=2))
    par(mar=c(1,1,1,1))
    thermo$basis <<- NULL
    # basis definition
    basis(c('CO2','H2O','NH3','H2S','H2','H+'))
    basis('pH',7)
    val <- function(text) round(as.numeric(
      thermo$basis$logact[rownames(thermo$basis)==text]))
    # now to plotting
    # low sulfur and relatively oxidizing
    basis(c('H2','H2S'),c(-9,-9))
    plot.pie.calc("low",T=T,
      main=paste("calculated for H2",val("H2"),"H2S",val("H2S")))
    label.plot('a')
    plot.pie.obs("low")
    label.plot('b')
    # high sulfur and relatively reducing
    basis(c('H2','H2S'),c(-5,-3))
    plot.pie.calc("high",T=T,
      main=paste("calculated for H2",val("H2"),"H2S",val("H2S")))
    label.plot('c')
    plot.pie.obs("high")
    label.plot('d')
  }
  ## now run it!
  plot.pie(80)
  # notably, 80 deg C gives a computed result that is more similar
  # to the observed species abundances than at 25 deg C ...
  # plot.pie(25)

  \dontshow{
    layout(matrix(1))
  }


}

\references{

  Spear, J. R., Walker, J. J., McCollom, T. M. and Pace, N. R., 2005. Hydrogen and bioenergetics in the Yellowstone geothermal ecosystem. \emph{Proc. Natl. Acad. Sci. U. S. A.}, 102, 2555-2560.

}

\keyword{misc}
