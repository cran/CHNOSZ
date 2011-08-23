\name{protein}
\alias{protein}
\alias{protein.residue}
\alias{protein.info}
\alias{residue.info}
\title{Properties of Proteins}
\description{
Retrieve the amino acid compositions or thermodynamic properties and equations of state parameters of proteins.
}

\usage{
  protein(protein,organism=NULL,online=thermo$opt$online,chains=1)
  protein.residue(proteins)
  protein.info(T=25)
  residue.info(T=25)
}

\arguments{
  \item{protein}{character, names of proteins, protein identifiers, or amino acid sequences, or numeric, indices of proteins (rownumbers of \code{\link{thermo}$protein}), or dataframe, protein compositions to sum into new protein.}
  \item{organism}{character, organism identifiers, or physical state.}
  \item{proteins}{character, names of proteins.}
  \item{online}{logical, try an online search if the specified protein(s) are not found locally?}
  \item{chains}{numeric, number of polypeptide chains in added proteins.}
  \item{T}{numeric, temperature in units specified by \code{\link{nuts}}.}
}

\details{

  \code{protein} is a function to query the protein database and to perform group additivity calculations of the standard molal thermodynamic properties and equations of state parameters of proteins. In CHNOSZ, the database of amino acid compositions of proteins is located at \code{thermo$protein} and is populated when the package loaded. See the help for \code{\link{thermo}} for more information. 

  To distinguish names of proteins from those of other species, protein names in CHNOSZ have an underscore ("_") somewhere in their name, as in \samp{LYSC_CHICK}. If a protein name is submitted as a single argument to \code{protein} it is searched for in \code{thermo$protein}; if matches are found, the selected rows are returned. If protein and organism identifiers (e.g. \samp{LYSC} and \samp{CHICK}, respectively) are provided, the rownumbers of matches in \code{thermo$protein} are returned. 

  If no match is found in \code{thermo$protein}, an online search is invoked, unless \code{online} is FALSE. (If \code{online} is NA, the default value of \code{thermo$opt$online}, the user is prompted whether the online search should be performed, and this response is stored in \code{thermo$opt$online}.) The function attempts a search of the SWISS-Prot database (Boeckmann et al., 2003). If the amino acid composition of the protein is successfully retrieved by the online search, that composition is stored in \code{thermo$protein}. 

  If \code{protein} is numeric, the compositional information found in that row(s) of \code{thermo$protein} is combined with sidechain and backbone group contributions to generate the standard molal thermodynamic properties and equations of state parameters of the proteins at 25 \eqn{^{\circ}}{degree }C and 1 bar (Dick et al., 2006), and a dataframe of these values returned. The physical state of the proteins in this calculation is controlled by the value of \code{organism} (\samp{aq} or \samp{cr}; \code{NULL} defaults to \samp{aq}). Note that the properties of aqueous (and crystalline) proteins calculated in this step are hypothetically completely nonionized proteins; the contributions by ionization to the chemical affinities of formation reactions of aqueous proteins can be calculated during execution of \code{\link{affinity}} if the basis species contain \samp{H+}.

  If \code{protein} is a \code{\link{data.frame}}, it is taken to contain the compositions of one or more proteins that are summed to define a new protein whose amino acid composition is returned. In this case, the argument \code{organism} should contain the name of the new protein, e.g. \samp{PROTEIN_NEW}.

  The \code{protein} function modifies the database if \code{organism} is a valid name for a protein (it contains an underscore). The function assumes that \code{protein} contains the amino acid sequence of the new protein to be added to \code{thermo$protein}. The \code{chains} argument specifies the number of polypeptide chains that are in the molecule (0 for amino acid residues, 1 for amino acids and most proteins).

  \code{protein.residue} generates average amino acid residue compositions of proteins. It takes the name(s) of one or more proteins (e.g. \samp{LYSC_CHICK}), retrieves their amino acid compositions from \code{thermo$protein}, and divides by the total number of amino acid residues in each of the proteins.

  \code{protein.info} is a utility to tabulate some properties of proteins. A dataframe is returned containing for each protein that is among the \code{\link{species}} of interest, the name of the protein, its length, formula, and values of the standard molal Gibbs energy of the neutral protein, net charge, standard molal Gibbs energy of the ionized protein, and average oxidation number of carbon. The value of \code{T} indicates the temperature at which to calculate the Gibbs energies and net charge. Net charge and standard molal Gibbs energy of the ionized protein become NA if \samp{H+} is not among the basis species. The values are rounded at a set number of digits for display, and the values of Gibbs energy are in kcal/mol.

  \code{residue.info} calculates the per-residue makeup of the proteins that have been loaded using \code{species}. This amounts to dividing the reaction coefficients in \code{thermo$species} by the length of the protein, but also takes into account the ionization state of the protein if \samp{H+} is one of the basis species. As with \code{protein.info}, the ionization state of the protein is calculated at the pH defined in \code{thermo$basis} and at the temperature specified by the \code{T} argument.

}

\value{
  If \code{protein} is one or more protein names, the matching row(s) of \code{thermo$protein}. If \code{protein} and \code{organism} are protein and organisms identifiers, rownumbers of \code{thermo$protein}. If \code{protein} is numeric, a dataframe with calculated thermodynamic properties and parameters of the neutral protein.
}

\seealso{
  Functions \code{\link{info}}, \code{\link{subcrt}}, \code{\link{species}} and others accept protein names as species arguments. Properties of ionized proteins can be calculated using \code{\link{ionize}} (usually implicitly called as part of a calculation of \code{\link{affinity}}). Compositions of proteins beyond those in \code{thermo$protein} (including yeast and \emph{E. coli}) are also provided in CHNOSZ; see \code{\link{get.protein}} for examples.
}

\examples{
  \dontshow{data(thermo)}
  ### Interaction with the 'protein' function

  ## Thermodynamic properties of proteins
  # get the composition of a protein
  protein("BPT1_BOVIN")
  # retrieve the rownumber of a protein in thermo$protein
  iprotein <- protein("LYSC", "CHICK")
  # calculate properties and parameters of aqueous protein
  protein(iprotein)
  # of crystalline protein
  protein(iprotein, "cr")
  # a call to info() causes the protein properties to
  # be appended to thermo$obigt				
  info("LYSC_CHICK")
  # thermodynamic properties can be calculated with subcrt()
  subcrt("LYSC_CHICK")				
\donttest{
  ### Table of properties of some proteins
  basis("CHNOS+")
  species(c("LYSC_CHICK", "CYC_BOVIN", "MYG_HORSE", "RNAS1_BOVIN"))
  protein.info()
  # the following gives the per-residue composition (i.e. formation 
  # reaction cofficients) for the ionized proteins
  residue.info()

  ## Protein Data from Online Sources
\dontrun{
  ## marked dontrun because it requires internet
  # this asks to search SWISS-Prot
  info("PRND_HUMAN")
  # an online search can also be started from the
  # "subcrt" function
  subcrt("SPRN_HUMAN") } 

  ## Inputting protein compositions
  # make a new protein
  protein("GGSGG", "PROTEIN_TEST")
  # a sequence can be pasted into the command line:
  # type this
  protein("
  # then paste the sequence (this is it)
  # and end the command by typing
  ","PROTEIN_NEW")
  # or use whatever name you want (with an underscore).

  ## Standard molal entropy of a protein reaction
  basis("CHNOS")
  # here we provide the reaction coefficients of the 
  # proteins (per protein backbone); 'subcrt' function calculates 
  # the coefficients of the basis species in the reaction
  s <- subcrt(c("CSG_METTL","CSG_METJA"), c(-1/530,1/530),
    T=seq(0, 350, length.out=50))
  thermo.plot.new(xlim=range(s$out$T), ylim=range(s$out$S),
    xlab=axis.label("T"), ylab=axis.label("DS0r"))
  lines(s$out$T, s$out$S)
  # do it at high pressure as well
  s <- subcrt(c("CSG_METTL","CSG_METJA"), c(-1/530,1/530),
    T=seq(0,350,length.out=50), P=3000)
  lines(s$out$T, s$out$S, lty=2)
  # label the plot
  title(main=paste("Standard molal entropy\n",
    "P = Psat (solid), P = 3000 bar (dashed)"))
  s$reaction$coeff <- round(s$reaction$coeff, 3)
  d <- describe(s$reaction,
    use.name=c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE))
  text(170, -3, c2s(s2c(d,sep="="),sep="\n"), cex=0.8)

  ### Metastability calculations

  ## subcellular homologs of yeast glutaredoxin
  ## as a function of logfO2 - logaH2O, after Dick, 2009
  basis("CHNOS+")
  protein <- c("GLRX1","GLRX2","GLRX3","GLRX4","GLRX5")
  loc <- c("(C)","(M)","(N)","(N)","(M)")
  species(protein,"YEAST")
  a <- affinity(H2O=c(-10,0), O2=c(-85,-60))
  diagram(a, names=paste(protein,loc))
  title(main=paste("Yeast glutaredoxins (black) and residues (blue)\n",
    describe(thermo$basis[-c(2,5),])))
  # note the difference when we set as.residue=TRUE to
  # plot stability fields for the residue equivalents of the
  # proteins instead of the proteins themselves ...
  # the residue equivalent for one of the larger proteins appears
  diagram(a, names=paste(protein,loc), as.residue=TRUE,
    add=TRUE, col="blue")

  ## surface-layer proteins from Methanococcus and others:
  ## a speciation diagram for surface layer proteins
  ## as a function of oxygen fugacity after Dick, 2008
  # make our protein list
  organisms <- c("METSC","METJA","METFE","HALJP","METVO",
    "METBU","ACEKI","BACST","BACLI","AERSA")
  proteins <- c(rep("CSG",6), rep("SLAP",4))
  proteins <- paste(proteins, organisms,sep="_")
  # set some graphical parameters
  lwd <- c(rep(3,6), rep(1,4))
  lty <- c(1:6,1:4)
  # load the basis species and proteins
  basis("CHNOS+")
  species(proteins)
  # calculate affinities
  a <- affinity(O2=c(-100,-65))
  # make diagram
  d <- diagram(a,ylim=c(-5,-1), legend.x=NULL, lwd=lwd,
    ylab=as.expression(quote(log~italic(a[j]))),yline=1.7)
  # label diagram
  text(-80,-1.9,"METJA")
  text(-74.5,-1.9,"METVO")
  text(-69,-1.9,"HALJP")
  text(-78,-2.85,"METBU",cex=0.8,srt=-22)
  text(-79,-3.15,"ACEKI",cex=0.8,srt=-25)
  text(-81,-3.3,"METSC",cex=0.8,srt=-25)
  text(-87,-3.1,"METFE",cex=0.8,srt=-17)
  text(-79,-4.3,"BACST",cex=0.8)
  text(-85.5,-4.7,"AERSA",cex=0.8,srt=38)
  text(-87,-4.25,"BACLI",cex=0.8,srt=30)
  # add water line
  abline(v=-83.1, lty=2)
  title(main=paste("Surface-layer proteins",
    "After Dick, 2008",sep="\n"))

  ## relative metastabilities of bovine proteins, 
  ## as a function of temperature along a glutathione redox buffer
  mod.buffer("GSH-GSSG",c("GSH","GSSG"), logact=c(-3,-7))   
  basis(c("CO2","H2O","NH4+","SO4-2","H2","H+"),
    c(-1,0,-4,-4,"GSH-GSSG",-7)) 
  basis("CO2","gas")
  prot <- c("CYC","RNAS1","BPT1","ALBU","INS","PRIO")
  species(prot,"BOVIN")
  a <- affinity(T=c(0,200))
  d <- diagram(a, as.residue=TRUE, ylim=c(-2,0.5))
  # add some text labels
  text(40,data.frame(d$logact)[25,],prot)
  title(main="Relative stabilities of bovine proteins on glutathione buffer")

  ## relative metastabilities of plasma proteins,
  ## using chemical activities of H2 and O2
  # clean up basis species, species from previous example
  data(thermo)
  basis(c("CO2","NH3","H2S","H2","O2","H+"))
  basis("O2","aq")
  basis(c("CO2","NH3","H2S","H+"),c(-3,-3,-10,-7))
  f <- system.file("extdata/abundance/AA03.csv", package="CHNOSZ")
  pdata <- read.csv(f, as.is=TRUE)
  notna <- !is.na(pdata$name)
  pname <- pdata$name[notna]
  # take out insulin C peptide to show more proteins
  pname <- pname[!pname \%in\% "INS.C"]
  species(pname,"HUMAN")
  a <- affinity(H2=c(-20,0), O2=c(-80,-60))
  diagram(a)
  title(main="Human Plasma Proteins")
  # note that the darker colors go with higher abundances
  # as reported by Anderson and Anderson, 2003
  # add lines showing equilibrium activity of H2O
  species(delete=TRUE)
  species("H2O")
  logaH2 <- seq(-20,0,length.out=128)
  for(logaH2O in c(-5,0,5)) {
    species("H2O",logaH2O)
    a <- affinity(H2=logaH2)
    logaO2 <- diagram(a,what="O2",do.plot=FALSE)$logact[[1]]
    lines(logaH2,logaO2,lty=2)
    itext <- 72 + 5 * logaH2O
    lab <- paste("logaH2O =",logaH2O)
    text(logaH2[itext]+0.4,logaO2[itext],lab,srt=-64)
  }
}
}

\references{

  Anderson, N. L. and Anderson, N. G. (2003) The human plasma proteome: History, character and diagnostic prospects (Vol. 1 (2002) 845-867). \emph{Molecular and Cellular Proteomics} \bold{2}, 50. \url{http://dx.doi.org/10.1074/mcp.A300001-MCP200}

  Boeckmann, B., Bairoch, A., Apweiler, R., Blatter, M.-C., Estreicher, A., Gasteiger, E., Martin, M. J., Michoud, K., Donovan, C., Phan, I., Pilbout, S. and Schneider, M. (2003) The SWISS-PROT protein knowledgebase and its supplement TrEMBL in 2003. \emph{Nucleic Acids Res.} \bold{31}, 365--370. \url{http://www.expasy.org}, accessed on 2007-12-19.

  Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \url{http://www.biogeosciences.net/3/311/2006/bg-3-311-2006.html}

  Dick, J. M. (2008) Calculation of the relative metastabilities of proteins using the CHNOSZ software package. \emph{Geochem. Trans.} \bold{9}:10. \url{http://dx.doi.org/10.1186/1467-4866-9-10}

}

\keyword{primary}
