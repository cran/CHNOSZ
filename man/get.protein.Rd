\name{get.protein}
\alias{get.protein}
\alias{yeastgfp}
\title{Amino Acid Compositions of Many Proteins}
\description{
Calculate the amino acid compositions of collections of proteins.}

\usage{
  get.protein(protein, organism, abundance = NULL, pname = NULL, 
    average = TRUE, digits = 1) 
  yeastgfp(location, exclusive = TRUE)
}

\arguments{
  \item{protein}{character, name of protein or stress response experiment.}
  \item{organism}{character, organism (\samp{ECO}, \samp{SGD}) or \samp{YeastGFP}.}
  \item{abundance}{numeric, stoichiometry of proteins applied to sums of compositions.}
  \item{pname}{character, names of proteins.}
  \item{average}{logical, return an average composition of the proteins?}
  \item{digits}{numeric, number of decimal places to round the amino acid counts.}
  \item{location}{character, name of subcellular location (compartment).}
  \item{exclusive}{logical, report only proteins exclusively localized to a compartment?}
}

\details{

  When \code{protein} contains one or more Ordered Locus Names (OLN) or Open Reading Frame names (ORF), \code{get.protein} retrieves the amino acid composition of the respective proteins in \emph{Escherichia coli} or \emph{Saccharomyces cerevisiae} (for \code{organism} equal to \samp{ECO} or \samp{SGD}, respectively). The calculation depends on presence of the objects \code{\link{thermo}$ECO} and \code{\link{thermo}$SGD}, which contain the amino acid compositions of proteins in these organisms. If \code{protein} is instead a name of one of the stress response experiments contained in \code{thermo$stress}, e.g. \samp{low.C} or \samp{heat.up}, the function returns the amino acid compositions of the corresponding proteins.

  If the abundances of the proteins are given in \code{abundance}, the individual protein compositions are multiplied by these values then summed into an overall composition; the average is taken if \code{average} is \code{TRUE}, then the amino acid frequencies are rounded to the number of decimal places specified in \code{digits}. The default value of \code{abundance} (1) means the protein compositions are simply summed together. The output of \code{get.protein} can be used as input to \code{\link{add.protein}} to add the proteins to the \code{thermo$protein} data frame in preparation for further calculations. Unless names for the new proteins are given in \code{pname}, they are generated using the values in \code{protein}.

  (NOTE: \code{get.protein} replaces the \code{proteome} function that was present in \pkg{CHNOSZ} up to version 0.7-2. \code{proteome} had a coding error that led to incorrect calculations of the average composition of proteins when \code{abundance} was not equal to 1.)

  The \code{yeastgfp} function returns the identities and abundances of proteins with the requested subcellular localization (specified in \code{location}) using data from the YeastGFP project that is stored in \code{\link{thermo}$yeastgfp}. The default value of \code{exclusive} (\code{FALSE}) tells the function to grab all proteins that are localized to a compartment even if they are also localized to other compartments. If \code{exclusive} is \code{TRUE}, only those proteins that are localized exclusively to the requested compartments are identified, unless there are no such proteins, then the non-exclusive localizations are used (applies to the \samp{bud} localization). The values returns by \code{yeastgfp} can be fed to \code{get.protein} in order to get the amino acid compositions of the proteins.

  \samp{HTCC1062.faa} is a FASTA file of 1354 protein sequences in the organism \emph{Pelagibacter ubique} HTCC1062 downloaded from the NCBI RefSeq collection on 2009-04-12. The specific search term was Protein: txid335992[Organism:noexp] AND "refseq"[Filter].


}

\value{
  For \code{get.protein}, returns the amino acid composition of the specified protein(s) summed together (if \code{single} is \code{TRUE}) or individually (if \code{single} is \code{FALSE}); or, if \code{add} is \code{TRUE}, the index of protein(s) which were added to the \code{thermo$protein} dataframe. \code{yeastgfp} returns a list with elements \code{yORF} and \code{abundance}.
}


\examples{
  \dontshow{data(thermo)}

  ## basic examples of get.protein
  # amino acid composition of two proteins
  get.protein(c('YML020W','YBR051W'),'SGD')
  # average composition of proteins
  get.protein(c('YML020W','YBR051W'),'SGD',
    abundance=1,pname='PROT1_NEW')
  # 1 of one and 1/2 of the other
  get.protein(c('YML020W','YBR051W'),'SGD',
    abundance=c(1,0.5),average=FALSE,pname='PROT2_NEW')
  # compositions of proteins induced in carbon limitation 
  get.protein('low.C','SGD')

  ## overall composition of proteins exclusively localized 
  ## to cytoplasm of S. cerevisiae with reported expression levels
  y <- yeastgfp('cytoplasm')
  p <- get.protein(y$yORF,'SGD',y$abundance,'cytoplasm')
  # add the proteolog and calculate its properties
  i <- add.protein(p)
  protein(i)

  ## Chemical activities of model subcellular proteins
  # (one-dimensional speciation diagram as a function of logfO2)
  basis('CHNOS')
  names <- colnames(thermo$yeastgfp)[6:28]
  # calculate amino acid compositions using 'get.protein' function 
  for(i in 1:length(names)) {
    y <- yeastgfp(names[i])
    p <- get.protein(y$yORF,'SGD',y$abundance,names[i])
    add.protein(p)
  }
  species(names,'SGD')
  res <- 200
  t <- affinity(O2=c(-77,-72,res))
  mycolor <- topo.colors(6)[1:4]
  mycolor <- rep(mycolor,times=rep(6,4))
  oldpar <- par(bg='black',fg='white')
  logact <- diagram(t,balance='PBB',names=names,ylim=c(-4,-1.9),legend.x=NULL,
    color=mycolor,lwd=2,cex.axis=1.5,residue=TRUE)$logact
  # so far good, but how about labels on the plot?
  for(i in 1:length(logact)) {
    imax <- which.max(as.numeric(logact[[i]]))
    adj <- 0.5
    if(imax > 180) adj <- 1
    if(imax < 20) adj <- 0
    text(seq(-77,-72,length.out=res)[imax],logact[[i]][imax],
      labels=names[i],adj=adj)
  }
  title(main=paste('Subcellular proteologs of S. cerevisiae\n',
    describe(thermo$basis[-5,])),col.main=par('fg'))
  par(oldpar)

  ## Oxygen fugacity - activity of H2O predominance 
  ## diagrams for proteomes in 23 YeastGFP localizations
  # arranged by decreasing metastability:
  # order of this list of locations is based on the 
  # (dis)appearance of species on the current set of diagrams
  names <- c('actin','early.Golgi','ER','vacuolar.membrane',
    'cell.periphery','nucleolus','Golgi','lipid.particle',
    'punctate.composite','peroxisome','bud','ER.to.Golgi',
    'nuclear.periphery','ambiguous','late.Golgi','cytoplasm',
    'nucleus','mitochondrion','endosome','vacuole',
    'spindle.pole','bud.neck','microtubule')
  nloc <- c(5,5,5,3,2,3)
  inames <- 1:length(names)
  # define the system
  basis('CHNOS+')
  # calculate amino acid compositions using 'get.protein' function 
  for(i in 1:length(names)) {
    y <- yeastgfp(names[i])
    p <- get.protein(y$yORF,'SGD',y$abundance,names[i])
    add.protein(p)
  }
  species(names,'SGD')
  t <- affinity(H2O=c(-5,0,256),O2=c(-80,-66,256))
  # the plot setup
  layout(matrix(c(1,1,2:7),byrow=TRUE,nrow=4),heights=c(0.7,3,3,3))
  # a title
  par(mar=c(0,0,0,0))
  plot.new()
  text(0.5,0.5,paste('Proteologs for subcellular locations of',
   'S. cerevisiae\n',describe(thermo$basis[-c(2,5),])),cex=1.5)
  opar <- par(mar=c(3,4,1,1),xpd=TRUE)
  for(i in 1:length(nloc)) {
    diagram(t,balance='PBB',names=names[inames],
      ispecies=inames,cex.axis=1.1)
    label.plot(letters[i])
    title(main=paste(length(inames),'locations'))
    # take out the stable species
    inames <- inames[-(1:nloc[i])]
  }
  layout(matrix(1))
  par(opar)

  ### examples for stress response experiments

  # coefficient of variation of relative 
  # abundances of proteins induced in heat 
  # response experiments (Richmond et al., 1999)
  # as a function of fO2 and temperature
  a <- get.protein("heat","ECO")
  add.protein(a)
  basis('CHNOS+')
  species(a$protein,"ECO")
  a <- affinity(T=c(0,150),O2=c(-90,-40))
  d <- diagram(a,residue=TRUE,do.plot=FALSE,mam=FALSE)
  draw.diversity(d)
  title(main="Coefficient of variation of relative abundances
    of proteins in E. coli observed at 50 degC heat shock",
    cex.main=0.9)


  # predominance fields for overall protein 
  # compositions induced by
  # carbon, sulfur and nitrogen limitation
  # (experimental data from Boer et al., 2003)
  expt <- c('low.C','low.N','low.S')
  for(i in 1:length(expt)) {
    p <- get.protein(expt[i],"SGD",abundance=1)
    add.protein(p)
  }
  # thermo set-up
  basis("CHNOS+") 
  basis("O2",-75.29)
  species(expt,"SGD")
  a <- affinity(CO2=c(-5,0),H2S=c(-10,0))
  diagram(a,balance="PBB",names=expt,color=NULL,residue=TRUE)
  title(main=paste("Metastabilities of proteins induced by",
    "carbon, sulfur and nitrogen limitation",sep="\n"))

  # predominance fields for overall protein 
  # compositions induced and repressed in 
  # an/aerobic carbon-limited experiments
  # (Tai et al., 2005)
  # the activities of glucose, ammonium and sulfate
  # are similar to the non-growth-limiting concentrations
  # used by Boer et al., 2003
  basis(c("glucose","H2O","NH4+","H2","SO4-2","H+"),
    c(-1,0,-1.3,999,-1.4,-7))
  # the names of the experiments in thermo$stress
  expt <- c("Clim.aerobic.down","Clim.aerobic.up",
    "Clim.anaerobic.down","Clim.anaerobic.up")
  # here we use abundance to indicate that the protein
  # compositions should be summed together in equal amounts
  for(i in 1:length(expt)) {
    p <- get.protein(expt[i],"SGD",abundance=1)
    add.protein(p)
  }
  species(expt,"SGD")
  a <- affinity(C6H12O6=c(-35,-20),H2=c(-20,0))
  diagram(a,residue=TRUE,color=NULL,as.residue=TRUE)
  title(main="Metastabilities of average protein residues in
    an/aerobic carbon limitation in yeast")

	  
}

\references{

  Boer, V. M., de Winde, J. H., Pronk, J. T. and Piper, M. D. W., 2003. The genome-wide transcriptional responses of \emph{Saccharomyces cerevisiae} grown on glucose in aerobic chemostat cultures limited for carbon, nitrogen, phosphorus, or sulfur. \emph{J. Biol. Chem.}, 278, 3265-3274.



  Richmond, C. S., Glasner, J. D., Mau, R., Jin, H. F. and Blattner, F. R., 1999. Genome-wide expression profiling in \emph{Escherichia coli} K-12. \emph{Nucleic Acids Res.}, 27, 3821-3835.


  Tai, S. L., Boer, V. M., Daran-Lapujade, P., Walsh, M. C., de Winde, J. H., Daran, J.-M. and Pronk, J. T., 2005. Two-dimensional transcriptome analysis in chemostat cultures: Combinatorial effects of oxygen availability and macronutrient limitation in \emph{Saccharomyces cerevisiae}. \emph{J. Biol. Chem.}, 280, 437-447.


}

\keyword{misc}
