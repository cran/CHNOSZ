\name{protein}
\alias{protein}
\alias{protein.residue}
\alias{protein.info}
\title{Properties of Proteins}
\description{
Retrieve the amino acid compositions or thermodynamic properties and equations of state parameters of proteins.
}

\usage{
  protein(protein,organism=NULL,online=thermo$opt$online)
  protein.residue(proteins)
  protein.info()
}

\arguments{
  \item{protein}{character, protein identifiers, or numeric, indices of protein (rownumbers of \code{thermo$protein}), or dataframe, protein compositions to sum into new protein.}
  \item{organism}{character, organism identifiers (required if \code{protein} is character), or physical state (optional if \code{protein} is numeric).}
  \item{proteins}{character, names of proteins.}
  \item{online}{logical, try an online search if the specified protein(s) are not found locally?}
}

\details{

  \code{protein} is a function to query the protein database and to perform group additivity calculations of the standard molal thermodynamic properties and equations of state parameters of proteins.

  The user will generally specify a protein by submitting the name of one in the argument to a function like \code{\link{species}} or \code{\link{subcrt}}. To distinguish names of proteins from those of other species, protein names in \pkg{CHNOSZ} have an underscore ("\_") somewhere in their name, as in \samp{LYSC\_CHICK}.

  If two character arguments are provided and neither one contains an underscore, a search in local and perhaps online data sources is initiated. For each \code{protein}--\code{organism} pair (the arguments should be the same length), the contents of \code{\link{thermo}$proteins} are searched for completely matching (both \code{protein} and \code{organism}) entries. 

  If no match is found in \code{thermo$proteins}, an online search is invoked, unless \code{online} is FALSE. (If \code{online} is NA, i.e., the default setting in \code{thermo$opt$online}, the user is prompted whether she or he wishes the search to be performed, and this response is stored in \code{\link{thermo}$opt$online}.) The function attempts a search of the SWISS-Prot database (Boeckmann et al., 2003). The search string in this case is formed by joining the corresponding elements of the two arguments with an intervening underscore character to make a name such as \samp{LYSC\_CHICK}. If the amino acid composition of the protein is successfully retrieved by the online searche, that composition is stored in \code{thermo$proteins}. For either local or online matches, the values returned by \code{protein} are the row numbers of the protein composition in \code{thermo$proteins}.

  If \code{protein} is numeric, the compositional information found in that row(s) of \code{thermo$proteins} is combined with sidechain and backbone group contributions to generate the standard molal thermodynamic properties and equations of state parameters of the proteins at 25 \eqn{^{\circ}}{degree }C and 1 bar (Dick et al., 2006), and a dataframe of these values returned. The physical state of the proteins in this calculation is controlled by the value of \code{organism} (\samp{aq} or \samp{cr}; \code{NULL} defaults to \samp{aq}). Note that the properties of aqueous (and crystalline) proteins calculated in this step refer to hypothetically completely nonionized proteins; the contributions by ionization to the chemical affinities of formation reactions of aqueous proteins can be calculated during execution of \code{\link{affinity}} if the basis species contain \samp{H+} (see \code{\link{ionize}}).

  If \code{protein} is character but \code{organism} is NULL, the function assumes that \code{protein} refers to the name of protein, that is searched for in \code{thermo$proteins}; if matches are found, the selected rows are returned. If \code{organism} looks like the name of a protein (it contains an underscore), the function assumes that \code{protein} contains the amino acid sequence and of a new protein, and the corresponding amino acid composition is added to \code{thermo$protein}, with the name given by \code{organism}. This allows for entry of protein compositions at the command line.

  If \code{protein} is a \code{\link{data.frame}}, it is taken to be representative of the compositions of one or more proteins that are summed to make a new protein. In this case, the argument \code{organism} should contain the name of the new protein, e.g. \samp{PROTEIN\_NEW}.

  \code{protein.residue} generates average residue compositions of proteins. It takes the name(s) of one or more proteins (e.g. \samp{LYSC\_CHICK}), retrieves their amino acid compositions from \code{\link{thermo}$protein}, and divides by the total number of amino acids in each protein.

  \code{protein.info} is a utility to tabulate some properties of proteins. A dataframe is returned containing for each protein that is among the \code{\link{species}} of interest, the name of the protein, its length, formula, and values of the standard molal Gibbs energy of the neutral protein, net charge, standard molal Gibbs energy of the ionized protein, and nominal carbon oxidation state. 


}

\value{
  If \code{protein} or \code{organism} contains an underscore, a row of the protein composition dataframe. Otherwise, if \code{protein} is numeric, a dataframe with calculated thermodynamic properties and parameters of the neutral protein, or if \code{protein} is character, an \code{\link{invisible}} return of numeric values (or NA for no match) representing rownumbers of \code{thermo$protein} that were matched to \code{protein}--\code{organism} identifiers. 
}

\seealso{
  \code{\link{get.protein}} for retrieving compositions of proteins in yeast and E. coli, including those identified in stress response experiments. \code{\link{add.protein}} for adding these compositions to the dataframe that is accessed by \code{protein}.
}

\examples{
  \dontshow{data(thermo)}
  
  ### Interaction with the 'protein'function

  ## Thermodynamic properties of proteins
  # get the composition of a protein
  protein('BPT1_BOVIN')
  # retrieve the rownumber of a protein in thermo$protein
  iprotein <- protein('LYSC','CHICK')
  # calculate properties and parameters from group additivity	
  protein(iprotein)
  # a call to info() causes the protein properties to
  # be appended to thermo$obigt				
  info('LYSC_CHICK')
  # the second time it is faster				
  info('LYSC_CHICK')
  # thermodynamic properties can be calculated with subcrt()
  subcrt('LYSC_CHICK')				

  ### Table of properties of some proteins
  basis('CHNOS+')
  species(c('LYSC_CHICK','CYC_BOVIN','MYG_HORSE','RNAS1_BOVIN'))
  # here, G in the Gibbs energy of a neutral protein, Z is the
  # charge of an ionized protein, G.Z is the Gibbs energy of the 
  # ionized protein, and Z.C is the nominal carbon oxidation state
  protein.info()
  
  ## Protein Data from Online Sources
  \dontrun{
    ## marked dontrun because it requires internet
    # this asks to search SWISS-Prot
    info('PRND_HUMAN')
    # an online search can also be started from the
    # 'subcrt' function
    subcrt('SPRN_HUMAN')
  }  ## end dontrun

  ## Inputting protein compositions
  # make a new protein
  protein('GGSGG','PROTEIN_TEST')
  # a sequence can be pasted into the command line:
  # type this
  protein('
  # then paste the sequence
  # and end the command by typing
  ','PROTEIN_NEW')
  # or use whatever name you want (with an underscore).

  ## Standard molal entropy of a protein reaction
  basis('CHNOS')
  # here we provide the reaction coefficients of the 
  # proteins (per protein backbone); 'subcrt' function calculates 
  # the coefficients of the basis species in the reaction
  t <- subcrt(c('CSG_METTL','CSG_METJA'),c(-1/558,1/530),
    T=seq(0,350,length.out=15))
  thermo.plot.new(xlim=range(t$out$T),ylim=range(t$out$S),
    xlab=axis.label('T'),ylab=axis.label('DS0r'))
  lines(t$out$T,t$out$S)
  # do it at high pressure as well
  t <- subcrt(c('CSG_METTL','CSG_METJA'),c(-1/558,1/530),
    T=seq(0,350,length.out=15),P=3000)
  lines(t$out$T,t$out$S,lty=2)
  # label the plot
  title(main=paste('Standard molal entropy\n',
    'P = Psat (solid), P = 3000 bar (dashed)'))
  t$reaction$coeff <- round(t$reaction$coeff,3)
  d <- describe(t$reaction,
    use.name=c(TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE))
  text(160,-8,c2s(s2c(d,sep='=',move.sep=TRUE),sep='\n'),cex=0.8)


  ### Metastability calculations

  ## sigma factors of E. coli as a function
  ## of logfO2 - logaNH3
  basis('CHNOS') 
  species(c('RPOE','RP32','RP54','RPOD'),'ECOLI')
  t <- affinity(NH3=c(-10,0),O2=c(-80,-75))
  diagram(t,balance='PBB',cex.axis=1.5)
  title(main=paste('Sigma factors of E. coli\n',
    describe(thermo$basis[-c(3,5),])))
  
  ## subcellular homologs of yeast glutaredoxin
  ## as a function of logfO2 - logaH2O
  basis('CHNOS')
  protein <- c('GLRX1','GLRX2','GLRX3','GLRX4','GLRX5')
  loc <- c('(C)','(M)','(N)','(N)','(M)')
  species(protein,'YEAST')
  t <- affinity(H2O=c(-10,0),O2=c(-85,-60))
  diagram(t,names=paste(protein,loc),cex.axis=1.5)
  title(main=paste('Subcellular homologs of yeast glutaredoxin\n',
    describe(thermo$basis[-c(2,5),])))


  ## surface-layer proteins from Methanococcus spp.:
  ## METVO (mesophile)
  ## METTL (thermophile)
  ## METJA (hyperthermophile)
  # a speciation diagram for surface layer proteins
  # as a function of oxygen fugacity
  # after Dick, 2008
  # make our protein list
  organisms <- c("METSC","METJA","METFE","HALJP","METVO",
    "METBU","ACEKI","BACST","BACLI","AERSA")
  proteins <- c(rep("CSG",6),rep("SLAP",4))
  proteins <- paste(proteins,organisms,sep="_")
  # set some graphical parameters
  lwd <- c(rep(3,6),rep(1,4))
  lty <- c(1:6,1:4)
  # load the basis species and proteins
  basis("CHNOS+")
  species(proteins)
  # calculate affinities
  a <- affinity(O2=c(-100,-65))
  # make diagram
  d <- diagram(a,ylim=c(-5,-1),residue=TRUE,legend.x=NULL,lwd=lwd,
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
  abline(v=-83.1,lty=2)
  title(main=paste("Surface-layer proteins","After Dick, 2008",sep="\n"))
  ## now, show the species richness 
  draw.diversity(d,"richness",logactmin=-4)
  title(main=paste("Surface-layer protein richness"))

  
}

\references{

  Boeckmann, B., Bairoch, A., Apweiler, R., Blatter, M.-C., Estreicher, A., Gasteiger, E., Martin, M. J., Michoud, K., Donovan, C., Phan, I., Pilbout, S. and Schneider, M., 2003. The SWISS-PROT protein knowledgebase and its supplement TrEMBL in 2003. \emph{Nucleic Acids Res.}, 31, 365-370. \url{http://www.expasy.org}, accessed on 2007-12-19.

  Dick, J. M., LaRowe, D. E. and Helgeson, H. C., 2006. Temperature, pressure and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences}, 3, 311-336.


}

\keyword{misc}
