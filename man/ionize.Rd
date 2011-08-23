\name{ionize}
\alias{ionize}
\title{Ionization Properties of Proteins}
\description{
  Calculate the charges of proteins and contributions of ionization to their chemical affinities of formation reactions or other thermodynamic properties. 
}

\usage{
  ionize(affinity = NULL, other = NULL)
}

\arguments{
  \item{affinity}{affinities of formation of ionizable groups and nonionized proteins.}
  \item{other}{other properties of ionizable groups and nonionized proteins.}
}

\details{

  Most of the time, the user doesn't need to call this function directly. It is called by \code{\link{affinity}} if proteins are among the species of interest, \samp{H+} is in the basis, and \code{thermo$opt$ionize} is TRUE. This function can also be used from the top level to calculate the ionization state and ionization contributions to thermodynamic properties other than the chemical affinity (see examples below). The ionization model used by this function takes account of the ionization constants of amino acid sidechain groups and the terminal groups in proteins as a function of temperature and pressure calculated using the properties and parameters taken from Dick et al., 2006.

  \code{\link{ionize}} responds to different argument configurations. \code{affinity} and \code{other} both refer to the list contained in the \code{values} element of the output of \code{\link{affinity}}. If both these arguments are missing, \code{\link{species}} is called to load the ionizable sidechain and backbone groups of proteins into the global species list. If \code{other} is missing, the function calculates and outputs the affinities of formation of the ionized proteins. If \code{other} is present, it contains the values of other properties, for example \samp{Cp}, for which ionization values are calculated. Finally, if \code{other} is identical to \code{affinity}, the function outputs the calculated net charges (overall ionization states) of the proteins. Note that \code{ionize} only computes values for proteins; those for other species are returned unaltered in the output.

}

\value{
  List, each element of which corresponds to the affinity of the formation reaction or other property of an ionized protein.
}

\seealso{
  \code{\link{protein.info}} and \code{\link{residue.info}} use this function to display ionization states of proteins.
}

\examples{
  \dontshow{data(thermo)}
  ### Direct Interaction with 'ionize'
 
  ## Charge of LYSC_CHICK as a function of pH and T
  # After Fig. 10 of Dick et al., 2006
  basis(c("CO2","H2O","NH3","H2S","O2","H+"),rep(999,6))
  species("LYSC_CHICK")
  # add the ionizable groups
  ionize()
  # get the affinities and charges (along equal temperature increments)
  T <- c(25,150,6); pH <- c(0,14,50)
  x <- affinity(pH=pH,T=T)
  z <- ionize(x$values,x$values)
  # plot charges at the temperatures we're interested in
  plot(z[[1]][,1],x=seq(0,14,length.out=50),type="l",xlab="pH",
    ylab="net charge (Z)")  # 25  deg C
  lines(z[[1]][,4],x=seq(0,14,length.out=50),col="red")    # 100 deg C
  lines(z[[1]][,6],x=seq(0,14,length.out=50),col="orange") # 150 deg C
  text(x=c(12,10,9),y=c(-15,-16,-18),labels=paste("T=",c(25,100,150),sep=""))
  # if cysteine is oxidized (to cystine disulfide bonds) it may not ionize. 
  # suppress its ionization by upping energy of the ionized group
  mod.obigt("[Cys-]",G=999999)
  x <- affinity(pH=pH,T=25)
  z <- ionize(x$values,x$values)
  lines(z[[1]],x=seq(0,14,length.out=50),lty=2)
  text(x=12,y=-7,"T=25, oxidized")
  # add experimental points
  RT71 <- read.csv(system.file("extdata/cpetc/RT71.csv",package="CHNOSZ"))
  points(RT71$pH,RT71$Z)
  title(main=paste("Ionization of unfolded LYSC_CHICK\n",
    "Experimental points at 25 degC from Roxby and Tanford, 1971"),
     cex.main=0.9)
  # forget our changes to thermo$obigt for next examples
  data(thermo)
\donttest{
  ## Heat capacity of LYSC_CHICK as a function of T
  basis("CHNOS+"); species("LYSC_CHICK")
  pH <- c(5,9,3); T <- c(0,100,25)
  T.values <- seq(T[1],T[2],length.out=T[3])
  # add the ionizable groups
  ionize()
  a <- affinity(pH=pH,T=T)
  # values for non-ionized protein
  c <- affinity(pH=pH,T=T,property="Cp.species")
  plot(T.values,c$values[[1]][1,],
    xlab=axis.label("T"),ylab=axis.label("Cp"),
    ylim=c(5000,8000),type="l",mgp=c(2.2,0.2,0)) 
  # values for ionized protein
  cp <- ionize(a$values,c$values)
  for(i in 1:3) {
    lines(seq(T[1],T[2],length.out=T[3]),cp[[1]][i,],lty=2)
    text(80,cp[[1]][i,][21],paste("pH=",pH[i],sep=""))
  }
  # Makhatadze and Privalov's group contributions
  T.MP <- c(5,25,50,75,100,125)
  points(T.MP,convert(MP90.cp(T.MP,"LYSC_CHICK"),"cal"))
  # Privalov and Makhatadze's experimental values
  e <- read.csv(system.file("extdata/cpetc/PM90.csv",package="CHNOSZ"))
  e <- e[e$protein=="LYSC_CHICK",]
  points(e$T,convert(e$Cp,"cal"),pch=16)
  title(main=paste("Calc'd heat capacity of LYSC_CHICK:",
    "non/ionized(solid/dashed);\n",
    "Makhatadze+Privalov 1990 (open, calc; filled, expt)"),
    cex.main=0.9)

  ### Metastability calculations using 'ionize'

  ## Eh-pH diagrams for intra/extracellular proteins
  organism <- c("PYRFU","ECOLI","YEAST")
  intracellular <- c("AMPM","AMPM","AMPM1")
  extracellular <- c("O08452","AMY1","PST1")
  basis("CHNOSe")  # for Eh we need electrons
  mycol <- c("red","green","blue")
  for(i in 1:3) {
    species(delete=TRUE)
    species(c(intracellular[i],extracellular[i]),organism[i])
    if(i == 1) add <- FALSE else add <- TRUE
    a <- affinity(pH=c(0,14),Eh=c(-1,0))
    diagram(a,add=add,color=NULL,names=species()$name,
      col=mycol[i],col.names=mycol[i])
  }
  title(main=paste("Intracellular (AMPM) and extracellular proteins\n",
    describe(thermo$basis[1:4,])))
  
  ## Buffer + ionization: Metastablilities of
  ## thiol peroxidases from model bactera
  ## (ECOLI, BACSU mesophile; AQUAE thermophile,
  ## THIDA acidophile, BACHD alkaliphile)
  basis("CHNOS+")
  organisms <- c("ECOLI","AQUAE","BACSU","BACHD","THIDA")
  species("TPX",organisms)
  # create a buffer with our proteins in it
  mod.buffer("TPX",paste("TPX",organisms,sep="_"))
  # set up the buffered activities
  basis(c("CO2","H2O","NH3","O2"),"TPX")
  a <- affinity(return.buffer=TRUE,T=50)
  basis(c("CO2","H2O","NH3","O2"),as.numeric(a[1:4]))
  a <- affinity(pH=c(0,14,200),T=c(25,80,200))
  diagram(a,color=NULL)
  title(main=paste("Thiol peroxidases from bacteria\n",
    describe(thermo$basis[-6,],T=NULL)),cex.main=0.9)

  ## Buffer + ionization: Metastable assemblage
  ## for E. coli sigma factors on a T-pH diagram
  # (sigma factors 24, 32, 38, 54, 70, i.e.
  # RpoE, RpoH, RpoS, RpoN, RpoD)
  proteins <- c("RPOE","RP32","RPOS","RP54","RPOD")
  basis("CHNOS+")
  basis("pH",7.4)
  # define and set the buffer
  change("_sigma",paste(proteins,"ECOLI",sep="_"))
  basis(c("CO2","NH3","H2S","O2"),"sigma")
  logact <- affinity(return.buffer=TRUE,T=25)
  # Set the activities of the basis species to constants 
  # corresponding to the buffer, and diagram the relative
  # stabilities as a function of T and pH
  basis(c("CO2","NH3","H2S","O2"),as.numeric(logact))
  species(paste(proteins,"ECOLI",sep="_"))
  a <- affinity(pH=c(5,10),T=c(10,40))
  diagram(a,balance="PBB",residue=FALSE)
  title(main=paste("Sigma factors in E. coli\n",
    describe(thermo$basis[-6,],T=NULL)),cex.main=0.95)
  }
}

\references{

  Dick, J. M., LaRowe, D. E. and Helgeson, H. C. (2006) Temperature, pressure, and electrochemical constraints on protein speciation: Group additivity calculation of the standard molal thermodynamic properties of ionized unfolded proteins. \emph{Biogeosciences} \bold{3}, 311--336. \url{http://www.biogeosciences.net/3/311/2006/bg-3-311-2006.html}

  Makhatadze, G. I. and Privalov, P. L. (1990) Heat capacity of proteins. 1. Partial molar heat capacity of individual amino acid residues in aqueous solution: Hydration effect. \emph{J. Mol. Biol.} \bold{213}, 375--384. \url{http://dx.doi.org/10.1016/S0022-2836(05)80197-4}

  Privalov, P. L. and Makhatadze, G. I. (1990) Heat capacity of proteins. II. Partial molar heat capacity of the unfolded polypeptide chain of proteins: Protein unfolding effects. \emph{J. Mol. Biol.} \bold{213}, 385--391. \url{http://dx.doi.org/10.1016/S0022-2836(05)80198-6}

}

\keyword{secondary}
