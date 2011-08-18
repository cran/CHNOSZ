\name{util.formula}
\alias{util.formula}
\alias{GHS}
\alias{element}
\alias{expand.formula}
\alias{ZC}
\alias{residue.formula}
\alias{protein.formula}
\title{Functions to Work with Chemical Formulas}

\description{
  Calculate the standard molal entropy of elements in a compound; calculate the standard molal Gibbs energy or enthalpy of formation, or standard molal entropy, from the other two; list coefficients of selected elements in a chemical formula; calculate the average oxidation number of carbon. Also, create a matrix having the chemical formulas of amino acid residues in proteins and calculate the chemical formulas of proteins from their amino acid composition.
}

\usage{
  GHS(species = NULL, DG = NA, DH = NA, S = NA, T = thermo$opt$Tr)
  element(compound, property = c("mass","entropy"))
  expand.formula(elements, makeup)
  ZC(x)
  residue.formula()
  protein.formula(proteins, as.residue = FALSE)
}

\arguments{
  \item{species}{character, formula of a compound from which to calculate entropies of the elements.}
  \item{DG}{numeric, standard molal Gibbs energy of formation.}
  \item{DH}{numeric, standard molal enthalpy of formation.}
  \item{S}{numeric, standard molal molal entropy.}
  \item{T}{numeric, temperature in Kelvin.}
  \item{compound}{character, name of element(s) or compound(s).}
  \item{property}{character, name(s) of thermodynamic properties.}
  \item{elements}{character, name(s) of elements.}
  \item{makeup}{dataframe, elemental composition of a compound returned by \code{\link{makeup}}.}
  \item{x}{character, object representing chemical formula.}
  \item{proteins}{dataframe, amino acid composition of one or more proteins in the same format as \code{\link{thermo}$protein}}.
  \item{as.residue}{logical, return the per-residue formula of the protein(s)?}
}

\details{

  \code{GHS} computes one of the standard molal Gibbs energy or enthalpy of formation from the elements (\code{DG}, \code{DH}) or entropy (\code{S}) at 298.15 K and 1 bar from values of the other two. If the \code{species} argument is present, it is used to calculate the entropies of the elements (\code{Se}) using \code{\link{element}}, otherwise \code{Se} is set to zero. The equation in effect can be written as \eqn{{\Delta}G^{\circ}={\Delta}H^{\circ}-T{\Delta}S^{\circ}}{DG = DH - T * DS}, where \eqn{{\Delta}S^{\circ}=S-S_e}{DS = S - Se} and \eqn{T} denotes the reference temperature of 298.15 K. If two of \code{DG}, \code{DH}, and \code{S} are provided, the value of the third is returned. If three are provided, the value of \code{DG} in the arguments is ignored and the calculated value of \code{DG} is returned. If none of \code{DG}, \code{DH} or \code{S} are provided, the value of \code{Se} is returned. If only one of the values is provided, an error results. Units of cal mol\eqn{^{-1}}{^-1} (\code{DG}, \code{DH}) and cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} (\code{S}) are assumed. It \code{T} is provided, it use used instead of the reference temperature.

  \code{element} returns a dataframe of the mass and entropy of one or more elements or formulas given in \code{compound}. The \code{property} can be \samp{mass} and/or \samp{entropy}.

  \code{expand.formula} converts a 1-column dataframe representing the elemental composition of a compound (see \code{\link{makeup}}) to a numeric vector, each value of which is the coefficient of the \code{elements} given in the argument. If any of these is not present in the makeup dataframe, its coefficient is set to zero. A non-zero coefficient of an element in the makeup dataframe does not appear in the output if that element is not one of \code{elements}.

  \code{ZC} returns the nominal carbon oxidation state for the chemical formula represented by \code{x}. (For discussion of nominal carbon oxidation state, see Hendrickson et al., 1970; Buvet, 1983.) If carbon is not present in the formula the result is \code{NaN}.

  \code{protein.formula} exists to quickly compute the chemical formulas of many proteins. The \code{proteins} argument contains the amino acid compositions of the proteins in the same format as the \code{\link{thermo}$protein} dataframe. \code{residue.formula} is called to calculate the chemical formulas of each of the 20 common amino acid residues (and the terminal H- and -OH). The amino acid compositions of the proteins and the output of \code{residue.formula} are multiplied using matrix multiplication to generate the result. 

}

\value{
  \code{GHS} and \code{ZC} return numeric values. \code{expand.formula} returns a numeric vector.
}

\seealso{
  \code{\link{makeup}} can be used to count the elements in formulas and display formulas in various formats.
}

\examples{
  \dontshow{data(thermo)}
  ## converting among Gibbs, enthalpy, entropy
  GHS("H") # entropy of H (element)
  # calculate enthalpy of formation of arsenopyrite 
  GHS("FeAsS",DG=-33843,S=68.5) 
  # return the value of DG calculated from DH and S
  # cf. -56687.71 from subcrt("water")
  GHS("H2O",DH=-68316.76,S=16.7123)  

  ## mass and entropy of compounds of elements
  element("CH4")
  element(c("CH4","H2O"),"mass")
  element("Z")   # charge
  # same mass, opposite entropy as charge
  element("Z-1") # i.e., electron
 
  ## count selected elements in a formula
  m <- makeup("H2O")
  expand.formula(c("H","O"),m)
  expand.formula(c("C","H","S"),m)

  ## calculate the average chemical formula of all of 
  ## the proteins in CHNOSZ' database
  ## this is much faster than a for-loop
  pf <- protein.formula(thermo$protein)
  colSums(pf)/nrow(pf)

  ## nominal carbon oxidation states
  ZC("CO2")  # 4
  ZC("CH4")  # -4
  ZC("CHNOSZ") # 7
  si <- info(info("LYSC_CHICK"))
  ZC(si$formula)  # 0.01631

  ## plot ZC of reference protein sequence
  ## for different organisms
  file <- system.file("extdata/refseq/protein_refseq.csv.xz",package="CHNOSZ")
  ip <- add.protein(file)
  # only use those organisms with a certain
  # number of sequenced bases
  ip <- ip[as.numeric(thermo$protein$abbrv[ip])>100000]
  pf <- protein.formula(thermo$protein[ip,])
  zc <- ZC(pf)
  # the organism names we search for
  # "" matches all organisms
  terms <- c("Streptomyces","Pseudomonas","Salmonella",
    "Escherichia","Vibrio","Bacteroides","Lactobacillus",
    "Staphylococcus","Streptococcus","Methano","Bacillus","Thermo","")
  tps <- thermo$protein$ref[ip]
  plot(0,0,xlim=c(1,13),ylim=c(-0.3,-0.05),pch="",
    ylab="average oxidation state of carbon in proteins",
    xlab="",xaxt="n",mar=c(6,3,1,1))
  for(i in 1:length(terms)) {
    it <- grep(terms[i],tps)
    zct <- zc[it]
    points(jitter(rep(i,length(zct))),zct,pch=20)
  }
  terms[13] <- paste("all organisms")
  axis(1,1:13,terms,las=2)
  title(main=paste("Average Oxidation State of Carbon:",
    "Total Protein per taxID in NCBI RefSeq",sep="\n"))
  
}

\references{

  Buvet, R. (1983) General criteria for the fulfillment of redox reactions, in \emph{Bioelectrochemistry I: Biological Redox Reactions}, Milazzo, G. and Blank, M., eds., Plenum Press, New York, p. 15--50. \url{http://www.worldcat.org/oclc/9282370}
 
  Hendrickson, J. B., Cram, D. J., and Hammond, G. S. (1970) \emph{Organic Chemistry}, 3rd ed., McGraw-Hill, New York, 1279 p. \url{http://www.worldcat.org/oclc/78308}

}

\keyword{util}
