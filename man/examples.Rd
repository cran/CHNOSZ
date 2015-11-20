\encoding{UTF-8}
\name{examples}
\alias{examples}
\alias{demos}
\title{Run Examples from the Documentation}
\description{
  Run the examples contained in each of the documentation topics. }

\arguments{
  \item{do.png}{logical, generate PNG files for the plots?}
  \item{which}{character, which example to run}
}

\usage{
  examples(do.png = FALSE)
  demos(which = c("sources", "NaCl", "density",
    "nucleobase", "ORP", "revisit", "findit", 
    "ionize", "buffer", "yeastgfp", "mosaic",
    "copper", "solubility", "wjd"), do.png=FALSE)
}

\details{

\code{examples} runs all the examples in the documentation for the package. 
\code{\link{example}} is called for each topic with \code{ask} set to \code{FALSE} (so all of the figures are shown without prompting the user).
If \code{do.png} is TRUE, the plots in the examples are saved as \code{\link{png}} files having names beginning with the name of each of the help topics.

\code{demos} is a function to run other examples that are provided as \code{\link{demo}s}.
\code{\link{demo}} is called with settings to not echo the source code and to not ask before making each plot.
The demo(s) to run is/are specified by \code{which}; the default is to run them in the order of the list below.
See the comments in the source code for more information about each demo.

  \tabular{ll}{
    \code{sources} \tab cross-check the reference list with the thermodynamic database \cr
    \code{NaCl} \tab equilibrium constant for aqueous NaCl dissociation (Shock et al., 1992) \cr
    \code{density} \tab density of H2O, inverted from IAPWS-95 equations (\code{\link{rho.IAPWS95}}) \cr
    \code{nucleobase} \tab relative stabilities of nucleobases and some amino acids \cr
    \code{ORP} \tab oxidation-reduction potential of redox standards as a function of temperature \cr
    \code{revisit} \tab detailed example of usage of \code{\link{revisit}} \cr
    \code{findit} \tab detailed example of usage of \code{\link{findit}} \cr
    \code{ionize} \tab ionize.aa(): contour plots of net charge and ionization properties of LYSC_CHICK \cr
    \code{buffer} \tab ionized proteins as a chemical activity buffer (1. thiol peroxidases 2. sigma factors) \cr
    \code{yeastgfp} \tab logfO2-logaH2O diagrams for model proteins based on YeastGFP localizations (Dick, 2009) \cr
    \code{mosaic} \tab Eh-pH diagram with two sets of changing basis species (Garrels and Christ, 1965) \cr
    \code{copper} \tab another example of \code{\link{mosaic}}: complexation of Cu with glycine (Aksu and Doyle, 2001) \cr
    \code{solubility} \tab solubility of calcite (cf. Manning et al., 2013) or CO2(gas) (cf. Stumm and Morgan, 1996) \cr
    \code{wjd} \tab run.wjd() with proteins: cell periphery of yeast \cr
  }

}

\references{
  Aksu, S. and Doyle, F. M. (2001) Electrochemistry of copper in aqueous glycine solutions. \emph{J. Electrochem. Soc.} \bold{148}, B51--B57. \url{http://dx.doi.org/10.1149/1.1344532}

  Dick, J. M. (2009) Calculation of the relative metastabilities of proteins in subcellular compartments of \emph{Saccharomyces cerevisiae}. \emph{BMC Syst. Biol.} \bold{3}:75. \url{http://dx.doi.org/10.1186/1752-0509-3-75}

  Garrels, R. M. and Christ, C. L. (1965) \emph{Solutions, Minerals, and Equilibria}, Harper & Row, New York, 450 p. \url{http://www.worldcat.org/oclc/517586}
  
  Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \url{http://dx.doi.org/10.2138/rmg.2013.75.5}

  Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \eqn{^{\circ}}{°}C and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \url{http://dx.doi.org/10.1039/FT9928800803}

  Stumm, W. and Morgan, J. J. (1996) \emph{Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters}, John Wiley & Sons, New York, 1040 p. \url{http://www.worldcat.org/oclc/31754493}
}


\examples{
\dontshow{data(thermo)}
demos(c("ORP", "NaCl"))
\dontshow{par(thermo$opar)}
\dontrun{
# use the following to run examples in all help topics
examples()}
}


\concept{Extra thermodynamic modeling}
