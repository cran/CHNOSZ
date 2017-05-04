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
  \item{to.file}{logical, generate image files for the plots?}
}

\usage{
  examples(do.png = FALSE)
  demos(which = c("sources", "protein.equil", "affinity", "NaCl",
    "density", "ORP", "revisit", "findit", "ionize", "buffer",
    "protbuff", "yeastgfp", "mosaic", "copper", "solubility",
    "wjd", "dehydration", "bugstab", "Shh", "activity_ratios"),
    to.file=FALSE)
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
    \code{protein.equil} \tab chemical activities of two proteins in metastable equilibrium (Dick and Shock, 2011) \cr
    \code{affinity} \tab affinities of metabolic reactions and amino acid synthesis \cr
    \code{NaCl} \tab equilibrium constant for aqueous NaCl dissociation (Shock et al., 1992) \cr
    \code{density} \tab density of \H2O, inverted from IAPWS-95 equations (\code{\link{rho.IAPWS95}}) \cr
    \code{ORP} \tab temperature dependence of oxidation-reduction potential for redox standards \cr
    \code{revisit} \tab coefficient of variation of metastable equilibrium activities of proteins \cr
    \code{findit} \tab minimize the standard deviation of logarithms of activities of sulfur species \cr
    \code{ionize} \tab ionize.aa(): contour plots of net charge and ionization properties of LYSC_CHICK \cr
    \code{buffer} \tab  minerals and aqueous species as buffers of hydrogen fugacity (Schulte and Shock, 1995) \cr
    \code{protbuff} \tab chemical activities buffered by thiol peroxidases or sigma factors \cr
    \code{yeastgfp} \tab subcellular locations: log fO2 - log aH2O and log a - log fO2 diagrams (Dick, 2009) \cr
    \code{mosaic} \tab Eh-pH diagram with two sets of changing basis species (Garrels and Christ, 1965) \cr
    \code{copper} \tab another example of \code{\link{mosaic}}: complexation of Cu with glycine (Aksu and Doyle, 2001) \cr
    \code{solubility} \tab solubility of calcite (cf. Manning et al., 2013) or \CO2 (cf. Stumm and Morgan, 1996) \cr
    \code{wjd} \tab Gibbs energy minimization: prebiological atmospheres and cell periphery of yeast \cr
    \code{dehydration} \tab log K of dehydration reactions; SVG file contains tooltips and links \cr
    \code{bugstab} \tab formation potential of microbial proteins in colorectal cancer (Dick, 2016) \cr
    \code{Shh} \tab affinities of transcription factors relative to Sonic hedgehog \cr
    \code{activity_ratios} \tab mineral stability plots with activity ratios on the axes \cr
  }

}

\references{
  Aksu, S. and Doyle, F. M. (2001) Electrochemistry of copper in aqueous glycine solutions. \emph{J. Electrochem. Soc.} \bold{148}, B51--B57. \url{https://doi.org/10.1149/1.1344532}

  Dick, J. M. (2009) Calculation of the relative metastabilities of proteins in subcellular compartments of \emph{Saccharomyces cerevisiae}. \emph{BMC Syst. Biol.} \bold{3}:75. \url{https://doi.org/10.1186/1752-0509-3-75}

  Dick, J. M. and Shock, E. L. (2011) Calculation of the relative chemical stabilities of proteins as a function of temperature and redox chemistry in a hot spring. \emph{PLoS ONE} \bold{6}, e22782. \url{https://doi.org/10.1371/journal.pone.0022782}

  Dick, J. M. (2016) Proteomic indicators of oxidation and hydration state in colorectal cancer. \emph{PeerJ} \bold{4}:e2238. \url{https://doi.org/10.7717/peerj.2238}

  Garrels, R. M. and Christ, C. L. (1965) \emph{Solutions, Minerals, and Equilibria}, Harper & Row, New York, 450 p. \url{http://www.worldcat.org/oclc/517586}
  
  Manning, C. E., Shock, E. L. and Sverjensky, D. A. (2013) The chemistry of carbon in aqueous fluids at crustal and upper-mantle conditions: Experimental and theoretical constraints. \emph{Rev. Mineral. Geochem.} \bold{75}, 109--148. \url{https://doi.org/10.2138/rmg.2013.75.5}

  Schulte, M. D. and Shock, E. L. (1995) Thermodynamics of Strecker synthesis in hydrothermal systems. \emph{Orig. Life Evol. Biosph.} \bold{25}, 161--173. \url{https://doi.org/10.1007/BF01581580}

  Shock, E. L., Oelkers, E. H., Johnson, J. W., Sverjensky, D. A. and Helgeson, H. C. (1992) Calculation of the thermodynamic properties of aqueous species at high pressures and temperatures: Effective electrostatic radii, dissociation constants and standard partial molal properties to 1000 \eqn{^{\circ}}{°}C and 5 kbar. \emph{J. Chem. Soc. Faraday Trans.} \bold{88}, 803--826. \url{https://doi.org/10.1039/FT9928800803}

  Stumm, W. and Morgan, J. J. (1996) \emph{Aquatic Chemistry: Chemical Equilibria and Rates in Natural Waters}, John Wiley & Sons, New York, 1040 p. \url{http://www.worldcat.org/oclc/31754493}
}


\examples{
\dontshow{data(thermo)}
\dontshow{opar <- par(no.readonly=TRUE)}
demos(c("ORP", "NaCl"))
\dontshow{par(opar)}
\dontrun{
# use the following to run examples in all help topics
examples()}
}


\concept{Extra thermodynamic modeling}
