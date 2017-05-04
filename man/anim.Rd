\encoding{UTF-8}
\name{anim.TCA}
\alias{anim}
\alias{anim.TCA}
\alias{anim.carboxylase}
\title{Functions to Make Animations}
\description{
  Make animated stability diagrams by creating a series of PNG files.
}

\usage{
  anim.TCA(redox = list(O2 = c(-95, -60)), high.T = FALSE,
    nframes = 100, pHlim = c(0,10), width = 420, height = 320)
  anim.carboxylase(T = 25:125, ntop = 5, lcex = 0.8, width = 420, height = 320)
}

\arguments{
  \item{redox}{list, redox variable and limits}
  \item{high.T}{logical, overlay high-temperature diagram?}
  \item{nframes}{numeric, number of frames to be animated}
  \item{pHlim}{numeric, pH limits to use for animation}
  \item{width}{numeric, width of plot device}
  \item{height}{numeric, height of plot device}
  \item{T}{numeric, temperature range for animation}
  \item{ntop}{numeric, number of names to show in legend}
  \item{lcex}{numeric, character expansion factor for legend}
}

\details{
These functions create a series of PNG figures that can be converted into an animated diagram.
The PNG files are created in the \samp{png} directory within the current working directory; the functions stop with an error if either this directory is not present or it is present but not empty.
After making the PNG files, they are converted to an animated GIF using the \samp{convert} tool from the ImageMagick software distribution (\url{http://www.imagemagick.org}), if it is available on the system.
The system command is called using \code{\link{system}}, or \code{shell} on Windows platforms.
When installing ImageMagick on Windows, be sure to leave the \sQuote{Add application directory to your system path} option checked; this will make the \samp{convert} command from ImageMagick available in the shell.

To ensure the calculations are performed as intended, each function first removes any existing system definition by calling \code{data(thermo)}.

\code{anim.TCA} creates a series of figures showing how a logaH2O - logfO2 activity diagram for various species involved in the tricarboxylic acid (TCA) cycle changes as a function of pH.
Alternatively, set \code{redox} to \code{list(H2=c(-20,0))} to draw a logaH2O - logaH2 diagram.
The diagrams are made at 25 \eqn{^{\circ}}{°}C unless \code{high.T} is TRUE, in which case high-temperature (100 \eqn{^{\circ}}{°}C) stability fields are overlaid.

\code{anim.carboxylase} animates equilibrium rank-activity diagrams along a combined temperature and logaH2 gradient, or makes a single plot on the default device (without conversion to animated GIF) if a single \code{temperature} is provided.
The proteins in the calculation are 24 carboxylases from a variety of organisms.
There are 12 ribulose phosphate carboxylase and 12 acetyl-coenzyme A carboxylase; 6 of each type are from nominally mesophilic organisms and 6 from nominally thermophilic organisms, shown as blue and red symbols on the diagrams.
The activities of hydrogen at each temperature are calculated using \eqn{\log a_{\mathrm{H_{2}}_{\left(aq\right)}}=-11+3/\left(40\times T\left(^{\circ}C\right)\right)}{logaH2 = -11 + 3/40 * T(degC)}; this equation comes from a model of relative stabilities of proteins in a hot-spring environment (Dick and Shock, 2011).
}

\examples{\dontshow{data(thermo)}
\dontrun{
# animate a stability diagram with a high-temperature overlay
anim.TCA(high.T=TRUE)
# using H2 instead of O2
anim.TCA(list(H2=c(-20,0))) }

# using anim.carboxylase in non-animation mode
anim.carboxylase(T=100)
}

\references{
Dick, J. M. and Shock, E. L. (2011) Calculation of the relative chemical stabilities of proteins as a function of temperature and redox chemistry in a hot spring. \emph{PLoS ONE} \bold{6}, e22782. \url{https://doi.org/10.1371/journal.pone.0022782}
}

\concept{Extra thermodynamic modeling}
