\name{anim.TCA}
\alias{anim.TCA}
\title{Functions to Make Animations}
\description{
  Make animated stability diagrams by creating a series of PNG files.
}

\usage{
  anim.TCA(redox = list(O2 = c(-95, -60)), high.T = FALSE,
    nframes.in = 15, nframes.out = 20, nframes.anim = 140,
    pHlim = c(0,10), width = 768, height = 576)
}

\arguments{
  \item{redox}{list, redox variable and limits.}
  \item{high.T}{logical, overlay high-temperature diagram?}
  \item{nframes.in}{numeric, number of lead-in frames.}
  \item{nframes.out}{numeric, number of lead-out frames.}
  \item{nframes.anim}{numeric, number of frames to be animated.}
  \item{pHlim}{numeric, pH limits to use for animation.}
  \item{width}{numeric, width of plot device.}
  \item{height}{numeric, height of plot device.}
}

\details{
  \code{anim.TCA} creates a series of PNG figures that can be converted into an animated diagram. The animation shows how a logaH2O - logfO2 activity diagram for various species involved in the tricarboxylic acid (TCA) cycle changes as a function of pH. Alternatively, set \code{redox} to \code{list(H2=c(-20,0))} to draw a logaH2O - logaH2 diagram. The diagrams are made at 25 degrees C unless \code{high.T} is TRUE, in which case the high-temperature (100 degrees C) stability fields are also drawn.
  The function fails with an error if either a directory named \samp{png} is not present in the \R working directory, or if the directory exist but is not empty. A number of PNG files are created in the directory. \code{nframes.in} and \code{nframes.out} specify the number of lead-in (before the start of the animation) and lead-out frames (at the end of the animation). The number of frames to be used for the actual animation (as pH increases ranges between the values specified in \code{pHlim}) is given by \code{nframes.anim}.
  After making the PNG files, they are converted to an animated GIF using the \samp{convert} tool from the ImageMagick software distribution (\url{http://www.imagemagick.org}), if it is available on the system.
}

\examples{
  \dontrun{
  # animate the stability diagram with a high-temperature overlay
  anim.TCA(high.T=TRUE)
  # to make an animation with fewer frames than the default
  anim.TCA(nframes.in=5,nframes.out=5,nframes.anim=20)
  # using H2 instead of O2
  anim.TCA(list(H2=c(-20,0)))
  }
}
