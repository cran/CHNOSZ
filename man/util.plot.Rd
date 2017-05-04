\encoding{UTF-8}
\name{util.plot}
\alias{util.plot}
\alias{thermo.plot.new}
\alias{label.plot}
\alias{usrfig}
\alias{label.figure}
\alias{water.lines}
\alias{mtitle}
\alias{ZC.col}
\title{Functions to Create and Modify Plots}

\description{
  Initialize a new plot window using preset parameters, add an axis or title to a plot, generate labels for axes and subplots, add stability lines for water, get colors for a set of numeric values.
}


\usage{
  thermo.plot.new(xlim, ylim, xlab, ylab, cex = par("cex"),
    mar = NULL, lwd = par("lwd"), side = c(1,2,3,4), 
    mgp = c(1.7, 0.3, 0), cex.axis = par("cex"), col = par("col"),
    yline = NULL, axs = "i", do.box = TRUE, las = 1,
    xline = NULL, ...)
  label.plot(x, xfrac = 0.05, yfrac = 0.95, paren = FALSE,
    italic = FALSE, ...)
  usrfig()
  label.figure(x, xfrac = 0.05, yfrac = 0.95, paren = FALSE,
    italic = FALSE, ...)
  water.lines(xaxis = "pH", yaxis = "Eh", T = 298.15, P = "Psat", 
    which = c("oxidation","reduction"), logaH2O = 0, lty = 2, lwd=1,
    col = par("fg"), xpoints = NULL, O2state="gas", plot.it = TRUE)
  mtitle(main, line=0, ...)
  ZC.col(z)
}

\arguments{
  \item{xlim}{numeric, limits of the \eqn{x}{x}-axis}
  \item{ylim}{numeric, limits of the \eqn{y}{y}-axis}
  \item{xlab}{character, \eqn{x}{x}-axis label}
  \item{ylab}{character, \eqn{y}{y}-axis label}
  \item{cex}{numeric, character expansion factor for labels}
  \item{mar}{numeric, width (number of lines) of margins on each side of plot}
  \item{lwd}{numeric, line width}
  \item{side}{numeric, which sides of plot to draw axes}
  \item{mgp}{numeric, sizes of margins of plot}
  \item{cex.axis}{numeric, character expansion factor for names of axes}
  \item{col}{character, line color}
  \item{yline}{numeric, margin line on which to plot \eqn{y}{y}-axis name}
  \item{axs}{character, setting for axis limit calculation}
  \item{do.box}{logical, draw a box around the plot?}
  \item{las}{numeric, style for axis labels}
  \item{xline}{numeric, margin line on which to plot \eqn{x}{x}-axis name}
  \item{...}{further arguments passed to \code{par} or \code{mtext}}
  \item{T}{numeric, temperature (K)}
  \item{x}{character, label to place on plot}
  \item{xfrac}{numeric, fractional location on \eqn{x}{x}-axis for placement of label}
  \item{yfrac}{numeric, fractional location on \eqn{y}{y}-axis for placement of label}
  \item{paren}{logical, add parentheses around label text?}
  \item{italic}{logical, italicize label text?}
  \item{xaxis}{character, description of \eqn{x}{x}-axis}
  \item{yaxis}{character, description of \eqn{y}{y}-axis}
  \item{P}{numeric, pressure (bar)}
  \item{which}{character, which of oxidation/reduction lines to plot}
  \item{logaH2O}{numeric, logarithm of the activity of \eqn{\mathrm{H_2O}}{H2O}}
  \item{lty}{numeric, line type}
  \item{xpoints}{numeric, points to plot on \eqn{x}{x} axis}
  \item{O2state}{character, state of O2}
  \item{plot.it}{logical, plot the lines?}
  \item{main}{character, text for plot title}
  \item{line}{numeric, margin line to place title}
  \item{z}{numeric, set of values}
}

\details{

  \code{thermo.plot.new} sets parameters for a new plot, creates a new plot using \code{\link{plot.new}}, and adds axes and major and minor tick marks to the plot. Plot parameters (see \code{\link{par}}) including \code{cex}, \code{mar}, \code{lwd}, \code{mgp} and \code{axs} can be given, as well as a numeric vector in \code{side} identifying which sides of the plot receive tick marks. \code{yline}, if present, denotes the margin line (default \code{\link{par}('mgp')[1]}) where the y-axis name is plotted.

\code{water.lines} plots lines representing the oxidation and reduction stability limits of water on \code{yaxis}-\code{xaxis} diagrams, where \code{yaxis} can be \samp{Eh} or \samp{O2}, and \code{xaxis} can be \samp{pH} or \samp{T}.
\code{which} controls which lines are drawn (\samp{oxidation}, \samp{reduction}, or both (the default)).
\code{logaH2O} denotes the logarithm of the activity of water.
With \code{O2state} set to \samp{gas} (the default), the logarithm of oxygen fugacity is plotted.
Change this to \samp{aq} to plot the logarithm of oxygen activity (do not change it if plotting \samp{Eh}).
\code{xpoints} is an optional list of points on the x axis to which to restrict the plotting (default of \code{NULL} refers to the axis limits).

\code{label.plot} and \code{label.figure} add identifying text within the plot region and figure region.
The value given for \code{x} is made into a label, optionally italicized and with parentheses (like \ifelse{latex}{\eqn{(a)}}{\ifelse{html}{\out{(<i>a</i>)}}{(a)}}).
The location of the label is controlled by \code{xfrac} and \code{yfrac} (the fractional coordinates of either the plot or figure region), and \code{...} can include other parameters such as \code{cex} and \code{adj} that are passed to \code{\link{text}}.

\code{usrfig} returns the limits of the figure region in \dQuote{user} coordinates (i.e. the limits of the plot region, from \code{\link{par}("usr")}).
It is a supporting function for \code{label.figure} but is also useful for other circumstances where information must be added at a particular location in a figure.

\code{mtitle} can be used to add a multi-line title to a plot.
It loops over each element of \code{main} and places it on a separate margin line using \code{\link{mtext}}.
The spacing of the last (bottom) line from the edge of the plot is specified by \code{line}.
This function exists to facilitate using \code{\link{expression}}s in multiline titles (see \code{\link{revisit}} for an example.) 

\code{ZC.col} uses \CRANpkg{colorspace}) to generate colors from a diverging palette (red - light grey - blue) corresponding to the values in \code{z}.
Red is associated with lower values of \code{z}.
This function is intended to generate colors for distinguishing average oxidation state of carbon \code{\link{ZC}}, but any numeric values can be supplied.

}

\keyword{utilities}
