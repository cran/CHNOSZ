\name{util.plot}
\alias{util.plot}
\alias{thermo.plot.new}
\alias{thermo.axis}
\alias{label.plot}
\alias{water.lines}
\alias{mtitle}
\alias{residualsplot}
\title{Functions to Create and Modify Plots}

\description{
  Initialize a new plot window using preset parameters, open a postcript file for plotting, add an axis or title to a plot, generate labels for plot axes and for identification of subplots and physical and chemical conditions, add stability lines for water to a diagram.
}


\usage{
  thermo.plot.new(xlim, ylim, xlab, ylab, cex = par("cex"),
    mar = NULL, lwd = par("lwd"), side = c(1,2,3,4), 
    mgp = c(1.5, 0.3, 0), cex.axis = par("cex"), col = par("col"),
    yline = NULL, axs = "i", do.box = TRUE, ticks = NULL, las = 1,
    xline = NULL)
  thermo.axis(lab = "x-axis", side = 1, line = 1.5, cex = par("cex"),
    lwd = par("lwd"), T = NULL, col = par("col"))
  label.plot(x, xfrac = 0.95, yfrac = 0.9, cex = 1, paren = TRUE, 
    adj = 1)
  water.lines(xaxis = "pH", yaxis = "Eh", T = 298.15, P = "Psat", 
    which = c("oxidation","reduction"), logaH2O = 0, lty = 2, lwd=1,
    col = par("fg"), xpoints = NULL, O2state="gas")
  mtitle(main, ...)
  residualsplot(residuals, property = "Cp", model = "big")
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
  \item{ticks}{numeric, same effect as \code{side} (retained for backwards compatibility)}
  \item{las}{numeric, style for axis labels}
  \item{xline}{numeric, margin line on which to plot \eqn{x}{x}-axis name}
  \item{lab}{character, description of axis label}
  \item{line}{numeric, margin line to plot axis label}
  \item{T}{numeric, temperature (K)}
  \item{x}{character, label to place on plot}
  \item{xfrac}{numeric, fractional location on \eqn{x}{x}-axis for placement of label}
  \item{yfrac}{numeric, fractional location on \eqn{y}{y}-axis for placement of label}
  \item{paren}{logical, add parentheses around label text?}
  \item{adj}{numeric, parameter for text alignment}
  \item{xaxis}{character, description of \eqn{x}{x}-axis}
  \item{yaxis}{character, description of \eqn{y}{y}-axis}
  \item{P}{numeric, pressure (bar)}
  \item{which}{character, which of oxidation/reduction lines to plot}
  \item{logaH2O}{numeric, logarithm of the activity of \eqn{\mathrm{H_2O}}{H2O}}
  \item{lty}{numeric, line type}
  \item{xpoints}{numeric, points to plot on \eqn{x}{x} axis}
  \item{O2state}{character, state of O2}
  \item{main}{character, text for plot title}
  \item{...}{further arguments passed to \code{mtext}}
  \item{residuals}{numeric, named vector of residuals to plot}
  \item{property}{character, name of property}
  \item{model}{character, name of model to use in plot title}
}

\details{

  \code{thermo.plot.new} sets parameters for a new plot, creates a new plot using \code{\link{plot.new}}, and adds axes and major and minor ticks to the plot. Plot parameters (see \code{\link{par}}) including \code{cex}, \code{mar}, \code{lwd}, \code{mgp} and \code{axs} can be given, as well as a numeric vector in \code{ticks} identifying which sides of the plot receive tick marks. \code{yline}, if present, denotes the margin line (default \code{\link{par}('mgp')[1]}) where the y-axis name is plotted. The very first time \code{thermo.plot.new} is called, it stores the existing value of \code{par(no.readonly=TRUE)} in \code{thermo$opar} so that the plot device can be reset to the previous state later on.

\code{water.lines} plots lines representing the oxidation and reduction stability limits of water on \code{yaxis}-\code{xaxis} diagrams, where \code{yaxis} can be \samp{Eh} or \samp{O2}, and \code{xaxis} can be \samp{pH} or \samp{T}.
\code{which} controls which lines are drawn (\samp{oxidation}, \samp{reduction}, or both (the default)).
\code{logaH2O} denotes the logarithm of the activity of water.
With \code{O2state} set to \samp{gas} (the default), the logarithm of oxygen fugacity is plotted.
Change this to \samp{aq} to plot the logarithm of oxygen activity (do not change it if plotting \samp{Eh}).
\code{xpoints} is an optional list of points on the x axis to which to restrict the plotting (default of \code{NULL} refers to the axis limits).

  \code{label.plot} adds identifying text to the plot; the value given for \code{x} is made into a label like \eqn{(a)}{(a)}. The location of the label is controlled by \code{xfrac} and \code{yfrac} (the fractional locations along the respective axes) as well as \code{adj} (the text alignment parameter, see \code{\link{text}}).

  \code{thermo.axis} is used to add axes and axis labels to plots, with some default style settings (rotation of numeric labels) and conversions between oxidation-reduction scales (called by \code{thermo.plot.new}). It also adds minor tick marks.

  \code{mtitle} can be used to add a multi-line title to a plot. It loops over each element of \code{main} and places it on a separate margin line using \code{\link{mtext}}. This function exists to facilitate using \code{\link{expression}}s in multiline titles (see \code{\link{revisit}} for an example.) 

  \code{residualsplot} produces a \code{\link{barchart}} with options useful for plotting residuals of group contribution models for thermodynamic properties. It plots horizontal bars stacked with largest on top. The names of the \code{residuals} argument (i.e., the names of model species) are plotted across from each respective bar. The axis title is taken from the \code{property} (probably \samp{Cp} or \samp{V}), and the plot title includes the \code{model} name. See the \sQuote{xadditivity} vignette for examples of these plots.

}

\section{Side Effects}{
  These functions create or modify a plot.
}

\keyword{utilities}
