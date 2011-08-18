\name{util.plot}
\alias{util.plot}
\alias{thermo.plot.new}
\alias{thermo.postscript}
\alias{thermo.axis}
\alias{label.plot}
\alias{axis.label}
\alias{species.label}
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
    mgp = c(1.2, 0.3, 0), cex.axis = par("cex"), col = par("col"),
    yline = NULL, axs = "i", do.box = TRUE, ticks = NULL, las = 1,
    xline = NULL)
  thermo.postscript(file, family = "Helvetica", width = 8, 
    height = 6, horizontal = FALSE)
  thermo.axis(lab = "x-axis", side = 1, line = 1.5, cex = par("cex"),
    lwd = par("lwd"), T = NULL, col = par("col"))
  label.plot(x, xfrac = 0.95, yfrac = 0.9, cex = 1, paren = TRUE, 
    adj = 1)
  axis.label(lab, opt = NULL, do.state = TRUE, oldstyle = FALSE,
    do.upper = FALSE, mol = "mol", state = NULL, as.expression = TRUE)
  species.label(formula, do.state = FALSE, state = "", do.log = FALSE,
    as.expression = TRUE)
  water.lines(xaxis = "pH", yaxis = "Eh", T = 298.15, P = "Psat", 
    which = c("oxidation","reduction"), logaH2O = 0, lty = 2, 
    col = par("fg"), xpoints = NULL)
  mtitle(main, ...)
  residualsplot(residuals, property = "Cp", model = "big")
}

\arguments{
  \item{xlim}{numeric, limits of the \eqn{x}{x}-axis.}
  \item{ylim}{numeric, limits of the \eqn{y}{y}-axis.}
  \item{xlab}{character, \eqn{x}{x}-axis label.}
  \item{ylab}{character, \eqn{y}{y}-axis label.}
  \item{cex}{numeric, character expansion factor for labels.}
  \item{mar}{numeric, width (number of lines) of margins on each side of plot.}
  \item{lwd}{numeric, line width.}
  \item{side}{numeric, which sides of plot to draw axes.}
  \item{mgp}{numeric, sizes of margins of plot.}
  \item{cex.axis}{numeric, character expansion factor for names of axes.}
  \item{col}{character, line color.}
  \item{yline}{numeric, margin line on which to plot \eqn{y}{y}-axis name.}
  \item{axs}{character, setting for axis limit calculation.}
  \item{do.box}{logical, draw a box around the plot?}
  \item{ticks}{numeric, same effect as \code{side} (retained for backwards compatibility).}
  \item{las}{numeric, style for axis labels}
  \item{xline}{numeric, margin line on which to plot \eqn{x}{x}-axis name.}
  \item{file}{character, path to a file.}
  \item{family}{character, font family.}
  \item{width}{numeric, width of plot.}
  \item{height}{numeric, height of plot.}
  \item{horizontal}{logical, create plot in landscape mode?}
  \item{lab}{character, description of axis label.}
  \item{line}{numeric, margin line to plot axis label.}
  \item{T}{numeric, temperature (K).}
  \item{x}{character, label to place on plot.}
  \item{xfrac}{numeric, fractional location on \eqn{x}{x}-axis for placement of label.}
  \item{yfrac}{numeric, fractional location on \eqn{y}{y}-axis for placement of label.}
  \item{paren}{logical, add parentheses around label text?}
  \item{adj}{numeric, parameter for text alignment.}
  \item{opt}{character or numeric, options for axis labels.}
  \item{do.state}{logical, append state abbreviation to label?}
  \item{oldstyle}{logical, use old style of axis labels?}
  \item{do.upper}{logical, use uppercase letters in axis label?}
  \item{mol}{character, string to use as the denominator of units in axis label.}
  \item{formula}{character, representation of chemical formula.}
  \item{state}{character, designation of physical state.}
  \item{as.expression}{logical, make the return value an expression?}
  \item{do.log}{logical, prepend a text indicating logarithm of activity or fugacity?}
  \item{xaxis}{character, description of \eqn{x}{x}-axis.}
  \item{yaxis}{character, description of \eqn{y}{y}-axis.}
  \item{P}{numeric, pressure (bar).}
  \item{which}{character, which of oxidation/reduction lines to plot.}
  \item{logaH2O}{numeric, logarithm of the activity of \eqn{\mathrm{H_2O}}{H2O}.}
  \item{lty}{numeric, line type.}
  \item{xpoints}{numeric, points to plot on \eqn{x}{x} axis.}
  \item{main}{character, text for plot title.}
  \item{...}{further arguments passed to \code{mtext}.}
  \item{residuals}{numeric, named vector of residuals to plot.}
  \item{property}{character, name of property.}
  \item{model}{character, name of model to use in plot title.}
}

\details{

  \code{thermo.postscript} calls \code{\link{postscript}} with some custom parameters used by the package author (and might be handy for other users of the package as well).

  \code{thermo.plot.new} sets parameters for a new plot, creates a new plot using \code{\link{plot.new}}, and adds axes and major and minor ticks to the plot. Plot parameters (see \code{\link{par}}) including \code{cex}, \code{mar}, \code{lwd}, \code{mgp} and \code{axs} can be given, as well as a numeric vector in \code{ticks} identifying which sides of the plot receive tick marks. \code{yline}, if present, denotes the margin line (default \code{\link{par}('mgp')[1]}) where the y-axis name is plotted. 

  \code{axis.label} returns an \code{\link{expression}} to be used for plotting an axis label, which may be the symbol for a thermodynamic properties, chemical activity or fugacity, or one of \samp{T}, \samp{P}, \samp{Eh}, \samp{pH}, \samp{pe} or \samp{logK}. If \code{as.expression} is \code{FALSE}, the quoted text of the label is returned instead. An expression for chemical activity or fugacity is returned if the first argument is the name of one of the basis species (e.g., \samp{O2}). The expression in this case includes italic and subscripted symbols, unless \code{oldstyle} is \code{TRUE}, when labels with a simpler format (e.g. \samp{O2 (log f)}) are returned. The default value of NULL of \code{opt} means to use the state this basis species is in, or if this basis species is not present to use the value in \code{\link{thermo}$opt$state}. Likewise, if \code{x} is \samp{T} or \samp{P} the units of temperature or pressure are determined using \code{\link{nuts}} (which also refers to \code{\link{thermo}$opt}). \code{do.upper}, if \code{TRUE}, tells the function to print the label using uppercase letters. Labels for properties can be generated by using e.g. \samp{DGf} or \samp{DG0r} as arguments. \code{mol} (default: \samp{mol}) refers to the denominator of the units (default: molality); this can be changed to represent e.g. specific units, by setting \code{mol} to \samp{g}. \code{opt} when generating labels for properties indicates the prefix to place on the units.

  \code{species.label} is like \code{axis.label} but is specifically intended for generating expressions for the chemical formulas of species. The state (in parentheses) is included in the expression only if \code{state} is not NULL. If \code{do.log} is TRUE, the expression will contain a prefix in front of the formula that indicates chemical activity (like \samp{log_a} or \samp{log_f}).

  \code{water.lines} plots lines representing the oxidation and reduction stability limits of water on \code{yaxis}-\code{xaxis} diagrams, where \code{yaxis} can be \samp{Eh} or \samp{O2}, and \code{xaxis} can be \samp{pH} or \samp{T}. \code{which} controls which lines (\samp{oxidation}, \samp{reduction}, or both (the default)) are drawn, \code{logaH2O} (default 0) denotes the logarithm of the activity of water, \code{lty} (default 2) the line type, \code{col} (default \code{\link{par}('fg')}, the foreground color), and \code{xpoints} an optional list of points on the x axis to which to restrict the plotting (default of \code{NULL} refers to the axis limits).

  \code{label.plot} adds identifying text to the plot; the value given for \code{x} is made into a label like \eqn{(a)}{(a)}. The location of the label is controlled by \code{xfrac} and \code{yfrac} (the fractional locations along the respective axes) as well as \code{adj} (the text alignment parameter, see \code{\link{text}}).

  \code{thermo.axis} is used to add axes and axis labels to plots, with some default style settings (rotation of numeric labels) and conversions between oxidation-reduction scales (called by \code{thermo.plot.new}). It also adds minor tick marks.

  \code{mtitle} can be used to add a multi-line title to a plot. It loops over each element of \code{main} and places it on a separate margin line using \code{\link{mtext}}. This function exists to facilitate using \code{\link{expression}}s in multiline titles (see \code{\link{revisit}} for an example.) 

  \code{residualsplot} produces a \code{\link{barchart}} with options useful for plotting residuals of group contribution models for thermodynamic properties. It plots horizontal bars stacked with largest on top. The names of the \code{residuals} argument (i.e., the names of model species) are plotted across from each respective bar. The axis title is taken from the \code{property} (probably \samp{Cp} or \samp{V}), and the plot title includes the \code{model} name. See the \sQuote{xadditivity} vignette for examples of these plots.

}

\section{Side Effects}{
  These functions create or modify a plot.
}

\keyword{util}
