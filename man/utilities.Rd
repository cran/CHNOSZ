\name{utilities}
\alias{utilities}
\alias{c2s}
\alias{s2c}
\alias{TP.args}
\alias{state.args}
\alias{eos.args}
\alias{GHS}
\alias{can.be.numeric}
\alias{expand.formula}
\alias{protein.length}
\alias{aminoacids}
\alias{nucleicacids}
\alias{thermo.plot.new}
\alias{thermo.postscript}
\alias{axis.label}
\alias{water.lines}
\alias{dPdTtr}
\alias{Ttr}
\alias{ZC}
\alias{label.plot}
\alias{thermo.axis}
\alias{describe}
\alias{basis.comp}
\alias{MP90.cp}
\alias{element}
\alias{mod.obigt}
\alias{which.pmax}
\alias{nonideal}
\alias{change}
\alias{add.protein}
\alias{add.obigt}
\alias{which.balance}
\alias{lsub}
\alias{lsum}
\alias{pprod}
\alias{psum}
\alias{mylapply}
\alias{spearman}
\alias{unitize}
\alias{grep.file}
\alias{read.fasta}
\alias{mtitle}
\alias{examples}
\alias{longex}
\title{Utility and Miscellaneous Functions}
\description{
  Provide various utilities for the user and for other functions in \pkg{CHNOSZ}. Convert between strings and character objects, calculate one of Gibbs energy, enthalpy or entropy from the other two, test for ability to become numeric, write and extract parts of chemical formulas and calculate nominal carbon oxidation states of formulas, handle arguments referring to temperature, pressure, states, and equations of state, calculate protein length, count amino acids in protein sequences, calculate \eqn{dP/dT}{dP/dT} and temperature of phase transitions, calculate heat capacities of unfolded proteins using an equation from the literature, initialize a new plot window using preset parameters, open a postcript file for plotting, add an axis or title to a plot, generate labels for plot axes and for identification of subplots and physical and chemical conditions, add stability lines for water to a diagram, add or alter properties of species in the thermodynamic database, calculate non-ideal contributions to apparent standard molal properties, identify a conserved basis species, perform arithmetic on lists, execute multicore calculations, calculate Spearman's rank correlation coefficient, scale logarithms of activity to a desired total activity, search a FASTA file for protein sequences, read protein sequences from a file and count numbers of amino acids in each sequence, and run all the examples provided in the package.
}

\usage{
  c2s(x, sep = " ")
  s2c(x, sep = NULL, keep.sep = TRUE, n = NULL, move.sep = FALSE)
  GHS(species = NULL, DG = NA, DH = NA, S = NA, T = thermo$opt$Tr)
  can.be.numeric(x)
  expand.formula(elements, makeup)
  ZC(x)
  eos.args(eos, property = NULL, T = NULL, P = NULL)
  TP.args(T = NULL, P = NULL)
  state.args(state = NULL)
  protein.length(protein)
  aminoacids(seq, nchar=1)
  nucleicacids(seq, type = "DNA", comp = NULL, comp2 = NULL)
  MP90.cp(T, protein)
  dPdTtr(x)
  Ttr(x, P = 1, dPdT = NULL)
  thermo.plot.new(xlim, ylim, xlab, ylab, cex = par("cex"),
    mar = NULL, lwd = par("lwd"), side = c(1,2,3,4), 
    mgp = c(1.2, 0.3, 0), cex.axis = par("cex"), col = par("col"),
    yline = NULL, axs = "i", do.box = TRUE, ticks = NULL)
  thermo.postscript(file, family = "Helvetica", width = 8, 
    height = 6, horizontal = FALSE)
  thermo.axis(lab = "x-axis", side = 1, line = 1.5, cex = par("cex"),
    lwd = par("lwd"), T = NULL, col = par("col"))
  axis.label(x, opt = NULL, do.state = TRUE, oldstyle = FALSE,
    do.upper = FALSE, mol = "mol")
  describe(x = NULL, T = NULL, P = NULL, use.name = FALSE, 
    as.reaction = NULL, digits = 1)
  basis.comp(basis)
  label.plot(x, xfrac = 0.95, yfrac = 0.9, cex = 1, paren = TRUE, 
    adj = 1)
  water.lines(xaxis = "pH", yaxis = "Eh", T = 298.15, P = "Psat", 
    which = c("oxidation","reduction"), logaH2O = 0, lty = 2, 
    col = par("fg"), xpoints = NULL)
  element(compound, property = c("mass","entropy"))
  mod.obigt(species, ..., missingvalues = NA)
  which.pmax(elts, na.rm = FALSE, pmin = FALSE)
  nonideal(species, proptable, IS, T)
  change(name, ...)
  add.protein(file="protein.csv")
  add.obigt(file="obigt.csv")
  which.balance(species)
  lsub(x, y)
  lsum(x, y)
  pprod(x, y)
  psum(x)
  mylapply(X, FUN, ...)
  spearman(a, b)
  unitize(logact = NULL, length = NULL, logact.tot = 0)
  mtitle(main, ...)
  grep.file(file, x = "", y = NULL, ignore.case = TRUE, startswith = ">", 
    lines = NULL, grep = "grep")
  read.fasta(file, i = NULL, ret = "count", lines = NULL, ihead = NULL)
  examples(do.png = FALSE, long = TRUE)
  longex(which)
}

\arguments{
  \item{x}{character object to convert (\code{s2c}, \code{c2s}, \code{axis.label}), or object to be tested (\code{can.be.numeric}), or numeric index of a mineral phase (\code{dPdTtr}, \code{Ttr}), or character object representing chemical formula (\code{ZC}), or character text to be added to plot (\code{plot.label}), or dataframe, definition of basis species or reaction (\code{describe}), or list (\code{lsub}, \code{lsum}, \code{psum}), or string to search in sequence headers (\code{grep.file}).}
  \item{sep}{character, the separator to insert or separator(s) to match (\code{c2s}, \code{s2c}).}
  \item{keep.sep}{logical, retain the separator in the output (TRUE) or discard it (FALSE) (\code{s2c}).}
  \item{n}{numeric, maximum number of items in the character object returned by \code{s2c}.}
  \item{move.sep}{logical, move the kept separator to the end of the preceding item.}
  \item{species}{character, formula of a compound from which to calculate entropies of the elements (\code{GHS}), or names of species to modify or add to the thermodynamic database (\code{mod.obigt}), or names or indices of species for which to calculate nonideal properties (\code{nonideal}), or dataframe, species definition such as that in \code{\link{thermo}$species} (\code{which.balance}).}
  \item{T}{numeric, temperature (K) (\code{TP.args}, \code{lines.water}, \code{describe}, \code{MP90.cp}, \code{nonideal}, \code{GHS}).}
  \item{P}{numeric, pressure (bar) (can also be character, \samp{Psat} in \code{TP.args}).}
  \item{eos}{character, name of equation of state (one of \samp{hkf}, \samp{cgl}, \samp{water}).}
  \item{property}{character, name(s) of thermodynamic properties (\code{eos.args}, \code{element}).}
  \item{state}{character, name(s) of states (e.g., \samp{cr}, \samp{aq}).}
  \item{DG}{numeric, standard molal Gibbs energy of formation (\code{GHS}).}
  \item{DH}{numeric, standard molal enthalpy of formation.}
  \item{S}{numeric, standard molal molal entropy.}
  \item{elements}{character, name(s) of elements (\code{expand.formula}).}
  \item{makeup}{dataframe, elemental composition of a compound returned by \code{\link{makeup}}.}
  \item{protein}{character, name of protein species; numeric, species index of protein (\code{protein.length}).}
  \item{seq}{character, amino acid sequence of a protein (\code{aminoacids}) or base sequence of a nucleic acid (\code{nucleicacids}).}
  \item{nchar}{numeric, \eqn{1} to return one-letter, \eqn{3} to return three-letter abbreviations for amino acids (\code{aminoacids}).}
  \item{type}{character, type of nucleic acid sequence (DNA or RNA) (\code{nucleicads}).}
  \item{comp}{character, type of complement sequence.}
  \item{comp2}{character, type of second complement sequence.}
  \item{dPdT}{numeric, values of (\eqn{dP/dT}{dP/dT}) of phase transitions (\code{Ttr}).}
  \item{xlim}{numeric, limits of the \eqn{x}{x}-axis (\code{thermo.plot.new}).}
  \item{ylim}{numeric, limits of the \eqn{y}{y}-axis.}
  \item{xlab}{character, \eqn{x}{x}-axis label.}
  \item{ylab}{character, \eqn{y}{y}-axis label.}
  \item{cex}{numeric, character expansion factor for labels, also in \code{plot.label}.}
  \item{mar}{numeric, width (number of lines) of margins on each side of plot.}
  \item{lwd}{numeric, line width.}
  \item{side}{numeric, which sides of plot to draw axes (\code{thermo.plot.new},\code{thermo.axis}).}
  \item{mgp}{numeric, sizes of margins of plot.}
  \item{cex.axis}{numeric, character expansion factor for names of axes.}
  \item{yline}{numeric, margin line on which to plot \eqn{y}{y}-axis name.}
  \item{axs}{character, setting for axis limit calculation.}
  \item{do.box}{logical, draw a box around the plot?}
  \item{ticks}{numeric, same effect as \code{side} (retained for back-compatibility).}
  \item{file}{character, path to a file (\code{thermo.postscript}, \code{add.protein}, \code{add.obigt}, \code{grep.file}, \code{read.fasta}).}
  \item{family}{character, font family.}
  \item{width}{numeric, width of plot.}
  \item{height}{numeric, height of plot.}
  \item{horizontal}{logical, create plot in landscape mode?}
  \item{opt}{character or numeric, options for axis labels (\code{axis.label}).}
  \item{oldstyle}{logical, use old style of axis labels?}
  \item{do.state}{logical, append state abbreviation to label?}
  \item{do.upper}{logical, use uppercase letters in axis label?}
  \item{mol}{character, string to use as the denominator of units in axis label.}
  \item{use.name}{logical, write names instead of formulas? (\code{describe}).}
  \item{as.reaction}{logical, interpret input as a reaction?}
  \item{digits}{numeric, how many digits to round logarithms of activities.}
  \item{basis}{numeric or character, species number or formula (\code{basis.comp}).}
  \item{xaxis}{character, description of \eqn{x}{x}-axis (\code{water.lines}).}
  \item{yaxis}{character, description of \eqn{y}{y}-axis.}
  \item{which}{character, which of oxidation/reduction lines to plot (\code{water.lines}), or which examples to run (\code{longex}).}
  \item{logaH2O}{numeric, logarithm of the activity of \eqn{\mathrm{H_2O}}{H2O}.}
  \item{lty}{numeric, line type.}
  \item{col}{character, line color (\code{water.lines}, \code{thermo.plot.new}, \code{thermo.axis}).}
  \item{xpoints}{numeric, points to plot on \eqn{x}{x} axis.}
  \item{xfrac}{numeric, fractional location on \eqn{x}{x}-axis for placement of label (\code{label.plot}).}
  \item{yfrac}{numeric, fractional location on \eqn{y}{y}-axis for placement of label.}
  \item{paren}{logical, add parentheses around label text?}
  \item{adj}{numeric, parameter for text alignment.}
  \item{lab}{character, axis label (\code{thermo.axis}).}
  \item{line}{numeric, line (distance from axis) to place axis label.}
  \item{compound}{character, name of element(s) or compound(s) (\code{element}).}
  \item{...}{character or numeric, properties of species to modify in the thermodynamic database (\code{mod.obigt}), or arguments to \code{change} that are passed to \code{mod.obigt} or \code{mod.buffer}, or additional arguments for \code{lappaly} or \code{mclapply} (\code{mylapply}), or arguments passed to \code{mtext} (\code{mtitle}).}
  \item{missingvalues}{numeric, values to assign to undefined properties.}
  \item{elts}{list, numeric vectors for which to find maximum values (in parallel) (\code{which.pmax}).}
  \item{na.rm}{logical, remove missing values?}
  \item{pmin}{logical, find minimum values instead of maximum ones?}
  \item{proptable}{list of dataframes of species properties (\code{nonideal}).}
  \item{IS}{numeric, ionic strength(s) used in nonideal calculations, mol kg\eqn{^{-1}}{^-1}.}
  \item{name}{character or numeric, name (or numeric index) of species or name of buffer to be modified (\code{change}).}
  \item{X}{vector, argument for \code{\link{lapply}} or \code{mclapply}.}
  \item{FUN}{function, argument for \code{lapply} or \code{mclapply}.}
  \item{y}{list (\code{lsub}, \code{lsum}), or numeric (\code{pprod}), or character, term to exclude in searching sequence headers (\code{grep.file}).}
  \item{ignore.case}{logical, ignore differences between upper- and lower-case? (\code{grep.file})}
  \item{startswith}{character, only lines starting with this expression are matched.}
  \item{lines}{list of character, supply the lines here instead of reading them from file.}
  \item{grep}{character, name of system \code{grep} command.}
  \item{i}{numeric, line numbers of sequence headers to read (\code{read.fasta}).}
  \item{ret}{character, specification for type of return (count, sequence, or FASTA format).}
  \item{ihead}{numeric, which lines are headers.}
  \item{a}{numeric values (\code{spearman}).}
  \item{b}{numeric values (\code{spearman}).}
  \item{logact}{numeric, logarithms of activity (\code{unitize}).}
  \item{length}{numeric, numbers of residues.}
  \item{logact.tot}{numeric, logarithm of total activity.}
  \item{main}{character, text for plot title.}
  \item{do.png}{logical, generate PNG plot files?}
  \item{long}{logical, run additional examples?}
  
}

\details{

  \code{c2s} joins the elements of a character object into a character object of length \eqn{1} (a string), and \code{s2c} splits a string into elements of a character object of length \eqn{n+1}, where \eqn{n} stands for the number of separators in the string. \code{sep} gives the separator to insert between successive items (in \code{c2s}) or the separator(s) to find in a string (in \code{s2c}). The default value of \code{sep} is a space (" ") in \code{c2s}. The default value for \code{sep} is NULL in \code{s2c}, indicating a separator at every position of \code{x} (the result in this case has length equal to \code{nchar(x)}). Argument \code{keep.sep} if \code{TRUE} (the default) instructs \code{s2c} to keep the separating values in the output, and \code{move.sep} if \code{TRUE} instructs \code{s2c} to append the kept separators to the preceding items instead of prepending them to the following ones. The maximum length of the object returned by \code{s2c} is determined by \code{n}; the default value of \code{NULL} indicates an unrestricted length.

  The \code{*.args} functions are used to normalize user-input arguments, which are case-insensitive. \code{eos.args} returns a list with elements named \code{props}, for all the properties available for the specified equations-of-state, \code{prop} for the lower-case version of \code{property}, and \code{Prop}, for the upper-case (of first letter) version of \code{property}. \code{eos.args} produces an error if one of the \code{property}s is not in the list of available properties. (See \code{\link{water}} and \code{\link{subcrt}} for the available properties for different species.) \code{TP.args} forces \code{T} and \code{P} to equal length. This function also looks for the keyword \samp{Psat} in the value of \code{P} and substitutes calculated values of the saturation vapor pressure (see \code{\link{water}}). \code{state.args} makes its argument lowercase, then transforms \samp{a}, \samp{c}, \samp{g}, and \samp{l} to \samp{aq}, \samp{gas}, \samp{cr}, and \samp{liq}, respectively.

  \code{GHS} computes one of the standard molal Gibbs energy or enthalpy of formation from the elements (\code{DG}, \code{DH}) or entropy (\code{S}) at 298.15 K and 1 bar from values of the other two. If the \code{species} argument is present, it is used to calculate the entropies of the elements (\code{Se}) using \code{\link{element}}, otherwise \code{Se} is set to zero. The equation in effect can be written as \eqn{{\Delta}G^{\circ}={\Delta}H^{\circ}-T{\Delta}S^{\circ}}{DG = DH - T * DS}, where \eqn{{\Delta}S^{\circ}=S-S_e}{DS = S - Se} and \eqn{T} denotes the reference temperature of 298.15 K. If two of \code{DG}, \code{DH}, and \code{S} are provided, the value of the third is returned. If three are provided, the value of \code{DG} in the arguments is ignored and the calculated value of \code{DG} is returned. If none of \code{DG}, \code{DH} or \code{S} are provided, the value of \code{Se} is returned. If only one of the values is provided, an error results. Units of cal mol\eqn{^{-1}}{^-1} (\code{DG}, \code{DH}) and cal K\eqn{^{-1}}{^-1} mol\eqn{^{-1}}{^-1} (\code{S}) are assumed. It \code{T} is provided, it use used instead of the reference temperature.

  \code{can.be.numeric} returns a value of TRUE or FALSE for each element of \code{x}.

  \code{expand.formula} converts a 1-column dataframe representing the elemental composition of a compound (see \code{\link{makeup}}) to a numeric vector, each value of which is the coefficient of the \code{elements} given in the argument. If any of these is not present in the makeup dataframe, its coefficient is set to zero. A non-zero coefficient of an element in the makeup dataframe does not appear in the output if that element is not one of \code{elements}.

  \code{ZC} returns the nominal carbon oxidation state for the chemical formula represented by \code{x}. (For discussion of nominal carbon oxidation state, see Hendrickson et al., 1970; Buvet, 1983.) If carbon is not present in the formula the result is \code{NaN}.

  The argument of \code{protein.length}, if it is character, refers to the name of protein(s) (e.g., \samp{LYSC_CHICK}) for which to calculate the length (number of amino acid residues). If the argument is numeric, it refers to the species index of a protein (value in \code{thermo$species$ispecies}). For a positive numeric argument to work, the protein information must have been previously loaded into the species list (using \code{\link{info}}). If the numeric value is negative, it refers to the rownumber of the protein in \code{thermo$protein}.
  
  \code{aminoacids} takes a character argument containing a protein sequence and counts the number of occurrences of each type of amino acid. The output is a dataframe with 20 columns, each corresponding to an amino acid, ordered in the same way as \code{thermo$protein}. If the first argument is NULL, the function returns the one-letter abbreviations (for \code{nchar} equal to 1) or the three-letter ones (if \code{nchar} is equal to 3) or the names of the amino acids (if \code{nchar} is NA) of twenty amino acids in the order used in \code{thermo$protein}.

  \code{nucleicacids} takes a DNA or RNA sequence and counts the numbers of bases of each type. Whether the sequence is DNA or RNA is specified by \code{type}. Setting \code{comp} to \samp{DNA} or \samp{RNA} tells the function to compute the base composition of that type of complement of the sequence. If \code{comp2} is specified, another complement is taken. The two rounds of complementing can be used in a single function call e.g. to go from a sequence on DNA minus strand (given in \code{seq}) to the plus strand (with \code{comp="DNA"}) and then from the DNA plus strand to RNA (with \code{comp2="RNA"}). The value returned by the function is a dataframe of base composition, which can be passed back to the function to obtain the overall chemical formula for the bases.

  \code{MP90.cp} takes \code{T} (one or more temperatures in \eqn{^{\circ}}{degrees }C) and \code{protein} (name of protein) and returns the heat capacity of the unfolded protein using values of heat capacities of the residues taken from Makhatadze and Privalov, 1990. Those authors provided values of heat capacity at six points between 5 and 125 \eqn{^{\circ}}{degrees }C; this function interpolates (using \code{\link{splinefun}}) values at other temperatures.

  \code{dPdTtr} returns values of \eqn{(dP/dT)_{Ttr}}{(dP/dT)Ttr}, where \eqn{Ttr}{Ttr} represents the transition temperature, of the phase transition at the high-\eqn{T}{T} stability limit of the \code{x}th species in \code{thermo$obigt} (no checking is done to verify that the species represents in fact one phase of a mineral with phase transitions). \code{dPdTtr} takes account of the Clapeyron equation, \eqn{(dP/dT)_{Ttr}}{(dP/dT)Ttr}=\eqn{{\Delta}S/{\Delta}V}{DS/DV}, where \eqn{{\Delta}S}{DS} and \eqn{{\Delta}V}{DV} represent the changes in entropy and volume of phase transition, and are calculated using \code{subcrt} at Ttr from the standard molal entropies and volumes of the two phases involved. Using values of \code{dPdT} calculated using \code{dPdTtr} or supplied in the arguments, \code{Ttr} returns as a function of \code{P} values of the upper transition temperature of the mineral phase represented by the \code{x}th species.

  \code{thermo.plot.new} sets parameters for a new plot, creates a new plot using \code{\link{plot.new}}, and adds axes and major and minor ticks to the plot. Plot parameters (see \code{\link{par}}) including \code{cex}, \code{mar}, \code{lwd}, \code{mgp} and \code{axs} can be given, as well as a numeric vector in \code{ticks} identifying which sides of the plot receive tick marks. \code{yline}, if present, denotes the margin line (default \code{\link{par}('mgp')[1]}) where the y-axis name is plotted. \code{thermo.postscript} calls \code{\link{postscript}} with some custom parameters. 

  \code{axis.label} returns an \code{\link{expression}} to be used for plotting an axis label, which may be the symbol for a thermodynamic properties, chemical activity or fugacity, or one of \samp{T}, \samp{P}, \samp{Eh}, \samp{pH}, \samp{pe} or \samp{logK}. An expression for chemical activity or fugacity is returned if the first argument is the name of one of the basis species (e.g., \samp{O2}). The expression in this case includes italic and subscripted symbols, unless \code{oldstyle} is \code{TRUE}, when labels with a simpler format (e.g. \samp{O2 (log f)}) are returned. The default value of NULL of \code{opt} means to use the state this basis species is in, or if this basis species is not present to use the value in \code{\link{thermo}$opt$state}. Likewise, if \code{x} is \samp{T} or \samp{P} the units of temperature or pressure are determined using \code{\link{nuts}} (which also refers to \code{\link{thermo}$opt}). \code{do.upper}, if \code{TRUE}, tells the function to print the label using uppercase letters. Labels for properties can be generated by using e.g. \samp{DGf} or \samp{DG0r} as arguments. \code{mol} (default: \samp{mol}) refers to the denominator of the units (default: molality); this can be changed to represent e.g. specific units, by setting \code{mol} to \samp{g}. \code{opt} when generating labels for properties indicates the prefix to place on the units.

  \code{water.lines} plots lines representing the oxidation and reduction stability limits of water on \code{yaxis}-\code{xaxis} diagrams, where \code{yaxis} can be \samp{Eh} or \samp{O2}, and \code{xaxis} can be \samp{pH} or \samp{T}. \code{which} controls which lines (\samp{oxidation}, \samp{reduction}, or both (the default)) are drawn, \code{logaH2O} (default 0) denotes the logarithm of the activity of water, \code{lty} (default 2) the line type, \code{col} (default \code{\link{par}('fg')}, the foreground color), and \code{xpoints} an optional list of points on the x axis to which to restrict the plotting (default of \code{NULL} refers to the axis limits).

  \code{label.plot} adds identifying text to the plot; the value given for \code{x} is made into a label like \eqn{(a)}{(a)}. The location of the label is controlled by \code{xfrac} and \code{yfrac} (the fractional locations along the respective axes) as well as \code{adj} (the text alignment parameter, see \code{\link{text}}).

  \code{thermo.axis} is used to add axes and axis labels to plots, with some default style settings (rotation of numeric labels) and conversions between oxidation-reduction scales (called by \code{thermo.plot.new}). It also adds minor tick marks.

  \code{describe} generates a textual representation of the temperature, pressure, and logarithms of activities of the basis species, given in the arguments by \code{x} (i.e. the dataframe in \code{thermo$basis}) and \code{T} and \code{P} (given in Kelvin and bar and converted by the function to those specified by \code{\link{nuts}}). The \code{digits} argument tells to what decimal place the logarithms of activities should be rounded. If any of the supplied arguments is NULL its specification is not printed in the output; T and P, if present, are prepended to the basis summary. If \code{x} instead is a dataframe representing a chemical reaction (as output by \code{\link{subcrt}} and identified by having a column named \code{coeff}), the function returns a textual summary of that reaction (i.e., showing reactants on the left, an equal sign, and products on the right; reactants and products are preceded by their reaction coefficient unless it is \eqn{1}). However, if only reactants or products can be found, or \code{as.reaction} is set to FALSE, the names or formulas of the species are printed with their coefficients and interceding plus or minus signs, as appropriate. Whether the names or formulas are printed is controlled by \code{use.name} (\code{FALSE} by default), a logical vector the length of which corresponds to the number of rows in \code{x} (but is expanded to the right length if needed). 

  \code{element} returns a dataframe of the mass and entropy of one or more elements or formulas given in \code{compound}. The \code{property} can be \samp{mass} and/or \samp{entropy}.
 
  \code{mod.obigt} changes one or more of the properties of one or more species or adds species to the thermodynamic database. These changes are lost if you reload the database by calling \code{\link{data}(thermo)} or if you quit the \R session without saving it. To modify the properties of species, give the names in the \code{species} argument and supply other arguments: if one of these arguments is \code{state}, species in those states will be updated. Additional arguments refer to the name of the property(s) to be updated and correspond to the column names of \code{\link{thermo}$obigt} (the names of the properties are matched to any part of compound column names, such as \samp{z.T}). The values provided should be in the units specifed in the documentation for the \code{thermo} data object. To add species, supply the new names in \code{species} and provide an argument named \code{formula} with the corresponding chemical formulas. Additional arguments refer to any of the properties you wish to specify. Properties that are not specified are assigned the value of \code{missingvalues} which is \code{NA} by default (however if \code{state} is missing it is set to the value of \code{thermo$opt$state}). The values returned (\code{\link{invisible}}-y) by \code{mod.obigt} are the rownumbers of the affected species.

  \code{which.pmax} takes a list of equal-length numeric vectors (or objects that can be coerced to numeric) in \code{elts} and returns the index of the vector holding the maximum value at each position. If \code{na.rm} is \code{TRUE}, values of \code{NA} are removed; if \code{pmin} is \code{TRUE} the function finds locations of the minimum values instead.

  \code{nonideal} takes a list of dataframes (in \code{proptable}) containing the standard molal properties of the identified \code{species}. For those species whose charge (determined by the number of Z in their \code{\link{makeup}}) is not equal to zero, the values of \code{IS} are combined with Alberty's (2003) equation 3.6-1 (Debye-Huckel equation) and its derivatives, to calculate apparent molal properties at the specified ionic strength(s) and temperature(s). The lengths of \code{IS} and \code{T} supplied in the arguments should be equal to the number of rows of each dataframe in \code{proptable}, or one to use single values throughout. The apparent molal properties that can be calculated include \code{G}, \code{H}, \code{S} and \code{Cp}; any columns in the dataframes of \code{proptable} with other names are left untouched. A column named \code{loggam} (logarithm of gamma, the activity coefficient) is appended to the output dataframe of species properties.

  \code{change} is a wrapper function to \code{\link{mod.obigt}} and \code{\link{mod.buffer}}. The \code{name} provided in the argument refers to the name or numeric index of the species to update or add using \code{\link{mod.obigt}}, unless the \code{name} begins with an underscore character, in which case the remaining part of the name (after the underscore) is passed to \code{\link{mod.buffer}}. The arguments in \code{...} are sent without change to the subordinate function.

  \code{add.protein} and \code{add.obigt} read data from the specified \code{file} and add it to either \code{\link{thermo}$protein} or \code{\link{thermo}$obigt}, as appropriate. Both of these are functions are run, with the default file names, when \pkg{CHNOSZ} is first loaded.


  \code{which.balance} returns, in order, which column(s) of \code{species} all have non-zero values. It is used by \code{\link{diagram}} and \code{\link{transfer}} to determine a conservant (i.e. basis species that are conserved in transformation reactions) if none is supplied by the user.

  \code{lsub} subtracts the elements of list \code{y} from the respective ones in list \code{x}. \code{lsum} sums the respective elements of lists \code{x} and \code{y}. \code{pprod} multiplies each element of list \code{x} by the respective numeric value in \code{y}. \code{psum} sums all elements of the list \code{x}.

 \code{mylapply} passes the given arguments to \code{\link{lapply}}, or to \code{mclapply} if the \pkg{multicore} package is loaded and the length of \code{X} is greater than 20. \code{mylapply} is used in \code{\link{affinity}} (in calculations for proteins activated by the \code{iprotein} argument), \code{\link{abundance.new}} (in parallel operations on list elements), and \code{\link{aminoacids}} and \code{\link{protein.length}} (in counting amino acids in sequences and determining lengths of proteins).

  \code{grep.file} is used to search for entries in a FASTA file. It returns the line numbers of the matching FASTA headers. It takes a search term in \code{x} and optionally a term to exclude in \code{y}. The \code{ignore.case} option is passed to \code{\link{grep}}, which does the work of finding lines that match. Only lines that start with the expression in \code{startswith} are searched; the default setting reflects the format of the header line for each sequence in a FASTA file. If \code{y} is NULL and a supported operating system is identified (right now, only Linux), the operating system's \samp{grep} function (or other specified in the \code{grep} argument) is applied directly to the file instead of R's \code{\link{grep}}; this mitigates the potential speed penalty of having to read the file into R using \code{\link{readLines}}. If the lines from the file were obtained in a preceding operation, they can be supplied to this function in the \code{lines} argument (this overrides the use of the OS's \samp{egrep}).

  \code{read.fasta} is used to retrieve entries from a FASTA file. The line numbers for the headers of the desired sequences are passed to the function in \code{i} (they can be generated using \code{grep.file}). The function returns various formats depending on the value of \code{ret}; the default \samp{count} returns a dataframe of amino acid counts (the dataframe can be given to \code{\link{add.protein}} in order to add the proteins to \code{\link{thermo}$protein}), \samp{seq} returns a list of sequences, and \samp{fas} returns a list of lines extracted from the FASTA file, including the headers (this can be used e.g. to generate a new FASTA file with only the selected sequences). Similarly to \code{grep.file}, this function utilizes the OS's \samp{grep} on supported operating systems in order to identify the header lines as well as \samp{cat} to read the file, otherwise \code{\link{readLines}} and \R's \code{\link{substr}} are used to read the file and locate the header lines. \code{lines}, if it is given, bypasses the reading of the file and also overrides the use of the OS's tools. If the line numbers of the header lines were previously determined, they can be supplied in \code{ihead}. When computing relative abundances of many proteins that might be found with \code{grep.fasta}, consider using the \code{iprotein} arugment of \code{\link{affinity}} to speed things up. Examples of these operations can be found in the documentation for \code{\link{revisit}}.

  \code{spearman} calculates Spearman's rank correlation coefficient for \code{a} and \code{b}.

  \code{unitize} scales the logarithms of activities given in \code{logact} so that the logarithm of total activity of residues is equal to zero (i.e. total activity of residues is one), or to some other value set in \code{logact.tot}.  \code{length} indicates the number of residues in each species. If \code{logact} is NULL, the function takes the logarithms of activities from the current species definition. If any of those species are proteins, the function gets their lengths using \code{protein.length}.

  \code{mtitle} can be used to add a multi-line title to a plot. It loops over each element of \code{main} and places it on a separate margin line using \code{\link{mtext}}. This function exists to facilitate using \code{\link{expression}}s in multiline titles (see \code{\link{revisit}} for an example.) 

  \code{examples} runs all the examples in the package using \code{\link{example}} (with \code{ask} set to \code{FALSE}). If \code{do.png} is TRUE, the plots in the examples are saved as \code{\link{png}} files having names beginning with the name of each of the help topics. If \code{long} is TRUE (the default), additional examples are run using \code{longex}. \code{longex} contains the text of many of the \samp{dontrun} examples in the documentation, which are marked as such in order to avoid long R CMD check timings.

}

\value{

  \code{s2c}, \code{c2s} and \code{axis.label} return character values. Numeric returns are made by \code{GHS}, \code{protein.length}, \code{dPdTtr}, \code{Ttr}, \code{ZC}, \code{MP90.cp}, \code{spearman}, \code{mod.obigt} and \code{grep.file}. A list is return by \code{eos.args} and \code{TP.args}, and character is returned by \code{state.args}. \code{can.be.numeric} returns logical. \code{aminoacids} and \code{nucleicacids} return character or dataframe. \code{lsub}, \code{lsum} and \code{pprod} return lists. \code{read.fasta} returns a list of sequences or lines (for \code{ret} equal to \samp{seq} or \samp{fas}, respectively), or a data frame with amino acid compositions of proteins (for \code{ret} equal to \samp{count}) with columns corresponding to those in \code{\link{thermo}$protein}. Functions with no (or unspecified) returns are \code{thermo.plot.new}, \code{thermo.postscript}, \code{label.plot} and \code{water.lines}.

}


\seealso{
  For some of the functions on which these utilities depend or were modeled, see  \code{\link{paste}}, \code{\link{substr}}, \code{\link{tolower}}, \code{\link{par}} and \code{\link{text}}.

}

\examples{
  \dontshow{data(thermo)}
  ## string to character
  s2c("hello world")
  s2c("hello world",sep=" ",keep.sep=FALSE)
  s2c("3.141592",sep=c(".","9"))
  s2c("3.141592",sep=c(".","9"),move.sep=TRUE)
  # character to string
  c2s(aminoacids())
  c2s(aminoacids(),sep=".")

  ## Spearman's rho
  spearman(c(1,2,3),c(2,3,4))  # 1
  spearman(c(1,2,3),c(4,3,2))  # -1

  ## argument processing
  eos.args("hkf",c("g","H","S","cP","V","kT","e"))
  ## produces an error because "Q" is not allowed in water.SUPCRT92
  \dontrun{eos.args("hkf",c("G","H","S","Cp","V","kT","E","Q"))}
  thermo$opt$water <- "IAPWS"  # needed for p and n in next line
  eos.args("water",c("p","u","cv","psat","rho","n","q","x","y","epsilon"))
  TP.args(c(273.15,373.15))
  TP.args(c(273.15,373.15),"Psat")
  TP.args(c(273.15,373.15),c(100,100,200,200))
  state.args(c("AQ","GAS"))
  state.args(c("a","l","liq"))

  ## converting among Gibbs, enthalpy, entropy
  GHS("H") # entropy of H (element)
  # calculate enthalpy of formation of arsenopyrite 
  GHS("FeAsS",DG=-33843,S=68.5) 
  # return the value of DG calculated from DH and S
  # cf. -56687.71 from subcrt("water")
  GHS("H2O",DH=-68316.76,S=16.7123)  
 
  ## count selected elements in a formula
  t <- makeup("H2O")
  expand.formula(c("H","O"),t)
  expand.formula(c("C","H","S"),t)

  ## count amino acids in a sequence
  aminoacids("GGSGG")
  aminoacids("WhatAmIMadeOf?")

  ## count nucleobases in a sequence
  nucleicacids("ACCGGGTTT")
  # the DNA complement of that sequence
  nucleicacids("ACCGGGTTT",comp="DNA")
  # the RNA complement of the DNA complement
  n <- nucleicacids("ACCGGGTTT",comp="DNA",comp2="RNA")
  # the formula of the RNA complement
  nucleicacids(n,type="RNA")

  ## calculate protein length
  protein.length("LYSC_CHICK")
  # another way to do it
  basis("CHNOS")
  species("LYSC_CHICK")
  protein.length(species()$ispecies)
  # another way to do it
  ip <- protein("LYSC","CHICK")
  protein.length(-ip)

  ## heat capacity as a function of temperature
  ## (Makhatadze & Privalov, 1990) units: J mol-1
  MP90.cp(c(5,25,50,75,100,125),"LYSC_CHICK")

  ## properties of phase transitions
  t <- info("enstatite")
  # (dP/dT) of transitions
  dPdTtr(t)  # first transition
  dPdTtr(t+1) # second transition
  # temperature of transitions (Ttr) as a function of P
  Ttr(t,P=c(1,10,100,1000))
  Ttr(t,P=c(1,10,100,1000))
  
  ## nominal carbon oxidation states
  ZC("CHNOSZ")
  t <- info(info("LYSC_CHICK"))
  ZC(t$formula)

  ## the basis stoichiometry of a made-up species
  # warns because P isn't in our basis
  basis("CHNOS")
  basis.comp("SPONCH")

  ## describing the basis species
  basis("CHNOSe")
  describe(thermo$basis)
  describe(thermo$basis,T=NULL,P=NULL)

  ## mass and entropy of compounds of elements
  element("CH4")
  element(c("CH4","H2O"),"mass")
  element("Z")   # charge
  # same mass, opposite entropy as charge
  element("Z-1") # i.e., electron

  ## scale logarithms of activity
  ## suppose we have two proteins whose lengths
  ## are 100 and 200; what are the logarithms
  ## of activity of the proteins that are equal to
  ## each other and that give a total activity of
  ## residues equal to unity?
  logact <- c(-3,-3)  # could be any two equal numbers
  length <- c(100,200)
  logact.tot <- 0
  loga <- unitize(logact,length,logact.tot)
  # the proteins have equal activity
  stopifnot(identical(loga[1],loga[2]))
  # the sum of activity of the residues is unity
  stopifnot(isTRUE(all.equal(sum(10^loga * length),1)))
  ## now, what if the activity of protein 2 is ten
  ## times that of protein 1?
  logact <- c(-3,-2)
  loga <- unitize(logact,length,logact.tot)
  # the proteins have unequal activity
  stopifnot(isTRUE(all.equal(loga[2]-loga[1],1)))
  # but the activities of residues still add up to one
  stopifnot(isTRUE(all.equal(sum(10^loga * length),1)))

  ## modify/add species
  info(t <- info("alanine","cr"))
  mod.obigt("alanine",state="cr",G=0,H=0,S=0)
  # now the values of G, H, and S are inconsistent
  # with the elemental composition of alanine
  info(t)
  # add a species
  mod.obigt("myname",formula="CHNOSZ",G=0,H=0)
  info(t <- info("myname"))
}

\references{

  Alberty, R. A., 2003. \emph{Thermodynamics of Biochemical Reactions}, John Wiley & Sons, Hoboken, New Jersey, 397 p. \url{http://www.worldcat.org/oclc/51242181}

  Buvet, R., 1983. General criteria for the fulfillment of redox reactions, in \emph{Bioelectrochemistry I: Biological Redox Reactions}, Milazzo, G. and Blank, M., eds., Plenum Press, New York, p. 15-50. \url{http://www.worldcat.org/oclc/9282370}
 
  Hendrickson, J. B., Cram, D. J., and Hammond, G. S., 1970. \emph{Organic Chemistry}, 3rd ed., McGraw-Hill, New York, 1279 p. \url{http://www.worldcat.org/oclc/78308}

  Makhatadze, G. I. and Privalov, P. L., 1990. Heat capacity of proteins. 1. Partial molar heat capacity of individual amino acid residues in aqueous solution: Hydration effect \emph{J. Mol. Biol.}, 213, 375-384. \url{http://dx.doi.org/10.1016/S0022-2836(05)80197-4}

 
}

\keyword{misc}
