\name{makeup}
\alias{makeup}
\title{Chemical Formulas and Compositions}
\description{
  Convert between different representations of chemical compositions of species.
}

\usage{
  makeup(compound = "", component = NULL)
}

\arguments{
  \item{compound}{character, the formula of a compound; other types are possible (see below).}
  \item{component}{numeric, coefficient used to sum makeups; other types are possible (see below).}
}

\details{

  This function performs operations on objects representing chemical compositions, or formulas, of compounds. In the simple case, when \code{compound} is a character object containing a chemical formula as a string, the return value is a dataframe object with a single column named \code{count} that lists the coefficients of each of the elements in the formula. With other combinations of arguments, this function can be used to sum compositions, count only specified elements, generate string representations of formulas from dataframes, or perform other operations listed in the table below.

  The elements in a formula are symbolized by an uppercase letter followed by zero or more lowercase letters. Multiple pairs of parentheses are allowed, but may not be nested. Coefficients after the elemental symbols or parentheses can be non-integer values and/or negative. Charge is denoted (at the end of a formula, or preceding an asterisk or colon) by a value which can be non-integer, optionally preceded by \samp{+} or \samp{-} (Hence, \samp{Fe+3} stands for ferric iron; \samp{Fe3+} is interpreted as three iron atoms with a net charge of plus one). An asterisk or a colon indicates to add the composition that follows, which may have a leading numerical coefficient. Some examples of valid formulas are \samp{Rh(SO4)3-4}, \samp{KAl3(OH)6(SO4)2} and \samp{Ca(Al2Si7)O18*6H2O}.

  Charge is treated here like a chemical element. The elemental symbol used for it, \samp{Z}, can appear in chemical formulas, and has atomic number zero and oxidation state +1. For example, the formula of the ammonium ion is written commonly as \samp{NH4+}; in \pkg{CHNOSZ} this same composition can also be identified with \samp{NH4Z}. Likewise, the electron is represented as \samp{Z-1}, and the common symbol for the electron, \samp{e-}, does not appear in chemical formulas in this package. Using this framework, any computations that are based on the complete stoichiometries of charged species, if they impose conservation of elements, account simultaneously for conservation of mass and charge. 

}

\value{
  When this function is called with a single argument that is a string representation of a formula, the elemental composition, or \emph{makeup} is returned as a dataframe of one column with \eqn{n} rows, where \eqn{n} is the number of elements present. If the single argument is other than character, or there is more than one argument, the action of the function is determined by the characteristics of the arguments. In the summary below, \emph{makeup} indicates a dataframe, and \emph{formula} refers to a string (character) representation of a formula. \emph{elements} represents a dataframe with a single row.

  (\emph{formula}) Return \emph{makeup}, the elemental composition of the formula.

  (\emph{makeup},\emph{makeup}) Return \emph{makeup}, the sum of the two arguments.

  (numeric,numeric) Return \emph{makeup}, the sum of the makeups of the species identified in the first argument (corresponding to rownumber of \code{\link{thermo}$obigt}), each multiplied by the coefficient given in the second argument (taken as \eqn{1} if missing).

  (\emph{formula},numeric) Return a list which contains the word (name of or coefficient on an element) beginning at the position in the formula string specified by the second argument, and the length of the word.

  (\emph{formula},\emph{formula}) Return \emph{makeup}, sum of the two arguments.

  (\emph{makeup},character) Return \emph{elements}: transpose \code{makeup}, using all the elements (if the second argument is ""), or only those specified in the second argument.

  (\emph{makeup},TRUE) Return \emph{elements}, \bold{not} the elements, but the stoichiometry of the basis species in the formation reaction for this compound. Requires prior definition of basis species; see \code{\link{basis}}.

  (\emph{elements},FALSE) Return \emph{elements}, change entries in the dataframe that are less than the value in \code{thermo$opt$cutoff} to zero.

  (\emph{formula},character) Return character, the unique elements contained in the \emph{formulas} (character argument is necessary but ignored).

  (list of \emph{makeup},character) Return character, the unique elements contained in the \emph{makeups} (character argument is necessary but ignored).

  (\emph{elements},character) Return character, a formula representation (i.e., "CO2") of the \emph{elements}, (character argument must be "" or the names of elements to consider).
  
}

\seealso{
  Functions used by \code{makeup} to parse character objects include \code{\link{substr}}, and \code{\link{c2s}}. Functions that use \code{makeup} include \code{\link{subcrt}} (to check reaction balance), \code{\link{species}} (to generate coefficients in formation reactions), and \code{\link{GHS}} (to calculate entropies of elements in formulas). The list of elements is contained in \code{\link{thermo}$element}.
}

\examples{
  \dontshow{data(thermo)}
  ## makeups (nrow = number of elements)
  # the composition of a single compound (ncol=1)
  makeup("CO2")
  # negative coefficients
  makeup("C-4O2") # minus four C's
  makeup("C-4O-2")  # interpreted as having charge -2
  makeup("C-4O-2+0") # no charge
  # sum the compositions			
  makeup(c("CO2","CH4"))
  # sum with coefficients	
  makeup(info(c("CO2","CH4")),c(-1,1))
  # this one adds up to zero
  makeup(c("C6H12O6","C2H5OH","CO2"),c(-1,2,2))

  ## formulas (ncol = number of elements)
  # as a dataframe	
  makeup("Zn(CH3CH2CH2CH2CO2)2","")
  print( t <- makeup(c("HCl","H2O"),"") )
  # as a character formula
  makeup(makeup(info("glycinium"),""),"")
  makeup(t,"")
  # C6H12O6
  makeup(makeup(makeup("CHO6H11C5"),""),"")

  ## charged species, electron
  # these return the same:
  makeup("NH4+")
  makeup("NH4Z")
  # constitution of the electron
  makeup(info("e-"))
  # this produce an INcorrect makeup for the electron, 
  # and gives a warning because "e" is not in the
  # table of elements.
  makeup("e-")
  # compositions of made-up compounds
  makeup("CHNOS")
  makeup("CHNOSZ")
  makeup("CHNOSe")  # Se is an element
}

\keyword{misc}
