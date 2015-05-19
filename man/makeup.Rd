\name{makeup}
\alias{makeup}
\alias{count.charge}
\alias{count.formulas}
\alias{count.elements}
\title{Parse Chemical Formulas}
\description{
  Count the charge and number of elements in a chemical formula.
}

\usage{
  makeup(formula, multiplier = 1, sum = FALSE, count.zero = FALSE)
  count.charge(formula)
  count.formulas(formula)
  count.elements(formula)
}

\arguments{
  \item{formula}{character, a chemical formula}
  \item{multiplier}{numeric, multiplier for the elemental counts in each formula}
  \item{sum}{logical, add together the elemental counts in all formulas?}
  \item{count.zero}{logical, include zero counts for elements?}
}

\details{
  \code{makeup} parses a chemical formula expressed in string notation, returning the numbers of each element in the formula. The formula may carry a charge, indicated by a + or - sign, possibly followed by a magnitude, after the uncharged part of the formula. The formula may have multiple subformulas enclosed in parentheses (but the parentheses may not be nested), each one optionally followed by a numeric coefficient. The formula may have one suffixed subformula, separated by \samp{*} or \samp{:}, optionally preceded by a numeric coefficient. All numbers may contain a decimal point.

  \code{makeup} calls a sequence of supporting functions depending on specific characters present in the formula. If the formula has a charge, it is first parsed using \code{count.charge}. If the formula has subformulas, in parentheses or as a suffix, they are separated and counted using \code{count.formulas}. Finally, the elements in each subformula are counted using \code{count.elements}.

 \code{count.elements} processes a simple chemical formula that must adhere to the following pattern: it starts with an elemental symbol; all elemental symbols start with an uppercase letter, and are followed by another elemental symbol, a number (possibly fractional, possibly signed), or nothing (the end of the formula).

  Any sequence of one uppercase letter followed by zero or more lowercase letters is recognized as an elemental symbol by \code{count.elements}, but \code{makeup} will issue a warning for elemental symbols that are not present in \code{\link{thermo}$element}.

  \code{makeup} can handle numeric and length > 1 values for the \code{formula} argument. If the argument is numeric, it identifies row number(s) in \code{thermo$obigt} from which to take the formulas of species. If \code{formula} has length > 1, the function returns a list containing the elemental counts in each of the formulas. If \code{count.zero} is TRUE, the elemental counts for each formula include zeros to indicate elements that are only present in any of the other formulas.

  The \code{multiplier} argument must have either length = 1 or length equal to the number of formulas. The elemental count in each formula is multiplied by the respective value. If \code{sum} is true, the elemental counts in all formulas (after any multiplying) are summed together to yield a single bulk formula. 
}

\value{
  \code{count.charge} returns a list with named elements \code{charge} and \code{uncharged}, indicating, respectively, the numeric value of the charge, and the original formula string excluding the charge. \code{count.formulas} returns a numeric vector with names refering to each of the subformulas or the whole formula if there are no subformulas. \code{count.elements} and \code{makeup} return numeric vectors with names refering to each of the elemental symbols in the formula. For \code{makeup}, if more than one formula is provided, a list of numeric vectors is returned, unless \code{sum} is TRUE.
}

\examples{
\dontshow{data(thermo)}
# the composition of a simple compound
makeup("CO2")     # 1 carbon, 2 oxygen
# the formula of lawsonite, with a parenthetical part and a suffix
makeup("CaAl2Si2O7(OH)2*H2O")
# fractional coefficients are ok
redfield <- c(106, 16, 1)
reddiv10 <- makeup("C10.6N1.6P0.1")
stopifnot(10*reddiv10 == redfield)

# the coefficient for charge is a number with a *preceding* sign
# e.g., ferric iron, with a charge of +3 is expressed as
makeup("Fe+3")
# transcribing the formula the way it appears in many 
# publications produces a likely unintended result:
# 3 iron atoms and a charge of +1
makeup("Fe3+")

# these all represent a single negative charge, i.e., electron
makeup("-1")
makeup("Z0-1")
makeup("Z-1+0")

# hypothetical compounds with negative numbers of elements
makeup("C-4(O-2)")   # -4 carbon, -2 oxygen
makeup("C-4O-2")     # -4 carbon,  1 oxygen, -2 charge
makeup("C-4O-2-2")   # -4 carbon, -2 oxygen, -2 charge

# the 'sum' argument can be used to check mass and charge
# balance in a chemical reaction
formula <- c("H2O", "H+", "Z0-1", "O2")
(mf <- makeup(formula, c(-1, 2, 2, 0.5), sum=TRUE))
stopifnot(all(mf==0))
}

\seealso{
  Many other functions in CHNOSZ rely on \code{makeup} for their operation: \code{\link{mass}} and \code{\link{entropy}} for calculating properties of chemical compounds from their elements; \code{\link{basis}} and \code{\link{i2A}} for constructing stoichiometric matrices (with \code{count.zero=TRUE}); \code{\link{subcrt}} for checking mass balance of chemical reactions; and others.
}

\concept{Secondary thermodynamic modeling}
