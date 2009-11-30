\name{eqdata}
\alias{eqdata}
\title{Read data from an EQ6 output file}
\description{
  This function extracts computational results for aqueous species at each step of reaction progress in an EQ6 output file. The results are written to a comma-separated value file that can be read by other programs.
}

\usage{
  eqdata(file, species, prop = "log act")
}

\arguments{
  \item{file}{character, path to EQ6 output file.}
  \item{species}{character, name(s) of aqueous species.}
  \item{prop}{character, property to get.}
}

\details{
  This function can be used to extract values from EQ6 (Wolery et al., 1990) output files and save them in a format readable by spreadsheets. It is independent of other functions in \pkg{CHNOSZ}. For more information, see the document in \url{../doc/eqdata.pdf}.
}

\value{
  Creates a comma-separated value file in the working directory, and invisibly returns the same data.
}

\references{
  Wolery, T. J., Jackson, K. J., Bourcier, W. L., Bruton, C. J., Viani, B. E., Knauss, K. G. and Delany, J. M., 1990. Current status of the EQ3/6 software package for geochemical modeling. \emph{ACS Symp. Ser.}, 416, 104-116. \url{http://dx.doi.org/10.1021/bk-1990-0416.ch008}
}

