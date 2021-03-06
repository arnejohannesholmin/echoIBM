% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/warnings2character.R
\name{warnings2character}
\alias{warnings2character}
\title{Returns the warning messages stored in 'x' which is an object of type "warnings" as retunred from warnings().}
\usage{
warnings2character(
  x = warnings(),
  numbering = TRUE,
  header = TRUE,
  width.cutoff = 50L,
  nlines = 2L
)
}
\arguments{
\item{x}{is the warnings-object to translate to character string vector.}

\item{numbering}{is TRUE if numbering should be adde to the warnings, as is done when warnings() is used.}

\item{header}{is TRUE the first string of the output is "Warning message" or "Warning messages".}

\item{width.cutoff}{is an integer in [20, 500] determining the cutoff at which line-breaking is tried.}

\item{nlines}{is an integer giving the maximum number of lines to produce. Negative values indicate no limit.}
}
\value{

}
\description{
Returns the warning messages stored in 'x' which is an object of type "warnings" as retunred from warnings().
}
\examples{
\dontrun{}

}
