% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_eps.R
\name{fun_eps}
\alias{fun_eps}
\title{Returns a function for the echo ability 'eps*' equal to the optimal acoustic cross sectional area of a target divided by the length of the target to the desired power (usually squared). The function fun_eps() has 2 methods depending on the type of the input object 'data'.}
\usage{
fun_eps(data = NULL, Nl = 1, method = c("closest", "linear", "constant"))
}
\arguments{
\item{data}{may be of 2 different types: (1) a function of 1 argument representing the acoustic cross sectional area of the target as a function of frequency, and (2) a list representing the empirical acoustic cross sectional area of a source in intervals or at points of the arguments. Names for the elements of the list adopted from read.TSD(): [['grff']] - A vector of arbitrary length representing the grid vector of acoustic frequency, and [['eps*']] - A vector no shorter than length(grff) holding the echo ability 'eps*' corresponding to the frequency grid vector 'grff'.}

\item{Nl}{is the number of targets from which the echo is simulated.}

\item{method}{defines the interpolation method.}
}
\value{

}
\description{
Returns a function for the echo ability 'eps*' equal to the optimal acoustic cross sectional area of a target divided by the length of the target to the desired power (usually squared). The function fun_eps() has 2 methods depending on the type of the input object 'data'.
}
\examples{
\dontrun{}

}
