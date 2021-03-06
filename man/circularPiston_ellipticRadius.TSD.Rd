% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/circularPiston_ellipticRadius.TSD.R
\name{circularPiston_ellipticRadius.TSD}
\alias{circularPiston_ellipticRadius.TSD}
\title{Returns the theoretical far field approximation of the acoustic intensity from a circular baffled source (circular piston) of radius as an elliptical function of azimuth angle 'data$dira' given the radii along the x-axis, data$wnsz[,1], and along the y-axis, data$wnsz[,2], as a function of incidence angle 'data$dire'.}
\usage{
circularPiston_ellipticRadius.TSD(data, appr = NULL)
}
\arguments{
\item{data}{is a list containing the elements:}

\item{appr}{is an approximating function used to replace the expression (2*besselJ(ang,1)/ang)^2.}
}
\value{

}
\description{
Returns the theoretical far field approximation of the acoustic intensity from a circular baffled source (circular piston) of radius as an elliptical function of azimuth angle 'data$dira' given the radii along the x-axis, data$wnsz[,1], and along the y-axis, data$wnsz[,2], as a function of incidence angle 'data$dire'.
}
\examples{
\dontrun{}

}
