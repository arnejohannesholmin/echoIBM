% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/echoIBM.calibrate.R, R/echoIBM.calibrate2.R
\name{echoIBM.calibrate}
\alias{echoIBM.calibrate}
\alias{echoIBM.calibrate2}
\title{Calibrates (and writes the calibration data to file) the simulation of observations from the MS70 sonar performed by echoIBM(). Not used in any of the other echoIBM code files.}
\usage{
echoIBM.calibrate(
  directory = NULL,
  cruise = "S2009116",
  event = 1,
  esnm = "MS70",
  type = c("sv", "ts"),
  dire = NULL,
  runs = 1,
  add = 30,
  N = 1e+05,
  j = 100,
  sgbs = 1e-06,
  max.memory = 1e+09,
  method = c("closest", "linear"),
  max.radius = 0.2,
  seed = 0,
  cores = 1,
  usemean = FALSE,
  nameadd = ""
)

echoIBM.calibrate2(
  directory = NULL,
  event = 1,
  type = c("sv", "ts"),
  dire = NULL,
  runs = 1,
  add = 30,
  N = 1e+05,
  j = 100,
  sgbs = 1e-06,
  max.memory = 1e+09,
  method = c("closest", "linear"),
  max.radius = 0.2,
  seed = 0,
  cores = 1,
  usemean = FALSE,
  nameadd = ""
)
}
\arguments{
\item{directory}{the path to the directory in which to put the calibration file (use TRUE for the default directory file.path(echoIBM_datasets_directory(), "Resources", "Calibration")).}

\item{cruise}{is either the idenfication number of the cruise, given as specified by the IMR (yyyynnn), or the path to the directory containing the event on which the calibration should be based.}

\item{event}{is the identifier of the event, either given as the number of the event, a string contained in the name of the event, or the path of the event directory.}

\item{esnm}{is the name of the acoustical instrument, given as a four character string. See sonR_implemented() for the implemented systems. May be given in 'data', and in lower case.}

\item{type}{is "sv" for calibrating the simulation model with respect to volume scattering 'sv' and "ts" for calibrating the simulation model with respect to backscattering cross section / target strength.}

\item{runs}{is an integer giving the number of simulated pings of which the mean is taken as the final calibration factors.}

\item{add}{is the angle in degrees to add on either end of the insonifying volume of the acoustical instrument, defining the volume in which point targets are distributed.}

\item{N}{is the number of point targets used in the calibration.}

\item{j}{is the sampling interval to use for the calibration (should be larger than 2).}

\item{sgbs}{is the backscattering cross section of the point targets.}

\item{max.memory}{is used in echoIBM to set the maximum memory for echoIBM to occupy.}

\item{filetag}{is a string to add to the file name.}
}
\description{
Calibrates (and writes the calibration data to file) the simulation of observations from the MS70 sonar performed by echoIBM(). Not used in any of the other echoIBM code files.

Calibrates (and writes the calibration data to file) the simulation of observations from the MS70 sonar performed by echoIBM(). Not used in any of the other echoIBM code files.
}
\examples{
\dontrun{}

\dontrun{}

}
