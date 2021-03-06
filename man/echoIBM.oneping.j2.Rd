% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/echoIBM.oneping.j2_old.R,
%   R/echoIBM.oneping.j2_temp.R, R/echoIBM_setB_T.R
\name{echoIBM.oneping.j2}
\alias{echoIBM.oneping.j2}
\alias{echoIBM_setB_T}
\title{Simulates one radial sphere (sample interval) 'j' of one echo sounder observation based on positions, orientations, sizes and other specifics of each fish in a known (simulated) school. Two for loops used (for each frequency and each beam).}
\usage{
echoIBM.oneping.j2(j, data)

echoIBM.oneping.j2(j, data)

echoIBM_setB_T(data, direction = 1, beamnr = 1)
}
\arguments{
\item{j}{is the sampling interval index.}

\item{data}{is the list of data for the simulation, containing all of the following elements:
\itemize{
  \item fishdirT
  \item fishposV
  \item dira
  \item dire
  \item wavenumber
  \item lthesel
  \item sllf
  \item indi
  \item bpt*
  \item rad*
  \item pbp*
}}

\item{direction}{Integer: 1 is emision and 2 is reception.}
}
\value{





}
\description{
Simulates one radial sphere (sample interval) 'j' of one echo sounder observation based on positions, orientations, sizes and other specifics of each fish in a known (simulated) school. Two for loops used (for each frequency and each beam).

Simulates one radial sphere (sample interval) 'j' of one echo sounder observation based on positions, orientations, sizes and other specifics of each fish in a known (simulated) school. Two for loops used (for each frequency and each beam).

Gets the beam pattern values from the transducer, either at emision or reception, depending on the value of direction.
}
\examples{
\dontrun{}

\dontrun{}

\dontrun{}

}
