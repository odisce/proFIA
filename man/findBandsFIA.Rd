% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cWrapper.R
\name{findBandsFIA}
\alias{findBandsFIA}
\title{Detect band in a FIA acquisition}
\usage{
findBandsFIA(
  xraw,
  firstScan = 1,
  lastScan = length(xraw@scantime),
  ppm = 2,
  sizeMin = NULL,
  dmz = 5e-04,
  beginning = NULL,
  end = NULL,
  nIso = 3,
  fracMin = 0.6
)
}
\arguments{
\item{xraw}{An xcmsRaw object as returned by \code{\link[xcms]{xcmsRaw}}.}

\item{firstScan}{The first scan to be considered, 1 for general use.}

\item{lastScan}{The last scan to be considered.}

\item{ppm}{The mass deviation in ppm for point in consecutives scans.}

\item{sizeMin}{The minimum size of a band.}

\item{dmz}{The minimum mass tolerance,useful for small masses}

\item{beginning}{The scan of the injection. May be determined using 
\code{\link{determiningInjectionZone}}.}

\item{end}{The end of injection, may be estimated using \code{\link{determiningInjectionZone}}.}

\item{nIso}{The minimum number of consecutive points for a signal to be detected.}

\item{fracMin}{The minimum fraction of points necessary in beginning:end for a signal to be detected.
contaminated by solvent.}
}
\value{
A vector contaning the inject peak
}
\description{
Detect bands of points with similar mass in conscutive scans.
Points may be moved if a better candidates is found.
}
\examples{
#Getting the path of a file.
if(require(plasFIA)){
  path_raw<-list.files(system.file(package="plasFIA","mzML"),full.names=TRUE)[2]

  #Opening the file with xcms
  xraw<-xcmsRaw(path_raw)

  #Getting the injection scan
  gp<-determiningInjectionZone(xraw)

  #performing band detection.
  tbands<-findBandsFIA(xraw,ppm = 2,sizeMin = gp[3]-gp[1],beginning=gp[1],end=gp[2])
  head(tbands)
}
}
