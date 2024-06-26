% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methodsContainer.R
\name{exportVariableMetadata,proFIAset-method}
\alias{exportVariableMetadata,proFIAset-method}
\alias{exportVariableMetadata}
\title{Export variable metadata.}
\usage{
\S4method{exportVariableMetadata}{proFIAset}(object, filename = NULL)
}
\arguments{
\item{object}{A proFIAset object.}

\item{filename}{If not NULL the result will be written
in filename}
}
\value{
A dataframe with the following columns :
\itemize{
    \item variableID an ID similar to the one of the peak table.
    \item mzMed the median value of group in the m/z dimension.
    \item mzMin the minimum value of the group in the m/z dimension.
    \item mzMax the maximum value of the group in the m/z dimension.
    \item scanMin the first scan on which the signal is detected.
    \item scanMax the last scan on which the signal is detected.
    \item nPeaks The number of peaks grouped in a group.
    \item meanSolvent The mean of solvent in the acquisition.
    \item signalOverSolventPvalue The mean p-value of the group.
    \item corMean The mean of the matrix effect indicator.
    \item SigSolMean The mean of ratio of the signal max
    intensity on the solvent max intensity.
}
}
\description{
Export the variable metadata of an experiment, to
be used for statistical analysis.
}
\examples{
if(require(plasFIA)){
  data(plasSet)
  vtab<-exportVariableMetadata(plasSet)
  head(vtab)
}
}
