% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methodsContainer.R
\name{plotFlowgrams,proFIAset-method}
\alias{plotFlowgrams,proFIAset-method}
\alias{plotFlowgrams}
\title{Plot raw temporal profiles of the selected group.}
\usage{
\S4method{plotFlowgrams}{proFIAset}(
  object,
  index = NULL,
  mz = NULL,
  subsample = NULL,
  ppm = 5,
  margin = 2e-04,
  posleg = c("topright", "bottomright", "bottom", "bottomleft", "left", "topleft",
    "top", "right", "center"),
  title = NULL,
  scaled = FALSE,
  area = FALSE,
  ...
)
}
\arguments{
\item{object}{A proFIAset object.}

\item{index}{The index of the group to be plotted.}

\item{mz}{An mz value to be looked for only used if index is null.
the research use the \code{\link{findMzGroup}} function.}

\item{subsample}{A subset of sample to be plotted.}

\item{ppm}{The tolerance for the research if mz is provided.}

\item{margin}{An area outer the EICs mz range on which the EIC may be extended.}

\item{posleg}{The position of the legend on the figure. See \code{\link[graphics]{legend}}.}

\item{title}{An optional vector of title for the plot. Need to be of the same}

\item{scaled}{Shall all the EIC be put on the same scale with maximum to 1.}

\item{area}{Shall the detectged area be plotted using transparency.}

\item{...}{Supplementary graphical parameters to be passed to lines.
length than index.}
}
\value{
No returned value
}
\description{
Plot raw temporal profiles from \code{\link{proFIAset}} object
corresponding to one or more molecules. The function will priorize
index, only using mz if index is set to \code{NULL}. As this require to come
back to the raw data, this can take some times, an we don't recommend using it
on more than 30 files. It is better to choose 30 files in the proFIAset object.
}
\examples{
if(require(plasFIA)){
    data(plasMols)
    data(plasSet)
    plotFlowgrams(plasSet,mz=plasMols[7,"mass_M+H"])
}
}
