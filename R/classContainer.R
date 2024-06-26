#####Annotate FIA
#An object to contain all the information from an FIA acquisition.


#' An S4 class to represent an FIA experiments.
#'
#' The S4 class also includes all the informations about processing, and the detected signals
#' are stored.
#'
#' @include noiseEstimator.R
#' @param object A proFIAset object.
#' @export
#' @slot peaks A matrix containg all the peaks which have been detected in each
#' individual files.
#' \itemize{
#'     \item mzmin the minimum value of the mass traces in the m/z dimension.
#'     \item mzmax the maximum value of the mass traces in the m/z dimension.
#'     \item scanMin the first scan on which the signal is detected.
#'     \item scanMax the last scan on which the signal is detected.
#'     \item areaIntensity the integrated area of the signal.
#'     \item maxIntensity the maximum intensity of the signal.
#'     \item solventIntensity the intensity of the solvent, 0 means that no significant
#'     solvent was detected.
#'     \item corSampPeak An idicator of matrix effect, if it's close to 1, the compound
#'     does not suffer from heavy matrix effect, if it is inferior to 0.5, the compound
#'     suffer from heavy matrix effect.
#'     \item signalOverSolventRatio The ratio of the signal max intensity on the oslvent max intensity.
#' }
#' @include noiseEstimator.R
#' @slot group  A matrix containing the information on the groups done between all the
#' acquisitions.
#' \itemize{
#'     \item mzMed the median value of group in the m/z dimension.
#'     \item mzMin the minimum value of the group in the m/z dimension.
#'     \item mzMax the maximum value of the group in the m/z dimension.
#'     \item scanMin the first scan on which the signal is detected.
#'     \item scanMax the last scan on which the signal is detected.
#'     \item nPeaks The number of peaks grouped in a group.
#'     \item meanSolvent The mean of solvent in the acquisition.
#'     \item pvalueMedian The median p-value of the group.
#'     \item corMean The mean of the matrix effect indicator.
#'     \item signalOverSolventMean The mean of ratio of the signal max
#'     intensity on the oslvent max intensity.
#'     \item corSd The standard deviation of the matrix effect indicator.
#'     \item timeShifted Is the peak shifted compared to the injection peak.
#' }
#' @slot groupidx The row of the peaks corresponding to each group
#' in \code{peaks}.
#' @slot step The step of processing of the experiment.
#' @slot path The path of the experiment.
#' @slot classes A table with two columns, "rname" the absolute path
#' of a file, and "group" the class to which the file belong.
#' @slot dataMatrix A matrix variables x experiments suitable for
#' statistical analysis.
#' @slot noiseEstimation A model of noise as estimated by
#' \link{estimateNoiseMS}
#' @slot injectionPeaks A list of the injection peak which have been
#' detected for each experiment.
#' @slot injectionScan A numeric vector giving the scan of injection of
#' sample.
#' @aliases proFIAset-class dataMatrix peaks groupMatrix phenoClasses injectionPeaks
setClass(
    "proFIAset",
    slot = list(
        peaks = "matrix",
        group = "matrix",
        groupidx = "list",
        step = "character",
        path = "character",
        classes = "data.frame",
        dataMatrix = "matrix",
        noiseEstimation = "noiseEstimation",
        injectionPeaks = "list",
        injectionScan = "numeric"


    ),
    prototype = list(
        peaks = matrix(nrow = 0, ncol = 0),
        group = matrix(nrow = 0, ncol = 0),
        groupidx = list(),
        step = "empty",
        path = character(),
        classes = data.frame(),
        dataMatrix = matrix(nrow = 0, ncol = 0),
        noiseEstimation = new("noiseEstimation"),
        injectionPeaks = list(),
        injectionScan = numeric(0)
    )
)
