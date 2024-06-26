#' Process FIA-HRMS datasets.
#'
#' Process FIA-HRMS datasets passing from raw data (mzMl, CDF,mzXML
#' format) to a peak table suitable for statistical anlysis.
#'
#' The full workflow is composed of the following chain of function.
#' \code{proFIAset}=>\code{group.FIA}=>\code{impute.FIA} and the full process is easy to automate using the
#' analyzeAcquisitionFIA which do all the steps easely. Resulting table may be easely
#' exported using the 3 exports function (\code{exportDataMatrix},\code{exportVariableMetadata}\code{exportSampleMetadata})
#' Groups may be visualised using \code{plotEICs} to plot all the EICs of a group. An overview of an acquisition may be otained using
#' the \code{plot} function.
#'
#' @docType package
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @importFrom pracma pinv
#' @importFrom grDevices colorRampPalette rainbow as.raster col2rgb rgb
#' @importFrom graphics abline legend lines matplot plot title barplot par
#'     points polygon layout rasterImage text plot.new segments mtext
#' @importFrom methods new
#' @importFrom stats aggregate coef convolve cor density
#'          dnorm lm loess median nextn pnorm sd
#'          predict quantile t.test integrate na.omit dexp
#' @importFrom utils combn object.size write.table
#' @importFrom xcms xcmsRaw plotTIC rawEIC plotRaw peaks
#' @importFrom Biobase AnnotatedDataFrame ExpressionSet
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom ropls opls
#' @importFrom missForest missForest
#' @importMethodsFrom ropls plot
#' @name proFIA-package
#' @aliases proFIA proFIA-package
NULL
