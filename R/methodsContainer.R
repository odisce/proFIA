#' Ge the class organization of directory.
#'
#' Find the classes organization of a directory, and return a
#' table. This function is called by proFIAset, and is useful
#' to check the structure which will be considered by a proFIAset object.
#'
#' @export
#' @param files The path to the experiment directory/
#' @return a table containing two columns, the absolute paths
#' of the files and the classes of the acquisition, as given by
#' subdirectories.
#' @aliases acquisitionDirectory
#' @examples
#' if(require(plasFIA)){
#'     path<-system.file(package="plasFIA","mzML")
#'     tabClasses<-acquisitionDirectory(path)
#'     tabClasses
#' }
acquisitionDirectory <- function(files = NULL) {
    if (!file.exists(files))
        stop("path must be a valid directory or file.")
    ###Case where it is a single path.
    if (!dir.exists(files))
        return(NA)

    filepattern <-
        c(
            "[Cc][Dd][Ff]",
            "[Nn][Cc]",
            "([Mm][Zz])?[Xx][Mm][Ll]",
            "[Mm][Zz][Dd][Aa][Tt][Aa]",
            "[Mm][Zz][Mm][Ll]"
        )
    filepattern <-
        paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    if (is.null(files))
        files <- getwd()
    info <- file.info(files)
    listed <- list.files(
        files[info$isdir],
        pattern = filepattern,
        recursive = TRUE,
        full.names = TRUE
    )
    if (length(listed) == 0)
        stop("impossible to find any supported data format.")
    groupf <- dirname(listed)
    allgroup <- unique(groupf)
    lg <- character()
    for (i in 1:length(allgroup)) {
        lg <- c(lg, listed[which(groupf == allgroup[i])])
    }
    message(
        length(allgroup),
        " directories have been found : ",
        paste(basename(allgroup), collapse = " "),
        "\n"
    )
    df <- data.frame(
    	rname <- as.character(lg),
    	group <- as.factor(basename(groupf)),
    	stringsAsFactors =
    		FALSE
    )
    colnames(df) <- c("path","class")
    return(df)
}



#' @export
 setGeneric("phenoClasses", function(object, ...)
     standardGeneric("phenoClasses"))
# #' Extract the classes table from a proFIAset object
# #'
# #' Extract the classes table from a proFIAset object.
# #'
# #' @export
# #' @param object A proFIAset object.
# #' @return a table containing two columns, the absolute paths
# #' of the files and the classes of the acquisition, as given by
# #' subdirectories.
# #' @aliases getPhenoClasses getPhenoClasses,proFIAset-method
# #' @examples
# #' if(require("plasFIA")){
# #'   data(plasSet)
# #'   tabClasses<-getPhenoClasses(plasSet)
# #' }

#' @describeIn proFIAset Extract the classes and the paths of the samples.
#' @export
setMethod("phenoClasses", "proFIAset", function(object) {
    return(object@classes)
})

#' @export
 setGeneric("dataMatrix", function(object, ...)
     standardGeneric("dataMatrix"))
# #' Extract the data matrix from a proFIAset object
# #'
# #' Extract the data matrix from a proFIAset object. If
# #' the matrix have not been created, it is created.
# #'
# #' @export
# #' @param object A proFIAset object.
# #' @param ... Supplementary arguments to be passed to makeDataMatrix
# #' if the matrix have not been created.
# #' @return a matrix with groups as rows and samples as columns.
# #' @aliases dataMatrix dataMatrix,proFIAset-method
# #' @seealso \code{\link{exportPeakTable}} for a more personnalized output.
# #' @examples
# #' if(require("plasFIA")){
# #'   data(plasSet)
# #'   dataMatrix<-dataMatrix(plasSet)
# #' }

#' @describeIn proFIAset Extract the dataMatrix containg variables as rows and samples as columns
#' @export
setMethod("dataMatrix", "proFIAset", function(object) {
    if (nrow(object@dataMatrix) == 0) {
        stop("Use makeDataMatrix to create the data matrix.")
    }
    return(object@dataMatrix)
})


#' @export
 setGeneric("groupMatrix", function(object, ...)
     standardGeneric("groupMatrix"))
# #' Extract the group matrix from a proFIAset object
# #'
# #' Extract the group matrix from a proFIAset object.
# #' We recommend to use the \code{\link{exportPeakTable}}
# #' function for a more formalized output.
# #'
# #' @export
# #' @param object A proFIAset object.
# #' @return A matrix with groups as rows and informations on group as columns.
# #' See the slot group of the proFIAset class \code{\link{proFIAset-class}} for more detail.
# #' @aliases groupMatrix groupMatrix,proFIAset-method
# #' @seealso \code{\link{exportPeakTable}} for a more personnalized output and
# #' \code{\link{proFIAset-class}} for a description of the columns of the group matrix.
# #' @examples
# #' if(require("plasFIA")){
# #'   data(plasSet)
# #'   groupMatrix<-groupMatrix(plasSet)
# #' }

#' @describeIn proFIAset Extract the matrix of group, see \code{\link{exportPeakTable}} for better output.
#' @export
setMethod("groupMatrix", "proFIAset", function(object) {
    return(object@group)
})


# #' Extract the peaks matrix from a proFIAset object
# #'
# #' Extract all the peak detected from an FIA acquisition.
# #' The output may be space consuming.
# #'
# #' @export
# #' @param object A proFIAset object.
# #' @return a matrix with peaks as rows and informations on peaks as columns.
# #' @aliases peaks peaks,proFIAset-method
# #' @seealso \code{\link{proFIAset-class}} to see the fieds on the peaks table.
# #' @examples
# #' if(require("plasFIA")){
# #'   data(plasSet)
# #'   peaksMatrix<-peaks(plasSet)
# #' }

#' @describeIn proFIAset Extract all the signals detected in individual samples.
#' @export
setMethod("peaks", "proFIAset", function(object) {
    return(object@peaks)
})

#' @export
 setGeneric("injectionPeaks", function(object, ...)
     standardGeneric("injectionPeaks"))
# #' Extract the injection peaks from a proFIAset object
# #'
# #' Extract all the regressed peak table from a proFIAset object.
# #'
# #' @export
# #' @param object A proFIAset object.
# #' @return a list with the injection peak correspoing to each sample.
# #' @aliases injectionPeaks injectionPeaks,proFIAset-method
# #' @seealso \code{\link{proFIAset-class}} to see the fieds on the peaks table.
# #' @examples
# #' if(require("plasFIA")){
# #'   data(plasSet)
# #'   lInjPeak<-injectionPeaks(plasSet)
# #' }

#' @describeIn proFIAset Extract all the regressed injection peaks.
#' @export
setMethod("injectionPeaks", "proFIAset", function(object) {
    return(object@injectionPeaks)
})


#' Process FIA experiment.
#'
#' Processes an experiment ordered as a tree of files, and return a
#' proFIAset object. TO check the funirshed structure you may use the
#' \code{\link{acquisitionDirectory}} function.
#' @export
#'
#' @param path The path of the files to be processed.
#' @param ppm The tolerance of the algorithms in deviation between scans.
#' @param parallel Shall parallelism using \pkg{BiocParallel} be used.
#' @param noiseEstimation Shall noise be estimated (recommended)
#' @param BPPARAM A BiocParallelParam object to be used for parallelism
#' if parallel is TRUE.
#' @param graphical Shall the plot of the regressed injection peak be shown.
#' @param ... Supplementary arguments to be passed to findFIAsignal, see
#' \code{\link{findFIASignal}}.
#' @aliases proFIAset
#' @seealso  To obtain more detail about the output see \code{\link{proFIAset-class}}.
#' @return A proFIAset object.
#' @examples
#' if(require("plasFIA")){
#'
#' pathplas<-system.file(package="plasFIA","mzML")
#'
#' #Parameters are defined.
#' ppm<-2
#'
#' plasSet<-proFIAset(pathplas,ppm=ppm,parallel=FALSE)
#' plasSet
#'
#' ###An example using parallelism with a snow cluster using BiocParallel package.
#' \dontrun{plasSet<-proFIAset(pathplas,ppm=ppm,parallel=TRUE,BPPARAM=bpparam("SnowParam"))}
#' }
proFIAset <-
    function(path,
             ppm,
             parallel = TRUE,
             BPPARAM = NULL,
             noiseEstimation = TRUE,
             f = c("regression","TIC"),
             graphical=FALSE,
             ...) {
        pFIA <- new("proFIAset")
        pFIA@classes <- acquisitionDirectory(path)
        f <- match.arg(f)
        #print(object@classes)
        pFIA@path <- path
        message(
            nrow(pFIA@classes),
            " files found in directories ",
            unique(as.character(pFIA@classes[, 2])),
            "\n",
            sep =
                " "
        )
        if(requireNamespace("BiocParallel")&is.null(BPPARAM)){
            BPPARAM <- bpparam()
        }
        ##Quick opening of two random files to evaluate intensity.
        vtab <- c(1,cumsum(table(pFIA@classes[,2])))
        tvmint <- apply(pFIA@classes[vtab,],1,function(x){
            requireNamespace("xcms")
            vxraw <- xcmsRaw(x[1])
            c(quantile(vxraw@env$intensity,probs=0.03),max(vxraw@env$intensity))
        })

        vmint <- 2*max(tvmint[2,])
        vminint <- mean(tvmint[1,])

        ###The nois emodel is estimated.
        nes <- NULL
        if (noiseEstimation) {
            if (parallel & requireNamespace("BiocParallel")) {
                message("Package BiocParallel used.")
                nes <- estimateNoiseMS(
                    pFIA@classes[, 1],
                    ppm = ppm,
                    parallel = parallel,
                    BPPARAM = BPPARAM,
                    maxInt=vmint,
                    minInt=vminint,
                    f=f,
                )
            } else{
                if (parallel) {
                    warning(
                        "BioCParallel package is not installed, impossible to use parallelism. Single core is used."
                    )
                }
                nes <- estimateNoiseMS(pFIA@classes[, 1],
                                             ppm = ppm,
                                             parallel =
                                                 FALSE,
                                       maxInt=vmint,
                                       minInt=vminint)
            }
        }
        # browser()
        nes <- fitModel(nes,absThreshold = NULL)
        # browser()
        if(is.na(nes)){
          nes <- estimateNoiseMS(pFIA@classes[, 1],
                                 ppm = ppm, minInt = 50, maxInt = vmint,
                                 parallel =
                                   FALSE)
          # browser()
          nes <- fitModel(nes,absThreshold = NULL)
          # browser()
        }


        if(graphical){
            plotNoise(nes)
        }
        pFIA@noiseEstimation <- nes
        ###We keep all the bands if possible.


        lRes <- list()
        if (parallel & requireNamespace("BiocParallel")) {
            if (is.null(BPPARAM))
                BPPARAM <- bpparam()
            lRes <- bplapply(
                as.list(pFIA@classes[,1]),
                openAndFindPeaks,
                ppm = ppm,
                es = pFIA@noiseEstimation,
                f=f,
                BPPARAM = BPPARAM,
                ... =
                    ...
            )
        } else{
            if (parallel) {
                warning(
                    "BioCParallel package is not installed, impossible to use parallelism."
                )
            }
            lRes <- sapply(
                pFIA@classes[, 1],
                openAndFindPeaks,
                simplify = FALSE,
                ppm = ppm,
                f = f,
                es = pFIA@noiseEstimation,
                ... =
                    ...
            )
        }
        ###The returned matrix is a res.

        matRes <- lapply(lRes, function(x) {
            x$signals
        })
        lPeaks <- lapply(lRes, function(x) {
            x$injectionPeak
        })
        lBegin <- sapply(lRes, function(x) {
            x$injectionScan
        })

        nPeaks <- lapply(matRes, nrow)
        vecClass <- rep(as.numeric(1:nrow(pFIA@classes)), times = unlist(nPeaks))
        matRes <- do.call(rbind, matRes)
        oldN <- colnames(matRes)
        matRes <- cbind(matRes, vecClass)
        colnames(matRes) = c(oldN, "sample")
        pFIA@peaks <- matRes
        pFIA@step <- "Peaks_extracted"
        pFIA@injectionPeaks <- lPeaks
        pFIA@injectionScan <- lBegin
        message(paste(
            "processus finished",
            nrow(pFIA@classes),
            "samples treated",
            nrow(pFIA@peaks),
            "peaks found."
        ))
        pFIA
    }


#' @include KNN_T.R
setMethod("show", "proFIAset", function(object) {
    if(nrow(object@classes)==0){
        cat("An empty proFIAset object.\n")
        return(invisible(NULL))
    }
    cat("A \"proFIAset\" object containing ",
        length(unique(object@classes[, 2])),
        " classes.\n")
    progress <- switch(object@step,
                      Peak_extracted = "The data have been peak picked.",
                      Samples_grouped = "Grouping has been done between samples.",
                      Matrix_created = "Data matrix has been created.",
                      Fillpeaks = "Missing values have been imputed.")
    cat(progress, "\n")
    if (nrow(object@peaks) > 0) {
        cat(nrow(object@peaks), " peaks detected.\n")
    }

    if (nrow(object@group) > 0) {
        cat(nrow(object@group), " features have been grouped.\n")
    }
    if (nrow(object@dataMatrix) > 0) {
        cat("The data matrix is avalaible.\n")
    }
    memsize <- object.size(object)
    cat("Memory usage:", signif(memsize / 2 ^ 20, 3), "MB\n")
})

switchSummaryPeak<-function(x){
	if(is.na(x)){
		return("Shifted in time")
	}
	if(x==0){
		return("Well-behaved peak")
	}
	if(x==1){
		return("Heavy shape distorsion (corSampPeak<0.2)")
	}

	###Thess case are never uspposed to happens.
	if(x==2){
		return("Shifted in time")
	}
	if(x==3){
		return("Shifted in time and matrix effect")
	}
}


#' Plot a summary of an FIA experiment.
#'
#' Plot a summary of an FIA acquisition. This summary aims to provides an overview of the
#' FIA acquisition and the processinf of FIA acquisition. It includes the following graphs :
#' \itemize{
#'     \item Barplot A barplot giving the number of peak detected in each sample. The number
#'     of shifted peak is indicated, which can indicate wrong FIA-acquisition, as well as the
#'     number of peak badly with a low correlation with the injection peak, which can indicate
#'     a strong matrix effect.
#'     \item Injection_Peaks A representation of all the injection peaks of the acquisition, a
#'     greatly different peak in th einjection may indicate an issue with the injection, or a problem
#'     of regression. If the peaks seems all different, try to use \link{proFIAset} with the f
#'     paramter on TIC, to avoid regression.
#'     \item Density A density plot of the m/z of all the features found. A missing interval at the end
#'     of the mz range may indicate a too low ppm parameter, while a missing interval at the beginning
#'     of the mz range may indicate a too low dmz parameter.
#'     \item PCA A PCA plot to obtain a quick diagnostic of the data. The PCA plot is supposed to be different
#'     before and after imputation.
#' }
#'
#' @export
#' @param x A proFIAset object.
#' @param type Shall the plotting be done by sample or by class for the barplot ?
#' @param ... Not used at the moment.
#'
#' @aliases plot.FIA plot,proFIAset-method plot,proFIAset,ANY,ANY-method
#' @return Nothing
#' @examples
#' if(require("plasFIA")){
#'    data(plasSet)
#'
#'
#'    ####Diagnostic plot after imputation
#'    plot(plasSet)
#'
#'    ####The same plot by classes.
#'    plot(plasSet,type="class")
#'
#'    ####Diagnostic plot before imputation
#'    plasSet <- makeDataMatrix(plasSet)
#'	  plot(plasSet)
#'
#' }
#' @rdname plot
setMethod("plot", "proFIAset",function(x, type=c("sample","class"),...){

	vstep <- 0
	title_pca <- NULL
	vomar <- par("mar")
	par(mar=c(3.1,2.1,2.1,1.1), font = 2, font.axis = 2, font.lab = 2)
	if(x@step=="Peak_extracted"){
		layout(matrix(1:2,ncol=2))
		vstep <- 1
	}
	if(x@step=="Samples_grouped"){
		layout(matrix(c(1,1,2,3),ncol=2,byrow=TRUE))
		vstep <- 2
	}
	if(x@step %in% ("Matrix_created")){
		layout(matrix(c(1,1,1,1,2,4,3,4),ncol=2,byrow=TRUE))
		vstep <- 4
		title_pca <- "PCA (non-imputed data set)"
	}
	if(x@step=="Fillpeaks"){
		layout(matrix(c(1,1,1,1,2,4,3,4),ncol=2,byrow=TRUE))
		vstep <- 4
		title_pca <- "PCA (imputed data set)"
	}

	### Barplot of the number of peaks found.
	type <- match.arg(type)
	###First plotting the number of peak detected by sample.
	#tsample <- table(object@peaks[,"sample"])
	swShiftCor <- 2*x@peaks[,"timeShifted"]+1*(x@peaks[,"corSampPeak"]<0.2)
	vclasses <- sapply(swShiftCor,switchSummaryPeak)
	if(type=="sample"){
		vaggreg <- x@peaks[,"sample"]
	}else{
		vaggreg <- x@classes[x@peaks[,"sample"],2]
	}
	valHeight <- aggregate(vclasses,by=list(vaggreg),FUN=function(x){
		tx <- table(x)
		vres <- as.numeric(tx)
		names(vres) <- names(tx)
		vres
	})
	valHeight <- as.matrix(valHeight[,-1])
	maxy <- apply(valHeight,1,sum)
	if(type=="sample"){
		rownames(valHeight) <- proFIA:::getRawName(x@classes[,1])
	}else{
		rownames(valHeight) <- unique(as.character(x@classes[,2]))
	}
	omar <- par("mar")
	par(mar=c(4.1,2.6,3.1,0.1))
	barVn <- barplot(t(valHeight),
					 ## legend=colnames(valHeight),
					 names.arg = rep("", nrow(valHeight)),
					 col=c("red","blue","green4"),
					 ## ylab="Number of peaks",
					 ## main="Number of peaks",
					 ## ylim=c(0,max(maxy)*1.25),
					 las = 3)
	mtext(substr(rownames(valHeight), 1, 10), side = 1, at = barVn, las = 2, cex = 0.6)
	mtext("Quality of feature flowgrams (in each sample)",
		  line = 1.7,
		  font = 2, col = "black", cex = 0.9)
	mtext("Well behaved",
		  line = 0.2,
		  at = par("usr")[1] + 0.2 * diff(par("usr")[1:2]), col = "green4", cex = 0.8)
	mtext("Shifted",
		  line = 0.2,
		  at = par("usr")[1] + 0.4 * diff(par("usr")[1:2]), col = "blue", cex = 0.8)
	mtext("Significant Matrix Effect",
		  line = 0.2,
		  at = par("usr")[1] + 0.7 * diff(par("usr")[1:2]), col = "red", cex = 0.8)

	par(mar=omar)

	###the injection peaks are plotted.
	plotSamplePeaks(x, diagPlotL = TRUE)

	if(vstep<2) return(invisible(NA))

	omar <- par("mar")
	par(mar=c(3.1,2.6,2.1,1.1))
	plot(density(x@group[,"mzMed"],bw=5),col="red",xlab="",
		 ## ylab="Density of features",
		 main="Density of m/z features",lwd=2)
	mtext("m/z", side = 1, line = 2, cex = 0.7)
	par(mar = omar)

	if(vstep<3) return(invisible(NA))

	tempMatrix <- t(x@dataMatrix)
	tempMatrix[which(is.na(tempMatrix))] <- 0
	res_pca <- suppressMessages(opls(tempMatrix, predI = 2,log10L = TRUE,crossvalI = min(7,nrow(tempMatrix) - 1), fig.pdfC = "none"))
	plot(res_pca,type="x-score",parAsColFcVn=x@classes[,2],parTitleL=FALSE)
	title(title_pca, line = 2.5)
	par(mar=vomar)
	return(invisible(NA))
})




setGeneric("group.FIA", function(object, ...)
    standardGeneric("group.FIA"))
#' Group the peaks of an FIA acquisition.
#'
#' Group the peaks in a FIA experiemnts by clustering under
#' an estimated density based on the accuracy in ppm
#' of the mass spectrometer.
#' @export
#'
#' @param object A proFIAset object.
#' @param ppmGroup A ppm parameter giving the size of the windows
#' considered, we recommend to use 0.5*ppm where ppm is the pp parameter
#' of band detection in the \code{\link{proFIAset}} function.
#' @param sleep If not 0 densities are plotted every \code{sleep} ms.
#' @param dmzGroup A minimum mz deviation value which can be considered over the deviation in ppm
#' if it is higher. This account for low mass deviation.
#' @param fracGroup The minimum fraction of samples of a class required to make
#' a group.
#' @param solvar Shall the group corresponding to solvent signal be kept. This in only useful if solvent
#' signal have been conserved in the \code{\link{findFIASignal}} function.
#' @return A proFIAset object with the \code{group} slot filled. See \code{\link{proFIAset-class}}.
#' @aliases group.FIA group.FIA,proFIAset-method
#' @examples
#' if(require("plasFIA")){
#' #proFIAset object is loaded
#' data(plasSet)
#'
#' #Parameters are defined.
#' ppm_group <- 1
#' dmz_group <- 0.0005
#' frac_group<-0.2
#'
#' plasSet<-group.FIA(plasSet,ppmGroup=ppm_group,fracGroup=frac_group,dmzGroup=dmz_group)
#' plasSet
#' }
setMethod("group.FIA", "proFIAset", function(object,
                                             ppmGroup,
                                             solvar = FALSE,
                                             sleep = 0,
											 dmzGroup=0.0005,
                                             fracGroup = 0.5) {
    ###CHecking that files have been peaks picked
    if (nrow(object@peaks) == 0) {
        stop("object does not contain any signals, impossible to group. See ?findFIAsignal")
    }
    tabP <- peaks(object)
    porder <- order(tabP[, "mz"])
    peaksL <- tabP[porder,]
    col_val <- NULL
    col_vec <- NULL
    uclasses <- NULL
    lty_vec <- NULL
    lwd_vec <- NULL
    nPoints <- 1024
    mzrange <- range(peaksL[, "mz"])
    #1500 is thet maximum mass of a metabolite. 1024 the number of points used in the density
    inter <- mzrange[2] * ppmGroup * nPoints / (20 * 10 ^ 6)
    message(
        "A mass interval of ",
        sprintf("%0.4f", inter),
        "m/z will be used for the density estimation."
    )
    if (sleep > 0) {
        uclasses <- unique(object@classes[, 2])
        col_val <- rainbow(length(uclasses))
        col_vec <- col_val[as.factor(object@classes[peaksL[, "sample"], 2])]
        lty_vec <- rep(1, length(col_val))
        lwd_vec <- rep(1, length(col_val))
    }
    mzrange <- c(floor(mzrange[1] / 50) * 50, ceiling(mzrange[2] / 50) * 50)
    massInter <- seq(mzrange[1], mzrange[2], inter / 2)
    masspos <- findEqualGreaterM(peaksL[, "mz"], massInter)
    num_group <- 0

    tabclasses <- table(as.character(object@classes[, 2]))
    legList <- c(
        "mzMed",
        "mzMin",
        "mzMax",
        "sizeSamp",
        "scanMin",
        "scanMax",
        "nPeaks",
        "meanSolvent",
        "timeShifted"
    )
    boolpval <- FALSE
    NamesPeakList <- colnames(tabP)
    boolpval <- ("signalOverSolventPvalue" %in% NamesPeakList)
    if (boolpval) {
        legList <- c(legList, "signalOverSolventPvalueMean")
    }
    boolcor = ("corSampPeak" %in% NamesPeakList)
    if (boolcor) {
        legList <- c(legList, "corSampPeakMean", "corSampPeakSd")
    }
    boolss = ("signalOverSolventRatio" %in% NamesPeakList)
    if (boolss) {
        legList <- c(legList, "signalOverSolventRatioMean")
    }
    groupmat <- matrix(nrow = 512, ncol = length(legList))
    colnames(groupmat) <- legList
    listgroup <- vector(mode = "list", 512)
    maxPeaks <- max(tabP[, "sample"])
    pos <- 0
    previousMz <- NULL
    imessage <- round(seq(1,(length(massInter) - 2),length=11))
    for (i in seq(length = length(massInter) - 2)) {
        if (i %in% imessage)
            message(floor(i/length(massInter)*10)*10, " :", num_group," ",appendLF=FALSE)
        mem <- i
        start <- masspos[i]
        end <- masspos[i + 2] - 1
        if (end - start <= 0)
            next
        subpeakl <- peaksL[start:end,]
        pos <- pos + 1

        ###Determining hte base of the signal to skip for a resolution.
        bw <- max(5 * massInter[i] * ppmGroup / (10 ^ 6),dmzGroup)
        if (bw > inter / 4) {
            warning(
                "Interval for grouping seems to short. Increase the inter parameter or modify the resolution parameter"
            )
        }
        den <-
            density(
                subpeakl[, "mz"],
                bw,
                from = massInter[i] -  bw,
                to = massInter[i + 2] +
                    bw,
                n = nPoints
            )

        ###Putting the close to 0 to 0
        maxy <- max(den$y)
        den$y[which(den$y <= bw * maxy)] = 0
        plim <- c(-1, 0, 0)
        oldMz <- previousMz
        previousMz <- NULL
        GraphPeakBegin <- NULL
        GraphPeakEnd <- NULL
        ###Peak are detected while density peaks remains in the selected interval.
        repeat {
            plim <- findLimDensity(den$y, plim[2] + 1, plim[3])

            if (plim[1] == plim[2])
                ###In this case the last peak is found.
                break
            selectedPGroup <- which(subpeakl[, 3] >= den$x[plim[1]] &
                                       subpeakl[, 3] <= den$x[plim[2]])
            if (length(selectedPGroup) == 0) {
                next
            }
            #       if (any(is.na(selectedPGroup))) {
            #         next
            #       }
            num_group <- num_group + 1

            ###Extending the maximal number of peaks if needed.
            if (num_group > nrow(groupmat)) {
                groupmat <-
                    rbind(groupmat, matrix(
                        nrow = nrow(groupmat),
                        ncol = ncol(groupmat)
                    ))
                listgroup <- c(listgroup, vector("list", length(listgroup)))
            }

            ###Checking if the group have not been added during th previous group.
            if (median(subpeakl[selectedPGroup, "mz"]) %in% oldMz ||
                (abs(median(subpeakl[selectedPGroup, "mz"]) - massInter[i]) < bw) ||
                (abs(median(subpeakl[selectedPGroup, "mz"]) - massInter[i + 2]) <
                 bw)) {
                num_group <- num_group - 1
                next
            }
            ###Checking the fraction of sampl found in the signal.
            posSample <- NULL
            if (solvar) {
                posSample <- selectedPGroup
            } else{
                ###Retaining the position of the peak on which sample have been detected.
                posSample <- selectedPGroup[which(subpeakl[selectedPGroup, "scanMin"] !=
                                                     1)]
            }
            ###CHecking that sample is detected in a high enough fraction of a class.
            vecLab <- object@classes[subpeakl[posSample, "sample"], 2]
            namesGroup <- names(tabclasses)
            found <- FALSE
            for (i in 1:length(namesGroup)) {
                numElem <- length(which(vecLab == namesGroup[i]))
                if (numElem >= fracGroup * tabclasses[i]) {
                    found <- TRUE
                    #groupmat[num_group, "class"] <- i
                    break
                }
            }

            ###If the peak is not reproductible, pass.
            if (!found) {
                num_group <- num_group - 1
                next
            }
            #       legList = c(
            #         "mzmed", "mzmin", "mzmax", "sizesamp", "scanmin", "scanmax","npeaks","solvent"
            #       )

            groupmat[num_group, "scanMin"] <-
                min(subpeakl[posSample, "scanMin"])
            groupmat[num_group, "scanMax"] <-
                max(subpeakl[posSample, "scanMax"])
            groupmat[num_group, "sizeSamp"] <-
                median(subpeakl[posSample, "scanMax"] - subpeakl[posSample, "scanMin"])
            groupmat[num_group, "mzMed"] <-
                median(subpeakl[selectedPGroup, "mz"])

            ##mz is stocked.
            previousMz <- c(previousMz, groupmat[num_group, 1])
            groupmat[num_group, c("mzMin", "mzMax")] <-
                range(c(subpeakl[selectedPGroup, "mzmin"], subpeakl[selectedPGroup, "mzmax"]))

            groupmat[num_group, "nPeaks"] <- length(posSample)
            groupmat[num_group, "timeShifted"] <- ifelse(sum(subpeakl[selectedPGroup, "timeShifted"])==0,0,1)

            groupmat[num_group, "meanSolvent"] <-
                mean(subpeakl[posSample, "solventIntensity"])
            if (boolpval) {
                groupmat[num_group, "signalOverSolventPvalueMean"] = median(subpeakl[posSample, "signalOverSolventPvalue"])
            }
            if (boolcor) {
                groupmat[num_group, "corSampPeakMean"] = mean(subpeakl[posSample, "corSampPeak"])
                groupmat[num_group, "corSampPeakSd"] = sd(subpeakl[posSample, "corSampPeak"])
            }
            if (boolss) {
                groupmat[num_group, "signalOverSolventRatioMean"] = mean(subpeakl[posSample,"signalOverSolventRatio"])
            }
            listgroup[[num_group]] <-
                sort(porder[(start:end)[posSample]])
            ###Adding the detected group to the current list.
            GraphPeakBegin <- c(GraphPeakBegin, plim[1])
            GraphPeakEnd <- c(GraphPeakEnd, plim[2])

        }
        oldMz <- previousMz
        if (sleep > 0 & !is.null(GraphPeakBegin)) {
            maxden <- max(den$y)
            maxint <- max(subpeakl[, "maxIntensity"])
            plot(den, main = paste(
                round(min(subpeakl[, 3]), 4),
                "-",
                round(max(subpeakl[, 3]), 4),
                " groups : ",
                length(GraphPeakBegin)
            ))
            for (gi in 1:length(GraphPeakBegin)) {
                lines(den$x[(GraphPeakBegin[gi]):(GraphPeakEnd[gi])],
                      den$y[(GraphPeakBegin[gi]):(GraphPeakEnd[gi])],
                      col =
                          630,
                      lwd = 1.5)
                abline(v = den$x[c(GraphPeakBegin[gi], GraphPeakEnd[gi])], col =
                           "darkgrey")
            }
            points(subpeakl[, "mz"],
                   subpeakl[, "maxIntensity"] * maxden / maxint,
                   col =
                       col_vec[start:end],
                   pch = 19)
            points(subpeakl[, "mz"],
                   subpeakl[, "maxIntensity"] * maxden / maxint,
                   col =
                       col_vec[start:end],
                   type = "h")

            ###Adding the solvent lvl.

            legend(
                "topright",
                c(as.character(uclasses), "peaks detected"),
                col = c(col_val, 630),
                lty =
                    c(lty_vec, 1),
                lwd = c(lwd_vec, 1.5)
            )
            Sys.sleep(sleep)
        }

    }
    #colnames(groupmat)<-c("mz","mzmin","mzmax","size","scan_min","scan_max","num_peaks","label")
    groupmat <- groupmat[1:num_group,]

    cat("\n", nrow(groupmat), " groups have been done .\n")
    object@group <- groupmat
    object@groupidx <- listgroup[1:num_group]
    object@step <- "Samples_grouped"
    object
})

getRawName <- function(filenames) {
    inte <- basename(filenames)
    fl <- strsplit(inte, split = ".", fixed = TRUE)
    return(unlist(lapply(fl, function(x) {
        x[1]
    })))
}


getUniqueIds <- function(ids){
	tid <- table(ids)
	pmul <- which(tid>1)

	lab <- names(tid)
	if(length(pmul)>0){
		for(i in 1:length(pmul)){
			posok <- which(ids == lab[pmul[i]])
			labvec <- paste(ids[posok],c('',paste('_',2:(tid[pmul[i]]),sep="")),sep = "")
			ids[posok] <- labvec
		}
	}
	ids
}

setGeneric("makeDataMatrix", function(object, ...)
    standardGeneric("makeDataMatrix"))
#' Construct the data matrix of a proFIAset object.
#'
#' Construct the data matrix of a proFIA set object, using the selected measure for the intensity,
#' area or maximum intensity. The choice of area or maximum intensity depends of your acquisition,
#' and your preferences.
#' @export
#' @param object A \code{\link{proFIAset}} object.
#' @param maxo Shall the intensity used to the area or the maximum
#' intensity
#' @aliases makeDataMatrix makeDataMatrix,proFIAset-method
#' @return A proFIAset object with the \code{dataMatrix} slot filled.
#' @seealso To obtain this data matrix see \code{\link{proFIAset}}.
#' @examples
#' if(require("plasFIA")){
#' #proFIAset object is loaded
#' data(plasSet)
#'
#' plasSet<-makeDataMatrix(plasSet)
#' plasSet
#' }

setMethod("makeDataMatrix", "proFIAset", function(object, maxo = FALSE) {
    if(!(object@step %in%c("Samples_grouped","Matrix_created","Fillpeaks"))){
        stop("Peaks needs to be grouped before creating the data Matrix, see ?group.FIA.")
    }
    unSample <- unique(peaks(object)[, "sample"])
    matResult  <-  matrix(NA_real_,
                         ncol  =  max(unSample),
                         nrow  =  nrow(object@group))
    for (i in 1:nrow(object@group)) {
        if (maxo) {
            matResult[i, object@peaks[object@groupidx[[i]], "sample"]] = object@peaks[object@groupidx[[i]], "maxIntensity"]
        } else{
            matResult[i, object@peaks[object@groupidx[[i]], "sample"]] = object@peaks[object@groupidx[[i]], "areaIntensity"]
        }
    }
    matResult <- matResult[, which(1:max(unSample) %in% unSample)]
    colnames(matResult) <- getUniqueIds(getRawName(object@classes[, 1]))
    rownames(matResult) <-
    	getUniqueIds(paste("M", sprintf("%0.4f", object@group[, 1]), sep = ""))
    object@dataMatrix <- matResult
    object@step <- "Matrix_created"
    object
})


makeSeparateColors <-
    function(rclasses,
             sizeSample = 2,
             sizeBlank = 2) {
        classes <- as.numeric(rclasses)
        tC <- table(classes)
        numClasses <- length(tC)
        inter <- seq(0, 1, length = ((sizeSample + sizeBlank) * numClasses))
        #cat(inter)
        vec_col <- rep(0, length(classes))
        for (i in 1:numClasses) {
            vec_col[which(classes == i)] = rainbow(tC[i], start = inter[(i - 1) * (sizeBlank +
                                                                                       sizeSample) + 1], end = inter[i * (sizeBlank + sizeSample) - sizeBlank])
        }
        vec_col
    }

makeSpaces <- function(rclasses,
                       sizeSample = 1,
                       sizeBlank = 4) {
    classes <- as.numeric(rclasses)
    tC <- table(classes)
    numClasses <- length(tC)
    tC <- c(0, cumsum(tC))
    #cat(inter)
    vec_space <- rep(0, length = (sizeSample * length(classes)))
    for (i in 1:(numClasses)) {
        vec_space[(tC[i] + 1):tC[i + 1]] <- (tC[i] + 1):tC[i + 1] + (i - 1) * sizeBlank
    }
    vec_space
}

setGeneric("findMzGroup", function(object, ...)
    standardGeneric("findMzGroup"))

#' find a group in a FIA experiment.
#'
#' Find a group corresponding to the given mass in a proFIAset
#' object, given a tolerance in ppm.The mz considered for a group
#' is the median fo the grouped signals between the various acquisition.
#' @export
#' @param object A proFIAset object.
#' @param mz A numeric vector of masses to be looked for.
#' @param tol The tolerance in ppm.
#' @param dmz The minimum tolerance in absolute mz compared to the tolerance in ppm.
#' @param closest Shall only the closest group be returned.
#' @return A vector of integer of the same length than mz giving the row
#' of the found group in the object group slot, or NA if the group is
#' not found.
#' @aliases findMzGroup findMzGroup,proFIAset-method
#' @return If closest is FALSE a list giving the index of the found group, if closest is TRUE a vector
#' giving the position. NA indicates that the mass signal have not bene found.
#' @seealso  You can visualize the group using \code{\link{plotFlowgrams}} function.
#' @examples
#' if(require("plasFIA")){
#' #proFIAset object is loaded
#' data(plasSet)
#'
#' #The table of spiked molecule is loaded
#' data(plasMols)
#'
#' #Mass to search and toolerance are defined
#' mass<-plasMols[22,"mass_M+H"]
#' tolppm <- 1
#'
#' plasSet<- makeDataMatrix(plasSet)
#'
#' index <- findMzGroup(plasSet,mass,tol=tolppm)
#'
#' plasMols[22,]
#' #We extract the corresponding group.
#' groupMatrix(plasSet)[index,]
#' }
setMethod("findMzGroup", "proFIAset", function(object, mz, tol,dmz = 0.005, closest = TRUE) {
    if(nrow(object@group)==0) stop("group needs to be created before using findMzGroup. See ?group.FIA.")

    vecRes <- vector(mode="list",length=length(mz))

    vectol <- tol*1e-6*object@group[,"mzMed"]
    vectol <- ifelse(vectol < dmz, dmz, vectol)
    bmin <- object@group[,"mzMed"]-vectol
    bmax <- object@group[,"mzMed"]+vectol
    for (i in 1:length(vecRes)) {
        pos <- which(bmin <= mz[i] &
        			 	bmax >= mz[i])
        if (length(pos) > 1 & closest) {
            best <- which.min(abs(object@group[pos, "mzMed"] - mz[i]) / mz[i])
            pos <- pos[best]
        }
        ##Checking if the found group got his med mz close to the peak.

        if (length(pos) == 0 | is.null(pos)) pos <- NA

       vecRes[[i]] <- pos
    }
    if(closest) vecRes <- unlist(vecRes)
    vecRes
})


setGeneric("impute.KNN_TN", function(object, ...)
	standardGeneric("impute.KNN_TN"))
#' Fill missing values in the peak table using K-nearest Neighbour for tuncated distributions.
#'
#' Impute the missing values in an FIA experiment using a Weighted
#' K-Nearest Neighbours on Truncated Distribution described by Jasmit S. Shah et al.
#' @export
#' @param object A proFIAset object.
#' @param k The number of neighbors considered, can be a fraction then in this case the
#' k will be taken for each class as the frac of the effective of the class. 3 at minima because comparison
#' is based on correlation.
#' @param classes how to handle imputation for different classes, if 'split', the classes
#' are taken separately, if 'unique', the imputation is done on the full data matrix.
#' @return A proFIAset object with the missing values imputated.
#' @aliases impute.KNN_TN impute.KNN_TN,proFIAset-method
#' @references Distribution based nearest neighbor imputation for truncated high dimensional data with applications to pre-clinical and
#'clinical metabolomics studies, J.S Shah 2017, BMC Bioinformatics.
#' @examples
#' if(require(plasFIA)){
#'     data(plasSet)
#'
#'     ###Reinitializing the data matrix
#'     plasSet<-makeDataMatrix(plasSet,maxo=FALSE)
#'     plasSet<-impute.KNN_TN(plasSet,2)
#' }

###We let the k being adaptative.
setMethod("impute.KNN_TN","proFIAset",function(object,k=0.6,classes=c("split","unique")){
	classes <- match.arg(classes)
	if(k!=round(k)){
		###
		if(k>1) stop("k shloud be an integer greater than 2 or a float inferior to 1.")
		rta <- table(object@classes[,1])
		prta <- which((as.numeric(rta)*k)<2)
		if(length(prta)!=0){
			stop(paste("The class(es)",names(rta)[prta],"got a too low number of samples to use with parameter",k))
		}
	}
	dm <- dataMatrix(object)
	dm[which(dm==0,arr.ind = TRUE)] <- NA
	dm <- log10(dm)
	allClasses <- unique(as.character(object@classes[,2]))
	b_imputation <- rep(TRUE,nrow(dm))
	resMatrix <- dm
	if(classes=="split" & length(allClasses)>1){
		####The matrix is splitted given the list sumbsample
		for(i in allClasses){
			psample <- which(as.character(object@classes[,2])==i)
			if(length(psample)<k) warnings(paste(
				"Using missing values imputation using",k,"while the class",i,"contains only",length(psample),"samples.\n"
			))
			###
			i_k <- ifelse(k==round(k),k,ceiling(k*length(psample)))
			temp_res <- imputeKNN(dm[,psample,drop=FALSE],k = i_k,distance="truncation")
			resMatrix[,psample] <- temp_res$imputedData

			cat(length(b_imputation),"vs",length(temp_res$problems))
			b_imputation <- b_imputation&temp_res$problems
		}

	}else{
		psample <- 1:ncol(dm)
		i_k <- ifelse(k==round(k),k,ceiling(k*length(psample)))
		temp_res <- imputeKNN(dm[,psample,drop=FALSE],k = i_k,distance="truncation")
		resMatrix[,psample] <- temp_res$imputedData
		b_imputation <- b_imputation&temp_res$problems
	}

	object@dataMatrix <- 10^resMatrix

	cnames <- colnames(object@group)
	object@group <- cbind(object@group,as.numeric(b_imputation))
	colnames(object@group) <- c(cnames,"imputationOk")
	return(object)
})



setGeneric("impute.randomForest", function(object, ...)
	standardGeneric("impute.randomForest"))
#' Fill missing values in the peak table using random forest.
#'
#' Impute the missing values in an FIA experiment using a random forest implemented
#' in the missForest package.
#' @export
#' @param object A proFIAset object.
#' @param parallel Shall parallelism be used.
#' @param ... supplementary arguements to be passed to missForest function.
#' @return A proFIAset object with the missing values imputated.
#' @aliases impute.randomForest impute.randomForest,proFIAset-method
#' @references Stekhoven, D.J. and Buehlmann, P. (2012), 'MissForest - nonparametric missing value imputation for mixed-type data',
#'  Bioinformatics, 28(1) 2012, 112-118, doi: 10.1093/bioinformatics/btr597
#' @examples
#' if(require(plasFIA)){
#'     data(plasSet)
#'     ###Reinitializing the data matrix
#'     plasSet<-makeDataMatrix(plasSet,maxo=FALSE)
#'     plasSet<-impute.randomForest(plasSet)
#' }
setMethod("impute.randomForest","proFIAset",function(object,parallel=FALSE,...){
	dm <- dataMatrix(object)
	if(! parallel){
		parallel <- "no"
	}else{
		parallel <- "variables"
	}

	dm <- missForest(dm,parallelize=parallel,...)$ximp

	object@dataMatrix <- dm
	return(object)
})



setGeneric("impute.FIA", function(object, ...)
	standardGeneric("impute.FIA"))
#' Fill missing values using the provided method.
#'
#' Impute the missing values in an FIA experiment using a Weighted
#' K-Nearest Neighbours on Truncated Distribution implemented in \code{\link{impute.KNN_TN}} or by random forest using the \code{\link{impute.randomForest}} function.
#' @export
#' @param object A proFIAset object.
#' @param method The method to be used for missing value imputation.
#' @param ... Arguments furnished to the imputation method. No argument is required for \code{\link{impute.randomForest}} and the number of neighbours should be furnished using the \code{\link{impute.KNN_TN}} function.
#' @aliases impute.FIA impute.FIA,proFIAset-method
#' @examples
#' if(require(plasFIA)){
#'     data(plasSet)
#'
#'     ###Reinitializing the data matrix an using KNN
#'     plasSet<-makeDataMatrix(plasSet,maxo=FALSE)
#'     plasSet<-impute.FIA(plasSet,method="KNN_TN",k=2)
#'
#'     ###Reinitializing the data matrix and using randomForest
#'     plasSet<-makeDataMatrix(plasSet,maxo=FALSE)
#'     plasSet<-impute.FIA(plasSet,method="randomForest")
#' }

###We let the k being adaptative.
setMethod("impute.FIA","proFIAset",function(object,method=c("KNN_TN","randomForest"),...){
	if(object@step=="Fillpeaks") stop("Missing value have already been imputed, redo it, use makeDataMatrix then
									  impute missing values.")
	method <- match.arg(method)
	if(method=="KNN_TN"){
		object <- impute.KNN_TN(object,...)
	}else if(method == "randomForest"){
		object <- impute.randomForest(object,...)
	}
	object@step <- "Fillpeaks"
	return(object)
})







setGeneric("peaksGroup", function(object, ...)
    standardGeneric("peaksGroup"))
#' Return the peaks corresponding to a group.
#'
#' Return the peaks corresponding ot a group given by his index.
#' @export
#' @param object A proFIAset object.
#' @param index A numeric vector r giving the group to be returned. NA are ignored.
#' @return The peaks in the given group, see \code{\link{proFIAset-class}}.
#' @aliases peaksGroup peaksGroup,proFIAset-method
#' @examples
#' if(require(plasFIA)){
#'     data(plasSet)
#'     data(plasMols)
#'
#'     #finding the molecules of plasMols
#'     vmatch<-findMzGroup(plasSet,mz=plasMols[,"mass_M+H"],tol=5)
#'
#'     mol_peaks<-peaksGroup(plasSet,index=vmatch)
#'     head(mol_peaks)
#' }
setMethod("peaksGroup", "proFIAset", function(object,
                                              index = NULL) {
    if (!is.numeric(index))
        stop("Index need to be an integer or a vector of numeric.")
    posna <- which(is.na(index))
    if(length(posna)!=0){
        index <- index[-posna]
    }
    if (length(index)==0)
        stop("Index is of size 0, or only NA are provided.")
    matGroup <- object@peaks[unlist(object@groupidx[index]),,drop=FALSE]
    return(matGroup)
})

colorKeys <- function(i) {
    if (i == 0)
        return("white")
    switch(EXPR = i,
           "red",
           "green3",
           "blue",
           "cyan",
           "orange",
           "black")
}


layoutMatrix <- function(n) {
    if (n == 1) {
        return(matrix(c(1)))
    }
    if (n == 2) {
        return(matrix(c(1, 2), nrow = (2)))
    }
    if (n == 3) {
        return(matrix(c(1, 2, 3), nrow = (3)))
    }
    if (n == 4) {
        return(matrix(c(1, 2, 3, 4), nrow = (2), byrow = TRUE))
    }
    if (n == 5) {
        return(matrix(c(1, 2, 3, 4, 5, 6), nrow = (2), byrow = TRUE))
    }
    if (n == 6) {
        return(matrix(c(1, 2, 3, 4, 5, 6), nrow = (2), byrow = TRUE))
    }
    if (n == 7) {
        return(matrix(
            c(1, 2, 3, 4, 5, 6, 7, 8, 9),
            nrow = (3),
            byrow = TRUE
        ))
    }
    if (n == 8) {
        return(matrix(
            c(1, 2, 3, 4, 5, 6, 7, 8),
            nrow = (2),
            byrow = TRUE
        ))
    }
    if (n == 9) {
        return(matrix(
            c(1, 2, 3, 4, 5, 6, 7, 8, 9),
            nrow = (3),
            byrow = TRUE
        ))
    }
    if (n == 10) {
        return(matrix(
            c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
            nrow = (4),
            byrow = TRUE
        ))
    }
    if (n == 11) {
        return(matrix(
            c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
            nrow = (4),
            byrow = TRUE
        ))
    }
    if (n == 12) {
        return(matrix(
            c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
            nrow = (4),
            byrow = TRUE
        ))
    }
}

toVoid <- function(n) {
    if (n == 5)
        return(5)
    if (n == 7)
        return(c(6, 8))
    if (n == 10)
        return(c(10, 12))
    if (n == 11)
        return(c(12))
    return(0)
}

setGeneric("exportVariableMetadata", function(object, ...)
    standardGeneric("exportVariableMetadata"))
#' Export variable metadata.
#'
#' Export the variable metadata of an experiment, to
#' be used for statistical analysis.
#'
#' @export
#' @param object A proFIAset object.
#' @param filename If not NULL the result will be written
#' in filename
#' @return A dataframe with the following columns :
#' \itemize{
#'     \item variableID an ID similar to the one of the peak table.
#'     \item mzMed the median value of group in the m/z dimension.
#'     \item mzMin the minimum value of the group in the m/z dimension.
#'     \item mzMax the maximum value of the group in the m/z dimension.
#'     \item scanMin the first scan on which the signal is detected.
#'     \item scanMax the last scan on which the signal is detected.
#'     \item nPeaks The number of peaks grouped in a group.
#'     \item meanSolvent The mean of solvent in the acquisition.
#'     \item signalOverSolventPvalue The mean p-value of the group.
#'     \item corMean The mean of the matrix effect indicator.
#'     \item SigSolMean The mean of ratio of the signal max
#'     intensity on the solvent max intensity.
#' }
#' @aliases exportVariableMetadata exportVariableMetadata,proFIAset-method
#' @examples
#' if(require(plasFIA)){
#'   data(plasSet)
#'   vtab<-exportVariableMetadata(plasSet)
#'   head(vtab)
#' }
setMethod("exportVariableMetadata", "proFIAset", function(object, filename =
                                                              NULL) {
    toExport = data.frame(object@group, stringsAsFactors = FALSE)
    toExport["variableID"] = getUniqueIds(row.names(object@dataMatrix))
    if (!is.null(filename)) {
        write.table(
            toExport,
            file = filename,
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t"
        )
    }
    return(invisible(toExport))
})

setGeneric("exportSampleMetadata", function(object, ...)
    standardGeneric("exportSampleMetadata"))
#' Export samples metadata.
#'
#' Export the samples metadata of an experiment, to
#' be used for statistical analysis.
#'
#' @export
#' @param object A proFIAset object.
#' @param filename If not NULL the result will be written
#' in filename
#' @return A dataframe with the following columns :
#' \itemize{
#'     \item sampleID an ID similar to the one of the peak table.
#'     \item class the group of the sample.
#' }
#' @aliases exportSampleMetadata exportSampleMetadata,proFIAset-method
#' @examples
#' if(require(plasFIA)){
#'    data(plasSet)
#'    tsample<-exportSampleMetadata(plasSet)
#'    head(tsample)
#' }
setMethod("exportSampleMetadata", "proFIAset", function(object, filename =
                                                            NULL) {
    if(nrow(object@dataMatrix)==0){
        stop("dataMatrix needs to be created before exporting sampleMetadata.")
    }
    toExport <- data.frame(
        sampleID = colnames(object@dataMatrix),
        class = as.character(object@classes[, 2]),
        stringsAsFactors =
            FALSE
    )
    if (!is.null(filename)) {
        write.table(
            toExport,
            file = filename,
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t"
        )
    }
    return(invisible(toExport))
})


setGeneric("plotFlowgrams", function(object, ...)
    standardGeneric("plotFlowgrams"))

#' Plot raw temporal profiles of the selected group.
#'
#' Plot raw temporal profiles from \code{\link{proFIAset}} object
#' corresponding to one or more molecules. The function will priorize
#' index, only using mz if index is set to \code{NULL}. As this require to come
#' back to the raw data, this can take some times, an we don't recommend using it
#' on more than 30 files. It is better to choose 30 files in the proFIAset object.
#'
#' @export
#' @param object A proFIAset object.
#' @param index The index of the group to be plotted.
#' @param mz An mz value to be looked for only used if index is null.
#' the research use the \code{\link{findMzGroup}} function.
#' @param subsample A subset of sample to be plotted.
#' @param ppm The tolerance for the research if mz is provided.
#' @param margin An area outer the EICs mz range on which the EIC may be extended.
#' @param posleg The position of the legend on the figure. See \code{\link[graphics]{legend}}.
#' @param title An optional vector of title for the plot. Need to be of the same
#' @param scaled Shall all the EIC be put on the same scale with maximum to 1.
#' @param area Shall the detectged area be plotted using transparency.
#' @param ...  Supplementary graphical parameters to be passed to lines.
#' length than index.
#' @aliases plotFlowgrams plotFlowgrams,proFIAset-method
#' @return  No returned value
#' @examples
#' if(require(plasFIA)){
#'     data(plasMols)
#'     data(plasSet)
#'     plotFlowgrams(plasSet,mz=plasMols[7,"mass_M+H"])
#' }
setMethod("plotFlowgrams", "proFIAset", function(object,
                                            index = NULL,
                                            mz = NULL,
                                            subsample = NULL,
                                            ppm = 5,
                                            margin = 0.0002,
                                            posleg=c("topright","bottomright", "bottom",
                                                     "bottomleft", "left", "topleft",
                                                     "top", "right", "center"),title=NULL,
                                            scaled=FALSE,area=FALSE,...) {
    posleg <- match.arg(posleg)
    if(!(object@step %in% c("Samples_grouped","Matrix_created","Fillpeaks"))){
        stop("Peaks needs to be grouped before plotting the EICs, see ?group.FIA.")
    }

    if (is.null(index)) {
        if (is.null(mz))
            stop("Index or mz need to be provided")
        index <- sapply(mz, findMzGroup, object = object, tol = ppm)
    }
    if(!is.null(title)){
    if(length(title)!=length(index)){
        stop("Title and index needs to have the same length, they are respectively ",
             length(title)," and ",length(index),"\n")
    }
    }
    tclasses <- object@classes
    matpres <- matrix(0,
                     nrow = length(index),
                     ncol = nrow(object@classes))

    for (i in 1:length(index)) {
        matpres[i, object@peaks[object@groupidx[[index[i]]],,drop=FALSE][, "sample"]] = 1
    }
    if (!is.null(subsample)) {
        matpres <- matpres[, subsample, drop = FALSE]
    } else{
        subsample <- 1:nrow(object@classes)
    }
    # print(matmz)
    # print(matpres)
    vall <- sapply(subsample, function(x, paths, peaks, idx, to_idx) {
        xraw  <-  xcmsRaw(paths[x])
        rEIC  <-  vector(mode="list",length=length(idx))
        for (j in seq_along(to_idx))
        {
          sel_idx <- idx[[to_idx[j]]]
          sel_idx <- sel_idx[peaks[sel_idx,"sample"]==x]
          if(length(sel_idx)==0){
            rEIC[[j]] <- NA
            next
          }else{
            mzrange <- c(min(peaks[sel_idx,"mzmin"]),max(peaks[sel_idx,"mzmax"]))
            scanlim <- c(min(peaks[sel_idx,"scanMin"]),max(peaks[sel_idx,"scanMax"]))
            rEIC[[j]] <- list(intensity=rawEIC(xraw, mzrange = c(mzrange[1], mzrange[2]))$intensity)
            if(scaled) rEIC[[j]][[1]] <- rEIC[[j]][[1]]/max(rEIC[[j]][[1]])
            sscan <- rep(FALSE,length(xraw@scantime))
            sscan[scanlim[1]:scanlim[2]] <- TRUE
            rEIC[[j]][[2]] <-sscan
          }
        }
        rEIC$scantime <- xraw@scantime
        rEIC
    }, paths = object@classes[, 1], peaks = object@peaks,
    to_idx = index,idx=object@groupidx, simplify =FALSE)
    #return(vall)
    ###Now we have a list of list[[list[[]]]] with expriment[acquisition]
    colvec <- NULL
    colleg <- NULL
    ileg <- NULL
    if (length(unique(object@classes[subsample, 2])) == 1) {
        colvec <- rainbow(length(subsample))
        colleg <- rainbow(length(subsample))
    } else{
        colvec  <-  makeSeparateColors(object@classes[, 2])[subsample]
        colleg <- unique(object@classes[subsample, 2])
        ileg <- match(unique(object@classes[subsample, 2]),object@classes[subsample, 2])
        ileg <- floor((ileg+c((ileg[-1]-1),length(subsample)))/2)
    }
    ###Now we create the area color if necessary
    colarea <- NULL
    if(area){
    	colarea <- col2rgb(colvec)
    	# colvec <- rep("black",length(colvec))
    	colarea <- apply(colarea,2,function(x){
    		rgb(x[1],x[2],x[3],80,maxColorValue = 255)
    	})
    }



    maxx  <-  max(unlist(sapply(vall, function(x) {
        x$scantime
    })))
    for (i in 1:length(index)) {
        pok  <-  which(matpres[i,]  !=  0)
        if (length(pok)  ==  0)
            next
        ###Getting the max
        vecval  <-  lapply(vall, "[[", i)
        maxy  <-  max(unlist(vecval[pok]))
        mtitle <- NULL
        if(is.null(title)){
            mtitle  <-  paste("Group :", index[i], sprintf("%0.4f", object@group[index[i],"mzMed"]), "m/z")
        }else{
            mtitle <- title[i]
        }


        plot(
            NULL,
            xlab  =  "Time (s)",
            ylab  =  "Intensity",
            main  =  mtitle,
            xlim  =  c(0, maxx),
            ylim  =  c(0, maxy)
        )
        for (j in 1:length(pok)) {

            lines(vall[[pok[j]]]$scantime, vall[[pok[j]]][[i]][[1]], col = colvec[pok[j]])
        	if(area){

        	  sx <- vall[[pok[j]]]$scantime[vall[[pok[j]]][[i]][[2]]]
        		sx <- c(sx,rev(sx))
        		sy <- vall[[pok[j]]][[i]][[1]]
        		sy <- sy[vall[[pok[j]]][[i]][[2]]]
        		sy <- c(sy,rep(0,length(sy)))
        		polygon(sx,sy,col=colarea[pok[j]])
        	}
        }
        ###Adding the legend.
        if (length(unique(object@classes[subsample, 2])) == 1) {
            legend(posleg,as.character(proFIA:::getRawName(object@classes[subsample, 1]))[pok],lty=rep(1,length(subsample))[pok],col=colvec[pok])
        } else{
            legend(posleg,as.character(unique(object@classes[subsample, 2])),lty=rep(1,length(unique(object@classes[subsample, 2]))),col =
                       colvec[ileg])
        }


    }
})

setGeneric("plotSamplePeaks", function(object, ...)
    standardGeneric("plotSamplePeaks"))
#' Plot the sample peak of a proFIAset object.
#'
#' Plot the sample peaks determined by regression for each sample.
#' If you want to check the various regression performed by
#' \emph{proFIA} in an individual sample we recommend you to use the
#' \code{\link{getInjectionPeak}}. A sample peak highly different form the
#' other may indicate a problem of the regression process, or an injection issue
#' in the acquisition.
#'
#' @export
#' @param object A proFIAset object.
#' @param subsample The subset of sample on which the sample may be plotted.
#' If it is numeric it will be viewed as sample row in the classes table,
#' if it is a character it will be viewed as a factor.
#' @param diagPlotL Boolean used by the usmmary plot function. Keep to FALSE in the general case.
#' @param ... SUpplementary arguments which will be passed to the \link{lines}
#' function.
#' @aliases plotSamplePeaks plotSamplePeaks,proFIAset-method
#' @return No value returned.
#' @examples
#' if(require(plasFIA)){
#'     data(plasSet)
#'     plotSamplePeaks(plasSet)
#' }
setMethod("plotSamplePeaks", "proFIAset", function(object, subsample = NULL, diagPlotL = FALSE, ...) {
	if (is.null(subsample)) {
		subsample <- 1:nrow(object@classes)
	} else if (is.character(subsample)) {
		message("subsample is a character trying to match it as a factor.")
		subsample <- which(object@classes[, 2] %in% subsample)
		if (length(subsample) == 0)
			stop(
				"Wrong subsample provided, subsample shall be a number
				or a class of the proFIAset object."
			)
	} else if (is.numeric(subsample)) {
		if (any(subsample > nrow(object@classes)))
			stop("Subsample shall be less than the number of  sample.")
	}
	toPlot <- object@injectionPeaks[subsample]
	vbeginning <- object@injectionScan[subsample]
	vl <- sapply(toPlot, length)
	rangex <- c(0, max(vbeginning + vl - 1))
	colvec <- makeSeparateColors(object@classes[subsample, 2])
	vleg <- table(object@classes[subsample, 2])
	## mtitle <- "Sample Injection Peaks"
	mtitle <- "Sample peak models"
	vcol <- NULL
	if (length(vleg) == 1) {
		## mtitle <- paste(mtitle, "class", object@classes[subsample[1], 2])
	} else{
		vleg <- cumsum(as.numeric(vleg))
		print(vleg)
		vcol <- colvec[vleg]
		vleg <- object@classes[subsample[vleg], 2]
	}
	if(diagPlotL) {
		omar <- par("mar")
		par(mar=c(2.6,2.6,2.1,1.1))
		plot(
			1:2,
			xlab = "",
			ylab = "Injection Peak",
			main = mtitle,
			xlim = rangex,
			ylim = c(0, 1),
			type="n"
		)
		mtext("Time (s)", side = 1, line = 2.5, cex = 0.7)
	}else{
	###To be usre to get the good plot.
	plot(
		1:2,
		 xlab = "Scan number",
		 ylab = "Injection Peak",
		main = mtitle,
		xlim = rangex,
		ylim = c(0, 1),
		type="n"
	)
	}
	if(diagPlotL) {
		par(mar = omar)
	}

	## mtext(mtitle, line = 1, cex = 0.8)
	for (i in 1:length(subsample)) {
		lines((vbeginning[i]):(vbeginning[i] + vl[i] - 1),
			  toPlot[[i]],
			  col = colvec[i],
			  ...)
	}
	if (length(vleg) != 1) {
		legend("topright",
			   as.character(vleg),
			   col = vcol,
			   lty = 1)
	}
	return(invisible(object@dataMatrix))
})



setGeneric("exportDataMatrix", function(object, ...)
    standardGeneric("exportDataMatrix"))
#' Export data matrix.
#'
#' Export the data matrix from a \code{\link{proFIAset}} object, to
#' be used for statistical analysis.
#'
#' @export
#' @param object A proFIAset object.
#' @param filename If not NULL the result will be written
#' in filename as a tabular separated values file.
#' @return A matrix with dimension samples x variables.
#' @aliases exportDataMatrix exportDataMatrix,proFIAset-method
#' @examples
#' if(require(plasFIA)){
#'     data(plasSet)
#'     dm<-exportDataMatrix(plasSet)
#'     head(dm)
#' }

setMethod("exportDataMatrix", "proFIAset", function(object, filename = NULL) {
    if (!is.null(filename)) {
        write.table(
            object@dataMatrix,
            file = filename,
            row.names = TRUE,
            col.names = TRUE,
            sep =
                "\t"
        )
    }
    return(invisible(object@dataMatrix))
})



#' Wrapper function for the full FIA analysis workflow.
#'
#' Perform the 4 steps pro fia workflow including :
#' \itemize{
#'  \item noise estimation. Noise is estimated.
#'  \item bands filtering. Bands are filtered using the \code{\link{findFIASignal}} function.
#'  \item peak grouping. Signals from different acquisition are grouped using \code{\link{group.FIA}} function.
#'  \item missing values imputations. Missing values are imputated using \code{\link{impute.KNN_TN}} function.
#' }
#' Minimal options to launch the workflow are provided, neithertheless if finer option tuning are
#' necessary, launching the workflow function by function is strongly advised.
#' @export
#' @param path The path to the directory of acquistion.
#' @param ppm The tolerance for deviations in m/z between scan in ppm \code{\link{findFIASignal}}
#' @param dmz The minimum tolerance for deviations in m/z between scan in Dalton \code{\link{findFIASignal}}
#' @param fracGroup The fraction of smaple form a class necessary to make a group.
#' @param ppmGroup A ppm tolerance to group signals between samples \code{\link{group.FIA}}.
#' @param dmzGroup The minimum tolerance in Dalton to group signals between samples \code{\link{group.FIA}}.
#' @param parallel A boolean indicating if parallelism is supposed to be used.
#' @param bpparam A BiocParallelParam object to be passed i BiocParallel is used.
#' @param noiseEstimation A boolean indicating in noise need to be estimated.
#' @param maxo Should the maximum intensity be used over the peak area.
#' @param SNT A value giving the SNT thrshold, used only if \code{noiseEstimation} is FALSE.
#' @param imputation The method to use for imputation 'randomForest' (\code{\link{impute.randomForest}})  or 'KNN_TN' (\code{\link{impute.KNN_TN}}). Put 'None' if  no imputation should be done.
#' @param ... Supplementary arguments to be passed to the imputation method, see \code{\link{impute.FIA}} for detail.
#' @return A filled proFIAset object ready for exportation.
#' @aliases analyzeAcquisitionFIA
#' @examples
#' if(require(plasFIA)){
#'     path<-system.file(package="plasFIA","mzML")
#'
#'     #Defining parameters for Orbitrap fusion.
#'     ppm<-2
#'     dmz <- 0.0005
#'     ppm_group<-1
#'     dmz_group <- 0.0003
#'     paral<-FALSE
#'     frac_group<-0.2
#'     k<-2
#'     maxo<-FALSE
#'     vimputation <- "randomForest"
#'     \dontrun{plasSet<-analyzeAcquisitionFIA(path,ppm=ppm,dmz=dmz,imputation=vimputation,fracGroup=frac_group,ppmGroup=ppm_group,dmzGroup=dmz_group,k=k,maxo=maxo,parallel=paral)}
#' }

analyzeAcquisitionFIA <-
    function(path,
             ppm,
    		 dmz = 0.001,
             fracGroup = 0.5,
             ppmGroup = NULL,
    		 dmzGroup = 0.0005,
             parallel = FALSE,
             bpparam = NULL,
             noiseEstimation = TRUE,
             SNT = NULL,
             maxo = FALSE,
    		 imputation = c("KNN_TN","randomForest","None"),
             ...) {
        if (!dir.exists(path)) {
            stop(
                "The specified must be a valid directory leading to an FIA experiment. Use findSignal
                if you wants to process in a single file."
            )
        }
    	imputation <- match.arg(imputation)
        if (ppm > 20) {
            warning(
                "proFIA have been designed for very high resolution data, choose a lower ppm
                value."
            )
        }
        ####Checking if noise is supposed ot be used.
        if (!noiseEstimation) {
            if (is.null(SNT)) {
                stop(
                    "No noise estimation and no SNT, peaks can't be filtered, please furnish one of the two."
                )
            } else{
                message("No noise estimation pure SNT threshold will be used.")
                maxo <- TRUE
            }
        } else{
            message("Noise estimation will be used.")
        }
        message("Step 1 : Noise Model Estimation and peak detection.")
        pset <- proFIAset(
            path,
            ppm,
            dmz=dmz,
            parallel = parallel,
            BPPARAM = bpparam,
            noiseEstimation = TRUE
        )
        if (is.null(ppmGroup)) {
            warning("ppmGroup is null default value is half of the ppm parameter value.")
            ppmGroup <- ppm / 2
        }
        mg <- "Step 2 : Peak Grouping."
        mg <- ifelse(maxo,
                    paste(mg, "Max Intensity will be used."),
                    paste(mg, "Area will be used."))
        message(mg)
        pset <- group.FIA(pset, ppmGroup, fracGroup = fracGroup, dmzGroup = dmzGroup)
        pset <- makeDataMatrix(pset,maxo=maxo)
        message("Step 3 : Missing values imputation.")
        if(imputation!="None"){
            pset <- impute.FIA(pset,method=imputation,...)
        }
        message(paste("Processing finished."))
        plot(pset)
        pset
        }

setMethod("cut", "proFIAset", function(x, subsample) {
    if (any(!(subsample %in% 1:nrow(x@classes)))) {
        stop(
            "Sample needs to be integer corresponding to a sample of the proFIAset object. To see
            a list of the sample use the getPhenoClasses method."
        )
    }
    x@classes <- x@classes[subsample,]
    dico <- list()
    for (i in 1:length(subsample)) {
        dico[subsample[i]] = i
    }
    x@peaks <- x@peaks[which(x@peaks[, "sample"] %in% subsample),]
    x@peaks[, "sample"] <- unlist(dico[x@peaks[, "sample"]])
    x@dataMatrix <- matrix(0, nrow = 0, ncol = 0)
    x@group <- matrix(0, nrow = 0, ncol = 0)
    x@groupidx <- list()
    x
    })



# #####BETA THose feature won(t be exported but will be disponible in near future)
# setGeneric("compareClasses", function(object, ...)
#     standardGeneric("compareClasses"))
#
# #' t-test of variables to check mean egality.
# #'
# #'Compare the mean of all the variables present in the
# #'data between mulltiples classes.
# #'
# #' @param object A proFIAset object.
# #' @param subsample If chracter a subset of the class to be plotted
# #' if numeric a subset of the sample to be plotted.
# #' @return A filled proFIAset object witht he pvalue and the ofld change for the
# #' classes. If a 0 remaines, the majority of the signal are im
# #' @aliases compareClasses compareClasses,proFIAset-method
# #' @examples
# #' print("examples to be put here")
# setMethod("compareClasses", "proFIAset", function(object, subsample = NULL) {
#     numGroup <- length(unique(object@classes))
#     if (is.null(subsample))
#         subsample <- 1:nrow(object@classes)
#     if (numGroup == 1)
#         stop("Multiple classes needs to be provided.")
#     if (numGroup > 2)
#         message("More than 2 each group will be tested against each other.")
#
#
#     #####Generating the list of sample to do
#     uSample <- unique(object@classes[, 2])
#     toCompare <- combn(uSample, 2)
#     toCompare <- toCompare[, which(toCompare[1,] != toCompare[2,]), drop =
#                               FALSE]
#     print(toCompare)
#     nameVec <- paste(toCompare[1,], toCompare[2,], sep = "_")
#     lRes <- list()
#     dMat <- dataMatrix(object)
#     for (i in 1:ncol(toCompare)) {
#         vx <- which(object@classes[, 2] == toCompare[1, i])
#         vy <- which(object@classes[, 2] == toCompare[2, i])
#         nx <- length(vx)
#         ny <- length(vy)
#         pvallab <- paste("p.value", nameVec[i], sep = "_")
#         abslab <- paste("p.value", nameVec[i], sep = "_")
#         ###Removing the 0 if there is some.
#         submat <- dataMatrix(object)[, c(vx, vy)]
#
#         ##Calculating the p-value.
#         pvalv <- apply(submat, 1, function(x, nx, ny) {
#             px0 <- which(x[1:nx] == 0)
#             if (length(px0) != 0) {
#                 nx <- nx - length(px0)
#                 x <- x[-px0]
#             }
#             py0 <- which(x[(nx + 1):ny] == 0)
#             if (length(py0) != 0) {
#                 ny <- ny - length(py0)
#                 x <- x[-py0]
#             }
#             if (nx == 0 | ny == 0)
#                 return(NA)
#             vt <- t.test(x[1:nx], x[(nx + 1):(nx + ny)])
#             c(vt[["p.value"]])
#         }, nx = nx, ny = ny)
#
#         ##Calculating the fold change.
#         foldv <- apply(submat, 1, function(x, nx, ny) {
#             px0 <- which(x[1:nx] == 0)
#             if (length(px0) != 0) {
#                 nx <- nx - length(px0)
#                 x <- x[-px0]
#             }
#             py0 <- which(x[(nx + 1):ny] == 0)
#             if (length(py0) != 0) {
#                 ny <- ny - length(py0)
#                 x <- x[-py0]
#             }
#             if (nx == 0)
#                 return(0)
#             if (ny == 0)
#                 return(Inf)
#             fold <- mean(x[1:nx]) / mean(x[(nx + 1):(nx + ny)])
#             fold
#         }, nx = nx, ny = ny)
#
#         #Cal
#
#
#         lRes[[i]] = pvalv
#     }
#     lRes <- do.call("cbind", lRes)
#     vgroup <- object@group
#     ngroup <- colnames(vgroup)
#     vgroup <- cbind(vgroup, lRes)
#     colnames(vgroup) <- c(ngroup, nameVec)
#     object@group <- vgroup
#     object
# })


plotVoidRaw <-
    function(mz, int, scanindex, scantime, title = NULL,bandl = NULL,legend=TRUE, ...) {
        vnum <- diff(c(scanindex, length(mz)))
        vtime <- rep(scantime, vnum)
        int <- log10(int)
        if (is.null(title))
            title <- paste("Raw data")
        xlim <- NULL
        ylim <- NULL
        ###Checking the xlim and ylim arg
        largs <- list(...)
        if ("xlim" %in% names(largs)) {
            xlim <- largs[["xlim"]]
            pok <- which(vtime <= xlim[2] & vtime >= xlim[1])
            if (length(pok) > 0) {
                vtime <- vtime[pok]
                mz <- mz[pok]
                int <- int[pok]
            }
        } else{
            xlim <- c(0, max(scantime))
        }

        if ("ylim" %in% names(largs)) {
            ylim <- largs[["ylim"]]
            pok <- which(mz <= ylim[2] & mz >= ylim[1])
            if (length(pok) > 0) {
                vtime <- vtime[pok]
                mz <- mz[pok]
                int <- int[pok]
            }
        }
        size <- NULL
        if ("size" %in% names(largs)) {
            size <- largs[["size"]]
        } else{
            size <- 0.4
        }

        ##Getting the color
        vpal <- colorRampPalette(c("green", "orange", "red"))

        ###Making the color vec
        rInt <- range(int)
        rInt[1] <- floor(rInt[1] - 0.01)
        rInt[2] <- ceiling(rInt[2] + 0.01)
        rcol <- seq(rInt[1], rInt[2], length = 50)
        vcol <- vpal(50)[.bincode(int, rcol)]
        if(legend){
        layout(matrix(1:2, ncol = 2),
               widths = c(3, 1),
               heights = c(1, 1))
        }
        legend_image <- as.raster(matrix(rev(vpal(50)), ncol = 1))
        plot(
            vtime,
            mz,
            main = title,
            xlab = "Time (s)",
            ylab = "m/z",
            cex = size,
            col = vcol,
            xlim = xlim,
            ylim = ylim,
            pch = 19
        )
        if(!is.null(bandl)){
        	abline(h=bandl,col="purple",lwd=0.5,lty=2)
        }
        if(legend){
        plot(
            c(0, 2),
            c(rInt[1], rInt[2]),
            type = 'n',
            axes = FALSE
            ,
            xlab = '',
            ylab = '',
            main = 'log10(Intensity)',
            cex.main = 0.9
        )
        text(
            x = 1.5,
            y = seq(rInt[1], rInt[2], l = rInt[2] - rInt[1] + 1),
            labels = seq(rInt[1], rInt[2], l = rInt[2] - rInt[1] + 1)
        )
        rasterImage(legend_image, 0, rInt[1], 1, rInt[2])
        layout(1)
        }
    }

# setGeneric("plotRaw", function(object, ...)
#     standardGeneric("plotRaw"))

#' Plotting of raw data
#'
#'Plot the raw data from a proFIAset object,the
#'the type of plot determines if the full raw data
#'needs to be plotted or only the data conressponding to
#'the detected peaks needs to be plottes. The path in the classes
#'table of the proFIAset object needs to be correct.
#'
#' @export
#' @param object A proFIAset object.
#' @param type \code{"raw"} indicate that raw data needs to be plotted
#' and \code{"peak"} indicate that only the filtered signals will be plotted.
#' @param sample The number of the sample in the object classes table
#' to be plotted. The classes table can be taken using the phenoClasses function.
#' @param bandl A boolean shall a line be added on the detected bands. Only used in type is set to \code{"peak"}
#' @param legend Shall the legend be plotted along the graph
#' @param ... xlim,ylim and size to be passed to plot functions.
#' @return No value is returned.
#' @aliases plotRaw plotRaw,proFIAset-method
#' @examples
#' if(require("plasFIA")){
#'     data(plasSet)
#'
#'     #Visualising the raw data
#'     plotRaw(plasSet,type="raw",ylim=c(215.9,216.2),sample=4)
#'
#'     #Plotting the filtered signals only.
#'     plotRaw(plasSet,type="peaks",ylim=c(215.9,216.2),sample=4)
#' }
setMethod("plotRaw", "proFIAset", function(object,
                                           type = c("raw", "peaks"),
										   sample = NULL,
										   bandl = TRUE,
                                           legend=TRUE,
                                           ...) {
    type <- match.arg(type)
    if (!(sample %in% 1:nrow(object@classes))) {
        stop("sample needs to be a row of the class table.")
    }
    xraw <- xcmsRaw(object@classes[sample, 1])
    if (type == "raw") {
        title <- paste("Raw data", getRawName(xraw@filepath))
        plotVoidRaw(xraw@env$mz,
                    xraw@env$intensity,
                    xraw@scanindex,
                    xraw@scantime,
                    title = title,
        			legend=legend,...)
    }
    if (type == "peaks") {
        ###Verifying that the peak list is there
        if (nrow(object@peaks) == 0)
            stop("Object needs to be peak picked.")
        pok <- which(object@peaks[, "sample"] == sample)
        if (length(pok) == 0)
            stop("no peaks have been detected for the chose sample.")
        peaklist <- object@peaks[pok, , drop = FALSE]
        largs <- list(...)
        if ("ylim" %in% names(largs)) {
            ylim <- largs[["ylim"]]
            pok <- which(peaklist[, "mzmin"] < ylim[2] &
                            peaklist[, "mzmax"] > ylim[1])
            if (length(pok) > 0) {
                peaklist <- peaklist[pok, , drop = FALSE]
            }
        }
        if (nrow(peaklist) == 0) {
            stop("Nothing to plot.")
        }
        title <- paste("Peaks data", getRawName(xraw@filepath))
        matAllP <- apply(peaklist, 1, function(x, xraw) {
            pok <- which(xraw@env$mz >= x["mzmin"] & xraw@env$mz <= x["mzmax"])
            vscan <- .bincode(pok, breaks = c(xraw@scanindex, length(xraw@env$mz)))
            ppok <- which(vscan >= x["scanMin"] & vscan <= x["scanMax"])
            return(c(pok[ppok]))
        }, xraw = xraw)
        matAllP <- unlist(matAllP)
        vscan <- .bincode(matAllP, breaks = c(xraw@scanindex, length(xraw@env$mz)))
        vmz <- xraw@env$mz[matAllP]
        vint <- xraw@env$intensity[matAllP]
        vindex <- table(vscan)
        tindex <- rep(0, length(xraw@scantime))
        tindex[as.numeric(names(vindex))] <- vindex
        vindex <- tindex
        vindex <- c(0, cumsum(vindex))
        vindex <- vindex[-length(vindex)]
        voo <- order(vscan)
        vmz <- vmz[voo]
        vint <- vint[voo]
        toReturn <- plotVoidRaw(vmz, vint, vindex, xraw@scantime,bandl =  peaklist[,"mz"],title = title,legend=legend,...)
    }
})

keysVariableMetadata<-function(x){
    switch(
        x,
        mzMed="Median mz of the signal in the different files",
        mzMin="Minimum mz of the signal in the different files",
        mzMax="Maximum mz of the signal in the different files",
        sizeSamp="Mean size of the mass traces in the samples",
        scanMin="The minimum scan of the mass trace in the samples",
        scanMax="The maximum  scan of the mass trace in the samples",
        nPeaks="The number of sample in whic the signal have been detected",
        meanSolvent="The mean level of solvent in the sample.",
        signalOverSolventPvalueMean="The mean p-value in the group",
        corSampPeakMean="The mean correlation with the injection peak in the sample",
        signalOverSolventRatioMean="mean of signal on solvent intensity ratio",
        timeShifted="indicator of time scaled signal",
        signalOverSolventPvalueMean="Mean p-value of the grouped signal."
    )
}

setGeneric("exportExpressionSet", function(object, ...)
    standardGeneric("exportExpressionSet"))

#' Export proFIAset to ExpressionSet
#'
#' Eport the data from a proFIAset object as an \code{ExpressionSet} object
#' from the \code{Biobase} package package.
#'
#' @export
#' @param object A proFIAset object.
#' @param colgroup Labels corresponding to the column names
#' of the group table.
#' @return An \code{ExpressionSet} object from the \code{Biobase} package
#' @aliases exportExpressionSet exportExpressionSet,proFIAset-method
#' @examples
#' if(require("plasFIA")&require("Biobase")){
#'     data(plasSet)
#'     eset<-exportExpressionSet(plasSet)
#'     eset
#' }
setMethod("exportExpressionSet", "proFIAset",
          function(object,colgroup=c("mzMed","scanMin","scanMax","nPeaks","corSampPeakMean","signalOverSolventRatioMean","timeShifted","signalOverSolventPvalueMean")){
    if(nrow(object@dataMatrix)==0){
        stop("Data matrix needs  to be created for the object to be converted in ExpressionSet")
    }
    if(!all(colgroup %in% colnames(object@group))) stop("Elements of colgroup are supposed to be the column
                                                        names of the group table of the proFIAset object.")
    fData <- as.data.frame(object@group[,colgroup])
    row.names(fData) <- getUniqueIds(paste("M", sprintf("%0.4f", object@group[, 1]), sep = ""))
    vecLab <- sapply(colgroup,keysVariableMetadata)
    metaData<-data.frame(labelDescription=vecLab)
    eset <- ExpressionSet(object@dataMatrix,featureData=AnnotatedDataFrame(data=fData, varMetadata=metaData))
    eset
})


setGeneric("exportPeakTable", function(object, ...)
    standardGeneric("exportPeakTable"))

#' Export proFIAset as a peak table.
#'
#'Export the data from a proFIAset object as
#'a peak table which containes the values
#'of measured for each variables for each samples
#'and supplementary information.
#'
#' @export
#' @param object A proFIAset object.
#' @param colgroup Labels corresponding to the column names
#' of the group table which will be added to the peak table.
#' @param mval How will missing values be treated, in default they
#' will be set to NA, or you can keep 0.
#' @param  filename The name of the file for the peak table to be exported.
#' @return A dataframe containg the datasets.
#' @aliases exportPeakTable exportPeakTable,proFIAset-method
#' @examples
#' if(require("plasFIA")){
#'     data(plasSet)
#'
#'     #Creating the peak table
#'     ptable<-exportPeakTable(plasSet)
#'     head(ptable)
#'
#'     #Directly in a file
#'     \dontrun{ptable<-exportPeakTable(plasSet,filename="peak_table.tsv")}
#' }
setMethod("exportPeakTable", "proFIAset",
          function(object,colgroup=c("mzMed","corSampPeakMean","meanSolvent","signalOverSolventRatioMean","timeShifted","signalOverSolventPvalueMean"),
                   mval=c("NA","zero"),filename=NULL){
              mval <- match.arg(mval)
              def <- ifelse(mval=="NA",NA,0)
              if(nrow(object@dataMatrix)==0){
                  stop("Data matrix needs  to be created for the object to be converted in peak table.")
              }
              if(!all(colgroup %in% colnames(object@group))) stop("Elements of colgroup are supposed to be the column
                                                                  names of the group table of the proFIAset object.")
              vcnames <- NULL
              vcnames <- c(colnames(object@dataMatrix),colgroup)
              peaktable <- cbind(object@dataMatrix,object@group[,colgroup])
              peaktable <- as.data.frame(peaktable)
              ##Getting the 0
              p0NA <- which(peaktable[,1:ncol(object@dataMatrix)]==0|is.na(peaktable[,1:ncol(object@dataMatrix)]),arr.ind=TRUE)
              if(length(p0NA)!=0){
                  peaktable[,1:ncol(object@dataMatrix)][p0NA]<-def
              }
              colnames(peaktable)<-vcnames
              rownames(peaktable)<-getUniqueIds(rownames(object@dataMatrix))
              if(!is.null(filename)){
                  cpeak <- c("id",colnames(peaktable))
                  ctable <- cbind(rownames(peaktable),peaktable)
                  colnames(ctable)<-cpeak
                  write.table(ctable,sep="\t",row.names = FALSE,quote =TRUE,file = filename)
              }
              return(invisible(peaktable))
          })
