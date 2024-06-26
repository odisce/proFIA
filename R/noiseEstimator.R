########
#An object created to estimate the noise of an ensemble on CFIA acquisition.
######

#' An S4 class to represent heteroscedasctic noise of MS.
#'
#' @slot bins The limit of the bins on which the noise have been estimated.
#' @slot varmean The estimate of the mean of the variance in each bin.
#' @slot size The size of each bins in number of elements.
#' @slot estimation The estimation function, if a model have been fitted.
#' @slot filelist The list of the fitted file.
#' @slot estimated A boolean indicating if a model have been fitted.
#' @slot frac A numeric giving the fraction of the point ot be conserved when estimating the noise model.
#' @slot intlim The maximum and minmum intensity on which the estimationis done.
#' @slot reglim The limits of the regression of the model on these data.
#' detected for each experiment.
setClass(
    "noiseEstimation",
    slot = list(
        bins = "numeric",
        varmean = "numeric",
        size = "integer",
        estimation = "function",
        filelist = "character",
        estimated = "logical",
        frac = "numeric",
        intlim = "numeric",
        reglim = "numeric"
    ),
    prototype = list(
        bins = numeric(0),
        varmean = numeric(0),
        size = integer(0),
        estimation = new("function"),
        filelist = character(0),
        estimated = FALSE,
        frac = 0.99,
        intlim = numeric(0),
        reglim = numeric(0)
    )
)


getMidBins <- function(bins) {
    return((bins[-length(bins)] + bins[-1]) / 2)
}

###Determine the first point on which to cut.
cuttingPoint <- function(osize, frac = 0.995) {
    osize <- osize[-1]
    vmax <- sum(osize)
    nsize <- cumsum(osize)
    vt <- frac * vmax
    poscut <- (which(nsize > vt))[1]
    poscut
}

interToRegress <- function(object,
                           quant = 0.95,
                           thresh = NULL) {
    vNoise <- NULL
    if (class(object) != "noiseEstimation") {
        vNoise <- object@noiseEstimation
    } else{
        vNoise <- object
    }
    trx <- getMidBins(vNoise@bins)
    try <- vNoise@varmean
    vsize <- vNoise@size
    trw <- log10(vsize)
    trw <- log10(vsize)
    vsize <- vsize / sum(vsize)
    vsize <- cumsum(vsize)
    plsup <- which(vsize > quant)[1] + 1
    
    if(plsup<=2){
      return(c(NA,NA))
    }
    plmin <- NULL
    ###Implement th enon loess method.
    if (is.null(thresh)) {
        mag <- movavg(try[1:plsup], max(min(7,floor(plsup/2)),2), type = "s")
        mag <- mag[3:(length(mag) - 2)]
        plmin <- which.min(mag) + 2
        
    } else{
        plmin <- which(try >= thresh)[1]
    }
    return(c(plmin, plsup))
}


setMethod("show", "noiseEstimation", function(object) {
    cat(
        "an \"noiseEstimator\" object with ",
        length(object@bins) - 1,
        " bins with mass range,",
        range(object@bins),
        ".\n"
    )
    
    if (is.null(object@estimation()))
        cat("No estimation have been done.\n")
    else
        cat("An estimation have been done.\n")
    memsize <- object.size(object)
    cat("Memory usage:", signif(memsize / 2 ^ 20, 3), "MB\n")
})



#' Plot the estimated noise from a proFIAset object.
#'
#' Plot an intensity vs variances plot for the noise estimated from a
#' MS acquisition. If a model is fitted to the data it will be plotted.
#' Only the interval used for the regression and on which the estimation will be used is shown.
#'
#' @export
#' @param object A noise estimation object or a \code{\link{proFIAset-class}} object.
#' @param xlim The xlim paramter to be passed to plot.
#' @param ylim The ylim paramter to be passed to plot.
#' @param ... SUpplementary arguments to be passed to plot
#' @param plotNoise plotNoise,proFIAset-method
#' @return  The plotted value as an x y list.
#' @examples
#' if(require(plasFIA)){
#'   data(plasSet)
#'   plotNoise(plasSet)
#' }
plotNoise <- function(object,
                      xlim = NULL,
                      ylim = NULL,
                      ...) {
    if (class(object) == "proFIAset") {
        object <- object@noiseEstimation
    } else if (class(object) != "noiseEstimation") {
        stop("plotNoise require a proFIAset or a noiseEstimation object ot be given as input.")
    }
    if (length(object@bins) == 0 |
        length(object@varmean) == 0 |
        length(object@size) == 0)
        stop("Nothing to plot, empty noise estimation object.\n")
    main <- "Noise variance estimation"
    toplot <- NULL
    dtp <- NULL
    
    ##Determining the position of cutting.
    
    vBins <- getMidBins(object@bins)
    if (object@estimated)
        main <- paste(main, "and fitted model")
    
    seqx <- NULL
    
    ###Plotting the points
    seqreg <- NULL
    if(object@estimated){
        seqreg <- object@reglim[1]:object@reglim[2]
    }else{
        seqreg <- 1:length(object@size)
    }
    xsp <- vBins[seqreg]
    ysp <- object@varmean[seqreg]
    colv <- log10(object@size[seqreg])
    #cat(colv)
    pos0 <- which(object@size[seqreg] == 0)
    if (length(pos0) != 0) {
        xsp <- xsp[-pos0]
        ysp <- ysp[-pos0]
        colv <- colv[-pos0]
    }
    
    ###Generating the gam of colors to be used.
    colv <- colorRampPalette(c("blue", "cyan"))(50)[as.numeric(cut(colv, breaks = 50))]
    colv <- c(colv[1], colv)
    
    ###Generating the oclor gradient
    if (is.null(xlim))
        xlim <- range(xsp)
    if (is.null(ylim))
        ylim <- range(ysp)
    
    plot(
        xsp,
        ysp,
        cex = 0.7,
        col = colv,
        pch = 19,
        xlim = xlim,
        ylim = ylim,
        xlab = "Intensity bins",
        ylab =
            "Estimated variance",
        main = main,
        ...
    )
    
    vtexleg <- c("Estimated variance(most numerous)",
                "Estimated variance(less numerous)")
    vcolleg <- colv[c(1, length(colv))]
    vltyleg <- c(FALSE, FALSE)
    vpchleg <- c(19, 19)
    
    
    if (object@estimated) {
        if (is.null(xlim)) {
            seqx <- seq(1, object@bins[length(object@bins)], length = 1000)
        } else{
            seqx <- seq(xlim[1], xlim[2], length = 1000)
            seqy <- sapply(seqx, object@estimation)
            lines(seqx, seqy, col = "red")
            vtexleg <- c(vtexleg, "Fitted model")
            vcolleg <- c(vcolleg, "red")
            vltyleg <- c(vltyleg, 1, 2)
            vpchleg <- c(vpchleg, NA, NA)
        }
    }
    legend("topleft",
           vtexleg,
           col = vcolleg,
           lty = vltyleg,
           pch = vpchleg)
    return(invisible(list(x = xsp, y = ysp)))
}



###Fit the modle using nls
makeBins <- function(n, minInt, maxInt) {
    ###Bins.
    res <- seq(log(minInt), log(maxInt), length = n)
    exp(res)
}




### The weigths used for the noise points correspond ot the standard deviation on the variance estimator.


##Maybe exported later.
setGeneric("fitModel", function(object, ...)
    standardGeneric("fitModel"))
setMethod("fitModel", "noiseEstimation", function(object,
                                                  modelType =
                                                      c("loess", "polynom", "plinear"),
                                                  weighted = FALSE,
                                                  threshold = NULL,
                                                  quant = 0.95, absThreshold = NULL) {
    modelType <- match.arg(modelType)
    xvalues <- getMidBins(object@bins)
    yvalues <- object@varmean
    vsize <- object@size
    wvalues <- sqrt(object@size)
    #Determining a limiting frequency to fit the data.
    ##The zero are removed, as the first points which is heavely biaised because of the low intensity of the filter.
    p0 <- which(object@size == 0)
    if (length(p0) != 0) {
        xvalues <- xvalues[-c(1, p0)]
        yvalues <- yvalues[-c(1, p0)]
        wvalues <- wvalues[-c(1, p0)]
        vsize <- vsize[-c(1, p0)]
    }
    ###Determining the injection model.
    beginning <- object@varmean[2]
    yvalues <- yvalues - beginning
    
    
    ###Cutting the vector on which do th regression
    limreg <- interToRegress(object, quant, threshold)
    if(is.na(limreg[1])) return(NA)
    
    if(!is.null(absThreshold)){
    	pLim <- which(xvalues>absThreshold)
    	if(length(pLim)==0){
    		stop("Absolute noise threshold is too high, no point considered for regression.")
    	}
    	pLim <- pLim[1]
    	limreg[1] <- max(limreg[1],pLim[1])
    }
    
    
    iilim <-  xvalues[c(limreg[1], limreg[2])]
    #print(limreg)
    
    xvalues <- xvalues[limreg[1]:limreg[2]]
    yvalues <- yvalues[limreg[1]:limreg[2]]
    wvalues <- wvalues[limreg[1]:limreg[2]]
    ###Checking is there is size 0 bins to remove.
    
    
    if (!weighted) {
        wvalues <- rep(1, length(wvalues))
    } else{
        wvalues <- wvalues / log2(yvalues)
    }
    ###The value of the interept is determined as the
    if (modelType == "polynom") {
        estimationModel <- lm(yvalues ~ xvalues + I(xvalues ^ 2) - 1, weights = wvalues)
        resCoef <- coef(estimationModel)
        object@estimation <- function(x) {
            if (x > iilim[2]) {
                return(0)
            }
            if (x < iilim[1]) {
                return(NA)
            }
            return(resCoef[2] * x ^ 2 + resCoef[1] * x + beginning)
        }
    } else if (modelType == "plinear") {
        object@estimation <- function(x) {
        	bmax <- ifelse(x>object@bins[limreg[2]],TRUE,FALSE)
        	bmin <- ifelse(x<object@bins[limreg[1]],TRUE,FALSE)
        	pok <- 1:length(x)
        	if(any(bmin)|any(bmax)){
        		pok <- pok[-c(bmin,bmax)]
        	}
        	
        	res <- numeric(length(x))
        	
        	
        	if(any(bmin)){
        		res[bmin] <- object@bins[limreg[1]]
        	}
        	if(any(bmax)){
        		res[bmax] <- object@bins[limreg[2]]
        	}
            cbin = .bincode(x[pok], object@bins)
            posSeg <- (x[pok] - object@bins[cbin]) / (object@bins[cbin + 1] - object@bins[cbin])
            res[pok] <- object@varmean[cbin]+(object@varmean[cbin+1]-object@varmean[cbin]) * posSeg
            return(res+beginning)
        }
    } else if (modelType == "loess") {
        m.lo <- loess(yvalues ~ xvalues, weights = wvalues)
        bins <- getMidBins(object@bins)
        object@estimation <- function(x,nullVal=yvalues[limreg[1]]) {
        	vreturn = predict(m.lo, x)
        	psup <- which(x > bins[object@reglim[2]])
            if (length(psup)>0) {
            	vreturn[psup] <- object@varmean[object@reglim[2]]
            }
        	pmin <- which(x < object@reglim[1])
            if (length(pmin)>0) {
            	vreturn[pmin] <- nullVal
            }
            return(vreturn+beginning)
        }
    }
    object@estimated <- TRUE
    object@reglim <- limreg
    object
})


estimateNoiseFile <-
    function(fname,
             ppm,
             nbin = 120,
             denoising = c("BWF", "SavGol"),
             f = c("TIC","regression"),
             graphical =
                 FALSE,
             maxInt = NULL,
             minInt = NULL,
    		 includeLowest=TRUE,
    		 dmz=0.0005,
             ...) {
    	if(length(fname)>1){
    		stop("To estimate the noise on more than one file use estimateNoiseMS.")
    	}
    	xraw <- NULL
    	if(class(fname)=="character"){
    		xraw <- xcmsRaw(fname)
    	}else{
    		xraw <- fname
    	}
        f <- match.arg(f)
        if(f=="regression"){
          sizepeak <-
            determiningInjectionZone(xraw,...)
        }else if(f=="TIC"){
          sizepeak <- determiningInjectionZoneFromTIC(xraw,...)
        }
        ###Checking if the sizepeak algorithm is bad.
        
        beginningInjection <- sizepeak[1]
        tabROI <- findBandsFIA(
            xraw,
            firstScan = 1,
            lastScan = length(xraw@scantime),
            ppm = ppm,
            sizeMin = (sizepeak[3] - sizepeak[1])*0.5,
            dmz = dmz,
            beginning = sizepeak[1],
            end=sizepeak[2],
            fracMin = 0.25
        )
        denoising <- match.arg(denoising)
        funcName <- paste("smooth.", denoising, sep = "")
        
        n <- length(xraw@scantime)
        ##Make the binning interval
        
        if (is.null(maxInt))
            maxInt <- max(xraw@env$intensity)
        if (is.null(minInt))
            minInt <- max(xraw@env$intensity)
        #cat("minInt ",minInt,"maxInt",maxInt,"\n")
        binlim <- makeBins(nbin + 1, minInt, maxInt)
        if(includeLowest){
        binlim[1] <- 0
        }
        binlim[length(binlim)] <- max(maxInt, binlim[length(binlim)]) + 1
        matres <-
            matrix(0,
                   nrow = 2 * length(xraw@scantime),
                   ncol = nrow(tabROI))
        if (graphical) {
            for (i in 1:nrow(tabROI)) {
                v <- tabROI[i,]
                rEIC <- rawEIC(xraw, mzrange = c(tabROI[i, "mzmin"], tabROI[i, "mzmax"]))
                ntitle <- paste(
                    "raw and smoothed EIC ",
                    sprintf("%0.2f", tabROI[i, "mzmin"]),
                    "-",
                    sprintf("%0.2f", tabROI[i, "mzmax"]),
                    sep =
                        ""
                )
                plot(xraw@scantime, rEIC$intensity, type = "l")
                fv <- putZero(rEIC$intensity, do.call(funcName, list(x = rEIC$intensity)))
                lines(xraw@scantime,
                      fv[1:n],
                      type = "l",
                      col = "red")
                noiseestimate <- rEIC$intensity - fv[1:n]
                bins <- .bincode(fv[1:n], binlim)
                matres[, i] <- c(noiseestimate, bins)
            }
        } else{
            matres <- apply(tabROI, 1, function(v, xr, binl) {
                rEIC <- rawEIC(xr, mzrange = c(v["mzmin"], v["mzmax"]))
                fv <- putZero(rEIC$intensity, do.call(funcName, list(x = rEIC$intensity)))
                noiseestimate <- rEIC$intensity - fv[1:n]
                bins = .bincode(fv[1:n], binl)
                c(noiseestimate, bins)
            }, xr = xraw, binl = binlim)
        }
        ###Obtaining the noise estimate
        matnoise <- as.numeric(matres[1:n,])
        matbin <- as.numeric(matres[(n + 1):(2 * n),])
        tbin <- table(matbin)
        vecmean <- aggregate(matnoise, by = list(bin = matbin), function(x) {
            mean(x ^ 2)
        })    ###Variance
        #vecvar=aggregate(matnoise,by=list(bin=matbin),function(x){return(sum((x^2-mean(x^2))^2)/(length(x)*(length(x)-1)))})  ###Standard deviation of the mean
        #print(vecvar)
        #vecvar[,2]=sapply(vecvar[,2],function(x){if(length(x)>1){return(x)}else if(is.nan(x)|is.na(x)){return(0)}else{return(x)}})
        
        vnum <- rep(0, nbin)
        vmean <- rep(0, nbin)

        vmean[vecmean[, 1]] <- vecmean[, 2]
        vmean <- sapply(vmean, function(x) {
            if (is.nan(x) | is.na(x)) {
                return(0)
            } else{
                return(x)
            }
        })
        vnum[as.numeric(names(tbin))] = tbin
        return(list(bins = binlim, noisevar = sqrt(vmean), size = as.integer(vnum)))
        #return(list(bins=binlim,noisevar=sqrt(vmean),size=as.integer(vnum),vecsd=vvar/sapply((2*sqrt(vmean)),function(x){max(x,1)})))
    }


#' Estimate the noise of a mass spectrometer using multiple MS acquisition.
#'
#' Determine the variances of the noise in function of the intensity
#' from multiples FIA acquisitions, using the method from Wentzell and Tarazuk(2014)
#' \emph{Characterization of heteroscedastic measurement noise in the absence of replicates}
#'
#' @export
#' @param list_files A list of files in which the noise should be estimated.
#' @param ppm The authorized deviation between scans in ppm, this parameter
#' will also be used to fuse the bands if there are close enough.
#' @param dmz The dmz parameters to be passed to findBandsFIA, see \link{findBandsFIA}
#' @param nBin The number of intensity bins to be used.
#' @param minInt The  minimum intensity expected in all the files.
#' @param maxInt The  maximum intensity expected in all the files.
#' @param f Which method should be used to find the injection peaks. Should be 'TIC' or'regression'
#' @param includeZero Should the left bin start a 0.
#' @param parallel Shall parrallelism be used.
#' @param BPPARAM A BiocParallelParam object to be used for parallelism
#' if parallel is TRUE.
#' @return A  noise estimator object.
#' @examples
#' ##Listing the files in plasFIA
#' if(require(plasFIA)){
#'   list_mzml <- list.files(system.file(package="plasFIA","mzML"),full.names=TRUE)
#' 
#'   ##For speed purpose
#'   list_mzml <- list_mzml[1:2]
#'   es <- estimateNoiseMS(list_mzml,2,parallel=FALSE)
#' }
estimateNoiseMS <-
    function(list_files,
             ppm,
             nBin = 500,
             minInt = 500,
             maxInt = 10 ^ 8,
             f=c("TIC","regression"),
             parallel,
    		 includeZero=TRUE,
    		 dmz=0.0005,
             BPPARAM = NULL) {
        f <- match.arg(f)
        res = NULL
        if (parallel & requireNamespace("BiocParallel")) {
            if (is.null(BPPARAM))
                BPPARAM <- bpparam()
            res <- bplapply(
                as.list(list_files),
                estimateNoiseFile,
                nbin = nBin,
                ppm = ppm,
                maxInt = maxInt,
                minInt =
                    minInt,
                includeLowest=includeZero,
                dmz=dmz,
                f=f,
                BPPARAM = BPPARAM
            )
            
            
        } else{
            if (parallel) {
                warning(
                    "BioCParallel package is not installed, impossible to use parallelism. Single core is used."
                )
            }else{
                res <- sapply(
                    list_files,
                    estimateNoiseFile,
                    nbin = nBin,
                    ppm = ppm,
                    maxInt = maxInt,
                    minInt =
                        minInt,
                    f=f,
                	includeLowest=includeZero,
                    simplify = FALSE
                )  
            }
        }
        
        
        if (is.null(minInt) |
            is.null(maxInt)) {
            stop("No intensities provided noise cannot be estimated")
        }

        vInt <- rep(0, length(res[[1]]$noise))
        vSize <- as.integer(rep(0, length(res[[1]]$noise)))
        for (i in 1:length(res)) {
            vInt <- (vSize * (vInt) + res[[i]]$size * res[[i]]$noise) / sapply((res[[i]]$size +
                                                                                   vSize), function(x) {
                                                                                       max(x, 1)
                                                                                   })
            vSize <- vSize + res[[i]]$size
        }
        if(!is.character(list_files)){
        	list_files <- sapply(list_files,function(x){
        	as.character(x@filepath)
        		})
        }
        
        est <- new("noiseEstimation")
        est@intlim <- c(minInt, maxInt)
        est@varmean <- vInt
        est@size <- as.integer(vSize)
        est@bins <- res[[1]]$bins
        est@filelist <- list_files
        est@estimated <- FALSE
        est
    }


#H0 test if the seq int seq and the sequence hseq may be a putati
calcPvalue <- function(nes, intseq, hseq) {
    if (length(intseq) != length(hseq))
        stop("length of the 2 seqs is different")
    if (!nes@estimated)
        stop("A model have not been fitted to the data")
    dist <- (intseq - hseq)
    #cat(dist,"\n")
    hvar <- sapply(hseq, nes@estimation, nullVal = NA)
    
    ##Points which haven't any values are not estimated.
    invpos <- which(sapply(hvar, is.na) | hvar < 0)
    if (length(invpos) != 0) {
        hvar <- hvar[-invpos]
    }
    if (length(hvar) == 0 &
        any(intseq > 2 * nes@bins[nes@reglim[1]]))
        return(0)
    
    if (length(hvar) == 0)
        return(1)
    
    testStat <- sum(dist)
    pvalue <- 2 * (1 - pnorm(testStat, 0, sum(hvar)))
    pvalue
}
