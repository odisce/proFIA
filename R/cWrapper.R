#' Detect band in a FIA acquisition
#'
#' Detect bands of points with similar mass in conscutive scans.
#' Points may be moved if a better candidates is found.
#'
#' @export
#' @param xraw An xcmsRaw object as returned by \code{\link[xcms]{xcmsRaw}}.
#' @param firstScan The first scan to be considered, 1 for general use.
#' @param lastScan  The last scan to be considered.
#' @param ppm The mass deviation in ppm for point in consecutives scans.
#' @param sizeMin The minimum size of a band.
#' @param dmz The minimum mass tolerance,useful for small masses
#' @param beginning The scan of the injection. May be determined using 
#' \code{\link{determiningInjectionZone}}.
#' @param end The end of injection, may be estimated using \code{\link{determiningInjectionZone}}.
#' @param nIso The minimum number of consecutive points for a signal to be detected.
#' @param fracMin The minimum fraction of points necessary in beginning:end for a signal to be detected.
#' contaminated by solvent.
#' @return A vector contaning the inject peak
#' @aliases findBandsFIA
#' @examples
#' #Getting the path of a file.
#' if(require(plasFIA)){
#'   path_raw<-list.files(system.file(package="plasFIA","mzML"),full.names=TRUE)[2]
#' 
#'   #Opening the file with xcms
#'   xraw<-xcmsRaw(path_raw)
#' 
#'   #Getting the injection scan
#'   gp<-determiningInjectionZone(xraw)
#' 
#'   #performing band detection.
#'   tbands<-findBandsFIA(xraw,ppm = 2,sizeMin = gp[3]-gp[1],beginning=gp[1],end=gp[2])
#'   head(tbands)
#' }
findBandsFIA <-
    function(xraw,
             firstScan = 1,
             lastScan = length(xraw@scantime),
             ppm = 2,
             sizeMin = NULL,
             dmz = 0.0005,
             beginning=NULL,
    		 end=NULL,
             nIso = 3,
    		 fracMin = 0.6) {
    	if(is.null(beginning)|is.null(end)|is.null(sizeMin)){
    		vp <- determiningInjectionZone(xraw)
    		beginning <- ifelse(is.null(beginning),vp[1],beginning)
    		end <- ifelse(is.null(end),vp[2],end)
    		sizeMin <- ifelse(is.null(sizeMin),3,sizeMin)
    	}

        mz <- as.numeric(xraw@env$mz)
        int <- as.numeric(xraw@env$intensity)
        sc_ind <- as.integer(xraw@scanindex)
        fS <- as.integer(firstScan)
        lS <- as.integer(lastScan)
        mS <- as.integer(length(xraw@scanindex))
        length_mz <- as.integer(length(mz))
        ppm_dev <- as.numeric(ppm)
        sM <- as.integer(sizeMin)
        beginning <- as.integer(beginning)
        end <- as.integer(end)
        nIso <- as.integer(nIso)
        Dmz <- as.numeric(dmz)
        Fmin <- as.numeric(fracMin)
        #cat("ppm :",ppm_dev,"dmz :",Dmz,"nmz",length(mz),"flmscans ",fS,lS,mS,"\n")
        message(paste("Beginning band detection for ",getRawName(xraw@filepath),"\n",sep=""))
        to_return <- .Call(
            "findBandsFIACentroids",
            mz,
            int,
            sc_ind,
            as.numeric(xraw@scantime),
            fS,
            lS,
            mS,
            length_mz,
            ppm_dev,
            nIso,
            beginning,
            end,
            sM,
            Dmz,
            Fmin,
            PACKAGE = "proFIA"
        )
        colnames(to_return) <-
            c(
                "mzmin",
                "mzmax",
                "mz",
                "scan min",
                "scan max",
                "size",
                "meanSolvent",
                "maxIntensity"
            )
        to_return
    }

#dyn.load("C:/Users/AD244905/Documents/fia-tandem/C/peaksGroup.dll")
findLimDensity <- function(dens, istart = 0, state) {
    if (!is.double(dens))
        dens <- as.double(dens)
    unlist(
        .C(
            "findLimDensity",
            dens,
            as.integer(length(dens)),
            as.integer(istart - 1),
            linflex = integer(1),
            rinflex = integer(1),
            state = as.integer(state),
            PACKAGE = "proFIA"
        )[4:6]
    ) + c(1, 1, 0)
}




#XCMS furnished function.
findEqualGreaterM <- function(x, values) {
    if (!is.double(x))
        x <- as.double(x)
    if (!is.double(values))
        values <- as.double(values)
    .C(
        "FindEqualGreaterM",
        x,
        length(x),
        values,
        length(values),
        index = integer(length(values)),
        PACKAGE = "proFIA"
    )$index + 1
}

findLinesAboveNoise <- function(x, noise, numPoints, first) {
    if (!is.double(x))
        x <- as.double(x)
    if (!is.double(noise))
        noise <- as.double(noise)
    if (!is.integer(numPoints))
        numPoints <- as.integer(numPoints)
    if (!is.integer(first))
        first <- as.integer(first)
    unlist(
        .C(
            "lineAboveNoise",
            first,
            as.integer(0),
            as.integer(0),
            x,
            as.integer(length(x)),
            numPoints,
            noise
            ,
            PACKAGE = "proFIA"
        )
    )[1:3]
}



findPeaksLimits <- function(sEIC, posL, posR) {
    if (!is.double(sEIC))
        sEIC <- as.double(sEIC)
    ##shift of 1 in C
    if(posR>=length(sEIC)) posR=length(sEIC)-1
    if(posL<=1) posL=2
    if (!is.integer(posL))
        posL <- as.integer(posL - 1)
    if (!is.integer(posR))
        posR <- as.integer(posR - 1)
    resL = .C(
        "findLim",
        sEIC,
        posL,
        posR,
        as.integer(length(sEIC)),
        integer(1),
        integer(1),
        PACKAGE = "proFIA"
    )
    c(resL[[5]], resL[[6]]) + 1
}


segmentCurve <- function(time, int, eps) {
    if (!is.double(time))
        time <- as.double(time)
    if (!is.double(int))
        int <- as.double(int)
    if (!is.double(eps))
       eps <- as.double(eps)
    segm = .Call(
        "segmentCurveW",
        time,
        int,
        eps,
        as.integer(length(time)),
        PACKAGE = "proFIA"
    )+1
    return(segm)
}



binaryClosest <- function(oseq, value, imin, imax) {
    if (!is.double(oseq))
        oseq <- as.double(oseq)
    if (!is.double(value))
        value <- as.double(value)
    if (imin < 1) {
        imin <- as.integer(1)
    } else{
        imin <- as.integer(imin)
    }
    if (imax > length(oseq)) {
        imax <- as.integer(length(oseq))
    } else{
        imax <- as.integer(imax)
    }
    .C("binarySearch",
       oseq,
       value,
       imin,
       imax,
       imid = integer(1),
       PACKAGE = "proFIA")$imid + 1
}

putZero <- function(oseq, filtseq) {
    if (!is.double(oseq))
        oseq <- as.double(oseq)
    if (!is.double(filtseq))
        filtseq <- as.double(filtseq)
    .C(
        "replaceByZero",
        rseq = filtseq,
        oseq,
        as.integer(length(oseq)),
        PACKAGE = "proFIA"
    )$rseq + 1
}

linearInterpolation <- function(oseq, imin, imax) {
    if (!is.double(oseq))
        oseq <- as.double(oseq)
    if (!is.integer(imin))
        imin <- as.integer(imin)
    if (!is.integer(imax))
        imax <- as.integer(imax)
    res <- .C(
        "linearInterpolation",
        rseq = oseq,
        imin,
        imax,
        ciso = integer(1),
        PACKAGE = "proFIA"
    )
    return(c(res$ciso, res$rseq))
}


checkIso <- function(intensity,
                     niso = 3,
                     smax = NULL) {
    if (is.null(smax))
        smax <- as.integer(length(intensity))
    if (!is.integer(smax))
        smax <- as.integer(length(smax))
    if (!is.double(intensity))
        intensity <- as.double(intensity)
    if (!is.integer(niso))
        niso <- as.integer(niso)
    .C("checkIso",
       intensity,
       niso,
       valid = as.integer(1),
       smax,
       PACKAGE = "proFIA")$valid
}