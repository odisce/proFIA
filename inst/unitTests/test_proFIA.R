test_findPeaksLimits <- function() {
    xsq=seq(-3,3,length=100)
    ysq=dnorm(xsq,mean=0,sd=0.5)
    checkEquals(proFIA:::findPeaksLimits(ysq,49,51),c(2,99))
    ysq=ysq+dnorm(xsq,mean=1.5,sd=0.5)
    checkEquals(proFIA:::findPeaksLimits(ysq,49,51),c(2,63))
}


test_segmentCurve <- function(){
    x=1:150
    y=c(rep(10,length=15),seq(10,100,length=30),seq(100,40,length=40),seq(40,10,length=65))
    checkEquals(proFIA:::segmentCurve(x,y,10000),c(1,150))
    checkEquals(proFIA:::segmentCurve(x,y,5),c(1,16,45,85,150))
}


test_binaryClosest <- function(){
    xsq=c(120.692,121.736,184.464,198.658,223.813,294.816,306.351,350.519,359.264,375.986)
    checkEquals(proFIA:::binaryClosest(xsq,350,0,length(xsq)),8)
}

test_smooth.BWF <-function(){
    xsq=seq(0,40,length=100)
    ysq=dnorm(xsq,mean=15,sd=5)
    nsq=c(0.002,-0.003,0,-0.003,-0.005,0,-0.005,0,0.001,-0.001,-0.001,-0.004,-0.006,0.004,0,
          0.004,-0.007,0,0.003,-0.001,-0.001,-0.004,0.01,0.003,0.009,0,-0.004,-0.008,
          -0.004,0.001,0.005,0.001,-0.006,0.003,-0.009,-0.002,-0.002,-0.01,0.004,-0.004,
          -0.002,0.006,0.003,-0.001,0.002,0.001,0.008,0.001,0.002,0.001,-0.002,0,0.008,
          -0.005,0,0.006,-0.002,-0.006,-0.005,-0.002,-0.005,-0.002,-0.003,-0.001,0.002,
          0.003,0.008,-0.001,-0.002,0.007,0.008,0,-0.002,-0.01,0.007,-0.005,-0.006,0.015,
          0.001,-0.001,-0.014,-0.001,0.004,0.002,0.004,0.009,0.005,-0.003,0.003,-0.001,
          -0.003,-0.015,-0.004,0.002,0.01,0.002,0.003,-0.005,0,-0.011)
    ysq=ysq+nsq
    checkEqualsNumeric(sum(proFIA:::smooth.BWF(ysq,freq = 0.2)),2.446505,tol=10^-5)
}

# test_findFIASignal<-function(){
#     if(require(plasFIA)){
#     data(plasSet)
#     xraw <- xcms:::xcmsRaw(plasSet@classes[2,1])
#     tp <- invisible(findFIASignal(xraw,ppm=2,es=plasSet@noiseEstimation,pvalthresh=0.01))
#     checkEquals(tp$injscan,17)
#     ##Differnece of minpack.lm on Linux.
#     checkEqualsNumeric(sum(tp$injpeak),36.57748,tolerance=0.2)
#     checkEqualsNumeric(nrow(tp$matrix),604,tolerance=0.05)
#     }
# }

test_group.FIA<-function(){
    if(require(plasFIA)){
    data(plasSet)
    tg<-proFIA:::group.FIA(plasSet,ppmGroup=0.5,fracGroup=1,dmz=0.0005)
    checkEquals(nrow(tg@group),219)
    }
}