## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=6, fig.path='figures/')

## ----loading, echo=TRUE, warning=FALSE, message=FALSE-------------------------
# loading the packages
library(proFIA)
library(plasFIA)

## ----paths--------------------------------------------------------------------
# finding the directory of the raw files
path <- system.file(package="plasFIA", "mzML")
list.files(path)

## ----profiaset,fig.show="hide",message=FALSE----------------------------------
# defining the ppm parameter adapted to the Orbitrap Fusion
ppm <- 2

# performing the first step of the workflow
plasSet <- proFIAset(path, ppm=ppm, parallel=FALSE)

## ----raw_plot-----------------------------------------------------------------
# loading the spiked molecules data frame
data("plasMols")

# plotting the raw region aroung the Diphenhydramine mass signal
plasMols[7,]
mzrange <- c(plasMols[7,"mass_M+H"]-0.1,plasMols[7,"mass_M+H"]+0.1)
plotRaw(plasSet, type="r", sample=3, ylim=mzrange, size=0.6)

## ----peaks_plot---------------------------------------------------------------
# plotting the filter Dipehnhydramine region.
plotRaw(plasSet, type="p", sample=3, ylim=mzrange, size=0.6)

## ----plot_injection-----------------------------------------------------------
# plotting the injection peak
plotSamplePeaks(plasSet)

## ----group,message=FALSE------------------------------------------------------
# selecting the parameters
ppmgroup <- 1

# due to the experimental design, sample fraction was set to 0.2
fracGroup <- 0.2

# grouping
plasSet <- group.FIA(plasSet, ppmGroup=ppmgroup, fracGroup=fracGroup)

## ----plotEICs-----------------------------------------------------------------
#plotting the EICs of the parameters.		
plotFlowgrams(plasSet,mz=plasMols[4,"mass_M+H"])

## ----find_group---------------------------------------------------------------
# Searching for match group with 2 ppm tolerance
lMatch <- findMzGroup(plasSet,plasMols[,"mass_M+H"],tol=3)

# index of the 40 molecules which may be used with plotEICs
molFound <- data.frame(names=plasMols[,"names"],found=lMatch)
head(molFound)

#Getting the molecules which are not detected
plasMols[which(is.na(lMatch)),]

## ----datamatrix---------------------------------------------------------------
# building the data matrix
plasSet <- makeDataMatrix(plasSet, maxo=FALSE)

## ----impute_fia, warning=FALSE, eval=FALSE------------------------------------
#  # k is supposed to be 3 at minimum, however here we have only 2 sample by class, the results of the imputation are therefore irrelevant.
#  k <- 3
#  
#  #Missing values  imputation using kNN for truncated distribution by default.
#  plasSet <- impute.FIA(plasSet,k=k)
#  
#  #Reinitializing the data matrix.
#  plasSet <- makeDataMatrix(plasSet)
#  
#  #Imputation using random forest.
#  plasSet <- impute.FIA(plasSet,method="randomForest")
#  
#  #As the dataset is ill-suited for missing value imputation we rebuild the data matrix.
#  plasSet <- makeDataMatrix(plasSet)

## ----plot,message=FALSE-------------------------------------------------------
plot(plasSet)

## ----analyzeAcquisitionFIA, eval=FALSE----------------------------------------
#  #selecting the parameters
#  ppm <- 2
#  ppmgroup <- 1
#  fracGroup <- 0.2
#  k <- 3
#  
#  # running the whole workflow in a single step
#  plasSet <- analyzeAcquisitionFIA(path, ppm=ppm, ppmGroup=ppmgroup, k=k,fracGroup = fracGroup,parallel=FALSE)
#  
#  # Running the wholoe workflow in a single step, using parallelism
#  # with the BiocParallel package
#  plasSet <- analyzeAcquisitionFIA(path, ppm=ppm, ppmGroup=ppmgroup, k=k,fracGroup = fracGroup,parallel=TRUE)
#  

## ----export-------------------------------------------------------------------
#Expression Set.
eset <- exportExpressionSet(plasSet)
eset

#Peak Table.
pt <- exportPeakTable(plasSet)

#3 Tables:
dm <- exportDataMatrix(plasSet)
vm <- exportVariableMetadata(plasSet)

## ----multivariate-------------------------------------------------------------
library(ropls)

data("plasSamples")
vconcentration <- plasSamples[,"concentration_ng_ml"]
#vconcentration=(c(100,100,1000,1000,10000,10000)*10^-10)
peakTable <- exportPeakTable(plasSet,mval="zero")

###Cutting the useless column
dataMatrix <- peakTable[,1:nrow(phenoClasses(plasSet))]

## ----plot_summary_opls, echo=FALSE--------------------------------------------
plasSet.opls <- opls(t(dataMatrix),log10(vconcentration),predI = 1,log10L = TRUE, orthoI = NA, devNewL = FALSE,crossvalI=5)

## ----plot_summary_opls_h, eval=FALSE------------------------------------------
#  plasSet.opls <- opls(t(peakTable),scale(log10(vconcentration)),predI = 1,log10L = TRUE, orthoI = NA)

## ----matrix_effect_plot-------------------------------------------------------
matEfInd <- peakTable$corSampPeakMean
nnaVl <- !is.na(matEfInd)
matEfInd <- matEfInd[nnaVl]
ordVi <- order(matEfInd)
matEfInd <- matEfInd[ordVi]
vipVn <- getVipVn(plasSet.opls)[nnaVl]
orthoVipVn <- getVipVn(plasSet.opls, orthoL = TRUE)[nnaVl]
colVc <- rev(rainbow(sum(nnaVl), end = 4/6))
plot(vipVn[ordVi], orthoVipVn[ordVi], pch = 16, col = colVc,
     xlab = "VIP", ylab = "VIP_ortho", main = "VIP_ortho vs VIP.",lwd=3)

##Adding the point corresponding to samples.
points(getVipVn(plasSet.opls)[lMatch],getVipVn(plasSet.opls, orthoL = TRUE)[lMatch], cex=1.2,pch=1,col="black",lwd=2)
legend("topright", legend = c(round(rev(range(matEfInd)), 2),"Spiked molecules."), pch=c(16,16,1),col = c(rev(colVc[c(1, length(colVc))]),1))

## ----plotRaw_exemple_1, message = FALSE,warning=FALSE, results='hide'---------
##Loading the plasFIA dataset
library(plasFIA)
library(proFIA)

data(plasSet)

###Selection of the first sample file
filepath <- phenoClasses(plasSet)[1,1]
filepath

###Loading the raw data
xraw <- xcmsRaw(filepath)

#proFIAset relies on the internal findBandsFIA function to detect m/z bands. The influence of ppm and dmz values can be visualized as follows:
band_list <- findBandsFIA(xraw, ppm = 15, dmz = 0.001)
mzlim <- c(233.067,233.082)
plotRaw(plasSet,sample=2,ylim=mzlim,type="r",legend=FALSE)
abline(h=band_list[,c("mzmin","mzmax")],lwd=0.5,lty=2,col="purple")

## ----plotRaw_exemple_2, message = FALSE,warning=FALSE, results='hide'---------
band_list <- findBandsFIA(xraw, ppm = 2, dmz = 0.0005)
plotRaw(plasSet,sample=2,ylim=mzlim,type="r",legend=FALSE)
abline(h=band_list[,c("mzmin","mzmax")],lwd=0.5,lty=2,col="purple")

## ----group_good_value, message=FALSE, results="hide", eval=FALSE--------------
#  plasSet <- group.FIA(plasSet,ppmGroup=5,dmzGroup=0.001,fracGroup=3/18,sleep=0.001)

## ----group_wrong_value, message=FALSE, results="hide",eval=FALSE--------------
#  plasSet <- group.FIA(plasSet,ppmGroup=1,dmz=0.001,fracGroup=3/18,sleep=0.001)

## ----missing_values_1---------------------------------------------------------
data(plasSet)

###You can reset the data matrix this way
plasSet <- makeDataMatrix(plasSet)

###Before imputation.
plot(plasSet)

## ----missing_values_2---------------------------------------------------------
plasSet <- impute.randomForest(plasSet)

###After imputation.
plot(plasSet)

## ----cheat_sheet--------------------------------------------------------------
system.file(package="proFIA")

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

