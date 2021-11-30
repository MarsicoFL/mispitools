#' Decision Threshold simulation: a function for computing likelihood ratio decision threshold.
#'
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#' @param weight The differential weight between false positives and false negatives. A value of 10 is suggested. 
#'
#' @export
#'
#' @examples
#' #library(forrel)
#' #x = linearPed(2)
#' #x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' #x = profileSim(x, N = 1, ids = 2)[[1]]
#' #datasim = makeLRsims(x, missing = 5, 1000)
#' #DTsim(datasim, 10)



DTsim = function(datasim, weight) {

nsims = nrow(datasim)
TPED = datasim$Related # We define the obtained values for H1
RPED = datasim$Unrelated # We define the obtained values for H2

ValoresLR = seq(1, nsims, length.out=nsims) 
FPs = 0
FNs = 0
for(i in 1:10000) { FPs[i] = sum(RPED > ValoresLR[i]); #False positives.
                    FNs[i] = sum(TPED < ValoresLR[i])} #False negatives.

#The differential weight between false positives and PMI.
Dis = sqrt(((FNs/nsims))^2+((weight*FPs/nsims)^2)) #Weighted Euclidean Distance
Tabla = base::data.frame(x=ValoresLR, y= Dis)
DT = which.min(Tabla$y)
print(paste("Decision threshold is:", DT))
}
