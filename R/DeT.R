#' Decision Threshold: a function for computing likelihood ratio decision threshold.
#'
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#' @param weight The differential weight between false positives and false negatives. A value of 10 is suggested. 
#'
#' @export
#' @return A value of Likelihood ratio suggested as threshold based on false positive-false negative trade-off.
#' @examples
#' library(forrel)
#' x = linearPed(2)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' datasim = simLRgen(x, missing = 5, 10, 123)
#' DeT(datasim, 10)



DeT = function(datasim, weight) {

nsims = nrow(datasim)
TPED = datasim$Related 
RPED = datasim$Unrelated 

ValoresLR = seq(1, nsims, length.out=nsims) 
FPs = 0
FNs = 0
for(i in 1:10000) { FPs[i] = sum(RPED > ValoresLR[i]);
                    FNs[i] = sum(TPED < ValoresLR[i])}

Dis = sqrt(((FNs/nsims))^2+((weight*FPs/nsims)^2)) 
Tabla = base::data.frame(x=ValoresLR, y= Dis)
DT = which.min(Tabla$y)
print(paste("Decision threshold is:", DT))
}
