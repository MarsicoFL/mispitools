#' Threshold rates: a function for computing error rates and Matthews correlation coefficient of a specific LR threshold.
#'
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#' @param threshold Likelihood ratio threshold selected for error rates calculation.
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#'
#' @export
#' @return Values of false positive and false negative rates and MCC for a specific LR threshold.
#' @examples
#' library(forrel)
#' x = linearPed(2)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' datasim = simLRgen(x, missing = 5, 10, 123)
#' Trates(datasim, 10)



Trates = function(datasim, threshold) {

nsims = nrow(datasim)
TPED = datasim$Related 
RPED = datasim$Unrelated 

FPR = sum(RPED > threshold)/nsims 
FNR = sum(TPED < threshold)/nsims 
TPR = sum(RPED < threshold)/nsims 
TNR = sum(TPED > threshold)/nsims 
MCC = (TPR*TNR-FPR*FNR)/(sqrt(TPR+FPR)*sqrt(TPR+FNR)*sqrt(TNR+FPR)*sqrt(TNR+FNR))

print(paste("FNR =", FNR, ";  FPR =", FPR,";  MCC =", MCC ))

}
