#' Threshold rates: a function for computing error rates and Matthews correlation coefficient of a specific LR threshold.
#'
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#' @param threshold Likelihood ratio threshold selected for error rates calculation.
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#'
#' @export
#'
#' @examples
#' #library(forrel)
#' #x = linearPed(2)
#' #x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' #x = profileSim(x, N = 1, ids = 2)[[1]]
#' #datasim = makeLRsims(x, missing = 5, 1000)
#' #Trates(datasim, 10)



Trates = function(datasim, threshold) {

nsims = nrow(datasim)
TPED = datasim$Related # We define the obtained values for H1
RPED = datasim$Unrelated # We define the obtained values for H2

FPR = sum(RPED > threshold)/nsims #False positives.
FNR = sum(TPED < threshold)/nsims #False negatives.
TPR = sum(RPED < threshold)/nsims #True positives.
TNR = sum(TPED > threshold)/nsims #True negatives.
MCC = (TPR*TNR-FPR*FNR)/(sqrt(TPR+FPR)*sqrt(TPR+FNR)*sqrt(TNR+FPR)*sqrt(TNR+FNR))

print(paste("FNR =", FNR, ";  FPR =", FPR,";  MCC =", MCC ))

}
