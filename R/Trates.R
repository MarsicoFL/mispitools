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

unrelated_values <- vector()
related_values <- vector()

list_length <- length(datasim[["Unrelated"]])

for (i in 1:list_length) {
  unrelated_value <- datasim[["Unrelated"]][[i]][["LRtotal"]][["H1:H2"]]
  related_value <- datasim[["Related"]][[i]][["LRtotal"]][["H1:H2"]]
  
  unrelated_values <- c(unrelated_values, unrelated_value)
  related_values <- c(related_values, related_value)
}

results_df <- data.frame(Unrelated = unrelated_values, Related = related_values)
datasim <- results_df

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
