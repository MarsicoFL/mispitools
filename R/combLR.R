#' Combine LRs: a function for combining LRs obtained from simulations.
#'
#' @param LRdatasim1 A data frame object with the results of simulations. Outputs from simLRgen or simLRprelim funcionts.
#' @param LRdatasim2 A second data frame object with the results of simulations. Outputs from simLRgen or simLRprelim funcionts.
#'
#' @return An object of class data.frame combining the LRs obtained from simulations (the function multiplies the LRs).  
#' @export
#' @examples
#' library(mispitools)
#' library(forrel) 
#' x = linearPed(2)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' LRdatasim1 = simLRgen(x, missing = 5, 10, 123)
#' LRdatasim2 = simLRprelim("sex")
#' combLR(LRdatasim1,LRdatasim2)

combLR = function(LRdatasim1, LRdatasim2) {
a <- as.numeric(unlist(LRdatasim1$Unrelated))*as.numeric(unlist(LRdatasim2$Unrelated))
b <- as.numeric(unlist(LRdatasim1$Related))*as.numeric(unlist(LRdatasim2$Related))

LRsimulated <- base::cbind(a,b)
base::colnames(LRsimulated) <- c("Unrelated", "Related")
base::structure(base::as.data.frame(LRsimulated))
}
