#' Make POIs gen: a function for obtaining a database with genetic information from simulated POIs or UHRs.
#'
#' @param reference Indicate the reference STRs/SNPs frequency database used for simulations.
#' @param numsims Number of simulations performed (numer of POIs or UHRs).
#' @param seed Select a seed for simulations. If it is defined, results will be reproducible. Suggested, seed = 123
#'
#' @return An object of class data.frame with genetic information from POIs (randomly sampled from the frequency database).
#' @export
#' @import forrel
#' @import pedtools
#' @examples
#' library(forrel) 
#' freqdata <- getfreqs(Argentina)
#' makePOIgen(numsims = 100, reference = freqdata, seed = 123)




makePOIgen = function(numsims = 100, reference ,seed = 123) {
set.seed(seed)
fid <- mid <- sex <- NULL

poi1 = pedtools::singleton("poi1")
poi1 = pedtools::setMarkers(poi1, locusAttributes = reference)
poi1 = forrel::profileSim(poi1, numsims)

poi2<- lapply(poi1, data.frame, stringsAsFactors = FALSE) # convertir los nested lists en dataframes.
poi <- do.call(rbind, poi2) # combinalos en uno solo
poi <- select(poi, -c(id, fid, sex, mid))
id <- seq(1:numsims)
poi <- cbind(id, poi)

base::structure(base::as.data.frame(poi))
}

