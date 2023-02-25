#' Simulate likelihoods ratio (LRs) based on genetic data: a function for obtaining expected LRs under relatedness and unrelatedness kinship hypothesis.
#'
#' @param reference Reference pedigree. It could be an input from read_fam() function or a pedigree built with pedtools.
#' @param missing Missing person ID/label indicated in the pedigree.
#' @param numsims Number of simulations performed.
#' @param seed Select a seed for simulations. If it is defined, results will be reproducible. Suggested, seed = 123
#'
#' @return An object of class data.frame with LRs obtained for both hypothesis, Unrelated where POI is not MP or Related where POI is MP.
#' @export
#' @import forrel
#' @import pedtools
#'
#' @examples
#' library(forrel)
#' x = linearPed(2)
#' plot(x)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' datasim = simLRgen(x, missing = 5, 10, 123)




simLRgen = function(reference, missing, numsims, seed) {
  st = base::Sys.time()

  if(pedtools::is.pedList(reference) && base::length(reference) == 1)
    reference = reference[[1]]

  if(!pedtools::is.ped(reference))
    base::stop("Expecting a connected pedigree as H1")

set.seed(seed)

poi1 = pedtools::singleton("poi1")
poi1 = pedtools::transferMarkers(from = reference, to = poi1)
poi1 = forrel::profileSim(poi1, numsims)

lr1 <- as.list(rep(NA, numsims))
for(i in 1:numsims) {
       lr1[[i]] = forrel::missingPersonLR(reference, missing, poi = poi1[[i]])
    }

poi2ped = forrel::profileSim(reference, numsims, ids = missing)

poi2 <- base::as.list(base::rep(NA, numsims))
for(i in 1:numsims) {
  poi2[[i]] = base::subset(poi2ped[[i]], missing)
}

base::rm(poi2ped)

lr2 <- base::as.list(rep(NA, numsims))

for(i in 1:numsims) {
  lr2[[i]] = forrel::missingPersonLR(reference, missing, poi = poi2[[i]])
}

LRsimulated <- base::cbind(base::sapply(lr1, function(x) {x[["LRtotal"]][["H1:H2"]]}), base::sapply(lr2, function(x) {x[["LRtotal"]][["H1:H2"]]}))
base::colnames(LRsimulated) <- c("Unrelated", "Related")
base::structure(base::as.data.frame(LRsimulated))
}
