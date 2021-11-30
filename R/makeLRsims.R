#' Make likelihoods ratio (LRs) simulations: a function for obtaining expected LRs under relatedness and unrelatedness kinship hypothesis.
#'
#' @param reference Reference pedigree. It could be an input from read_fam() function or a pedigree built with pedtools.
#' @param missing Missing person ID/label indicated in the pedigree.
#' @param numsims Number of simulations performed.
#' @param seed Select a seed for simulations. If it is defined, results will be reproducible. Suggested, seed = 123
#'
#' @return A dataframe with LRs obtained for both hypothesis, Unrelated where POI is not MP or Related where POI is MP.
#' @export
#' @import forrel
#' @import pedtools
#'
#' @examples
#' # A grandchild Search
#' # library(forrel)
#' # library(mispitools)
#' # x = linearPed(2)
#' # x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' # x = profileSim(x, N = 1, ids = 2)[[1]]
#' # plot(x)
#' # datasim = makeLRsims(x, missing = 5, 10, 123)




makeLRsims = function(reference, missing, numsims, seed) {
  st = base::Sys.time()

  if(!pedtools::is.ped(reference))
    base::stop("Expecting a connected pedigree as H1")

set.seed(seed)

### Person of interest 1: Unrelated
poi1 = pedtools::singleton("poi1")
# Transfer (empty) markers and simulate genotypes
poi1 = pedtools::transferMarkers(from = reference, to = poi1)
poi1 = forrel::profileSim(poi1, numsims)

lr1 <- as.list(rep(NA, numsims))
# Compute LR
for(i in 1:numsims) {
       lr1[[i]] = forrel::missingPersonLR(reference, missing, poi = poi1[[i]])
    }

### Person of interest 2: The true MP
# Simulate MP conditional on reference, and extract as singleton
poi2ped = forrel::profileSim(reference, numsims, ids = missing)

poi2 <- base::as.list(base::rep(NA, numsims))
# Extract MP as singleton
for(i in 1:numsims) {
  poi2[[i]] = base::subset(poi2ped[[i]], missing)
}

# free memory
base::rm(poi2ped)

#poi2 = subset(poi2, "MP")
#'
# Compute LR
lr2 <- base::as.list(rep(NA, numsims))

for(i in 1:numsims) {
  lr2[[i]] = forrel::missingPersonLR(reference, missing, poi = poi2[[i]])
}

# generating the dataframe
LRsimulated <- base::cbind(base::sapply(lr1, function(x) {x[["LRtotal"]][["H1:H2"]]}), base::sapply(lr2, function(x) {x[["LRtotal"]][["H1:H2"]]}))
base::colnames(LRsimulated) <- c("Unrelated", "Related")
base::structure(base::as.data.frame(LRsimulated))
}
