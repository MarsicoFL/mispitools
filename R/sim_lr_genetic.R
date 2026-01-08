#' Simulate Likelihood Ratios from Genetic Data
#'
#' @description
#' Simulates likelihood ratio (LR) distributions based on genetic (DNA) marker
#' data. This function generates expected LR distributions under two hypotheses:
#' \itemize{
#'   \item H1 (Related): The unidentified person IS the missing person
#'   \item H2 (Unrelated): The unidentified person is NOT the missing person
#' }
#'
#' This function wraps functionality from the \pkg{forrel} package to perform
#' missing person LR calculations using pedigree structures.
#'
#' @param reference A pedigree object with attached genetic markers. Can be
#'   created using \pkg{pedtools} functions like \code{linearPed()},
#'   \code{nuclearPed()}, etc., with markers attached via \code{setMarkers()}.
#' @param missing Character or numeric. The ID/label of the missing person
#'   in the pedigree.
#' @param numsims Integer. Number of simulations to perform. Default: 100.
#' @param seed Integer. Random seed for reproducibility. Default: 123.
#' @param numCores Integer. Number of CPU cores for parallel processing.
#'   Default: 1 (no parallelization).
#'
#' @return A list with two components:
#'   \itemize{
#'     \item \code{Unrelated}: List of LR objects from simulations where POI
#'       is unrelated to the pedigree (H2 simulations)
#'     \item \code{Related}: List of LR objects from simulations where POI
#'       is the actual missing person (H1 simulations)
#'   }
#'   Use \code{\link{lr_to_dataframe}} to convert this to a data.frame for
#'   further analysis.
#'
#' @details
#' The function performs two types of simulations:
#' \enumerate{
#'   \item \strong{H2 (Unrelated)}: Generates random genetic profiles for
#'     unrelated individuals using population allele frequencies, then
#'     calculates the LR for each profile.
#'   \item \strong{H1 (Related)}: Simulates genetic profiles for the actual
#'     missing person based on the pedigree structure, then calculates
#'     the LR for each profile.
#' }
#'
#' The LR is computed using \code{forrel::missingPersonLR()}, which calculates
#' the ratio of likelihoods: P(data | POI is MP) / P(data | POI is unrelated).
#'
#' @seealso
#' \code{\link{lr_to_dataframe}} for converting output to dataframe,
#' \code{\link{sim_lr_prelim}} for non-genetic LR simulations,
#' \code{\link{plot_lr_distribution}} for visualizing LR distributions,
#' \code{\link{decision_threshold}} for computing optimal thresholds.
#'
#' @references
#' Marsico FL, Vigeland MD, Egeland T, Herrera Pinero F (2021). "Making
#' decisions in missing person identification cases with low statistical
#' power." \emph{Forensic Science International: Genetics}, 52, 102519.
#' \doi{10.1016/j.fsigen.2021.102519}
#'
#' Vigeland MD, Egeland T (2021). "Joint DNA-based disaster victim
#' identification." \emph{Forensic Science International: Genetics}, 52, 102465.
#'
#' @export
#' @import forrel
#' @import pedtools
#'
#' @examples
#' library(forrel)
#' library(pedtools)
#'
#' # Create a simple pedigree: grandparent-parent-child
#' x <- linearPed(2)
#' plot(x)
#'
#' # Add genetic markers (using Norwegian frequencies as example)
#' x <- setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#'
#' # Simulate a profile for the reference person (ID 2)
#' x <- profileSim(x, N = 1, ids = 2)
#'
#' # Simulate LRs (person 5 is missing)
#' lr_sims <- sim_lr_genetic(x, missing = 5, numsims = 50, seed = 123)
#'
#' # Convert to dataframe for analysis
#' lr_df <- lr_to_dataframe(lr_sims)
#' head(lr_df)
#'
#' # Visualize distributions
#' plot_lr_distribution(lr_df)

sim_lr_genetic <- function(reference, missing, numsims = 100, seed = 123, numCores = 1) {

  st <- base::Sys.time()

  # Handle pedList input
if (pedtools::is.pedList(reference) && base::length(reference) == 1) {
    reference <- reference[[1]]
  }

  # Validate input
  if (!pedtools::is.ped(reference)) {
    base::stop("Expecting a connected pedigree as H1")
  }

  set.seed(seed)

  # === H2 Simulations (Unrelated) ===
  # Create singleton POI and transfer marker information
  poi1 <- pedtools::singleton("poi1")
  poi1 <- pedtools::transferMarkers(from = reference, to = poi1)
  poi1 <- forrel::profileSim(poi1, numsims, numCores = numCores)

  # Calculate LR for each unrelated profile
  lr1 <- as.list(rep(NA, numsims))
  for (i in 1:numsims) {
    lr1[[i]] <- forrel::missingPersonLR(reference, missing, poi = poi1[[i]])
  }

  # === H1 Simulations (Related) ===
  # Simulate profiles for the actual missing person
  poi2ped <- forrel::profileSim(reference, numsims, ids = missing, numCores = numCores)

  # Extract the missing person as singleton
  poi2 <- base::as.list(base::rep(NA, numsims))
  for (i in 1:numsims) {
    poi2[[i]] <- base::subset(poi2ped[[i]], missing)
  }

  base::rm(poi2ped)

  # Calculate LR for each related profile
  lr2 <- base::as.list(rep(NA, numsims))
  for (i in 1:numsims) {
    lr2[[i]] <- forrel::missingPersonLR(reference, missing, poi = poi2[[i]])
  }

  # Return as list structure
  base::structure(list(Unrelated = lr1, Related = lr2))
}
