#' Simulate Genetic Profiles for Persons of Interest
#'
#' @description
#' Generates a database of simulated genetic profiles for persons of interest
#' (POIs) or unidentified human remains (UHRs). Profiles are randomly sampled
#' from population allele frequencies.
#'
#' @param numsims Integer. Number of genetic profiles to simulate.
#'   Default: 100.
#' @param reference A named list of allele frequencies in pedtools format.
#'   Can be created using \code{\link{get_allele_freqs}}.
#' @param seed Integer. Random seed for reproducibility. Default: 123.
#'
#' @return A data.frame where:
#'   \itemize{
#'     \item First column \code{id}: POI identifier (1 to numsims)
#'     \item Subsequent columns: Genetic marker data with allele pairs
#'   }
#'   Each row represents one simulated individual.
#'
#' @details
#' This function uses \pkg{pedtools} and \pkg{forrel} to generate random
#' genetic profiles based on the provided allele frequency database. The
#' profiles represent unrelated individuals sampled from the population.
#'
#' This is useful for:
#' \itemize{
#'   \item Creating test databases for simulation studies
#'   \item Generating H2 (unrelated) profiles for LR calculations
#'   \item Educational demonstrations of genetic variation
#' }
#'
#' @seealso
#' \code{\link{get_allele_freqs}} for preparing frequency data,
#' \code{\link{sim_lr_genetic}} for genetic LR simulations.
#'
#' @references
#' Marsico FL, Vigeland MD, Egeland T, Herrera Pinero F (2021). "Making
#' decisions in missing person identification cases with low statistical
#' power." \emph{Forensic Science International: Genetics}, 52, 102519.
#' \doi{10.1016/j.fsigen.2021.102519}
#'
#' @export
#' @import forrel
#' @import pedtools
#' @import dplyr
#' @examples
#' library(forrel)
#'
#' # Get frequency data
#' freqdata <- get_allele_freqs(Argentina)
#'
#' # Simulate 50 POI profiles
#' poi_db <- sim_poi_genetic(numsims = 50, reference = freqdata, seed = 123)
#' head(poi_db)
#'
#' # Check available markers
#' names(poi_db)

sim_poi_genetic <- function(numsims = 100, reference, seed = 123) {

  set.seed(seed)

  # Avoid R CMD check notes
  fid <- mid <- sex <- id <- NULL

  # Create singleton POI and set markers
  poi1 <- pedtools::singleton("poi1")
  poi1 <- pedtools::setMarkers(poi1, locusAttributes = reference)
  poi1 <- forrel::profileSim(poi1, numsims)

  # Convert nested lists to dataframes
  poi2 <- lapply(poi1, data.frame, stringsAsFactors = FALSE)
  poi <- do.call(rbind, poi2)
  poi <- dplyr::select(poi, -c(id, fid, sex, mid))

  # Add ID column
  id <- seq(1, numsims)
  poi <- cbind(id, poi)

  as.data.frame(poi)
}
