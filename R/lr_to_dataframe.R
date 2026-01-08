#' Convert Genetic LR Simulations to Data Frame
#'
#' @description
#' Converts the list output from \code{\link{sim_lr_genetic}} into a tidy
#' data frame suitable for analysis and visualization. Extracts the total
#' LR values from each simulation.
#'
#' @param datasim A list object returned by \code{\link{sim_lr_genetic}},
#'   containing \code{Unrelated} and \code{Related} components with LR objects.
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{Unrelated}: Numeric LR values from H2 simulations
#'     \item \code{Related}: Numeric LR values from H1 simulations
#'   }
#'   The number of rows equals the number of simulations.
#'
#' @details
#' The function extracts \code{LRtotal[["H1:H2"]]} from each LR object
#' in the simulation lists. This represents the overall likelihood ratio
#' across all genetic markers.
#'
#' @seealso
#' \code{\link{sim_lr_genetic}} for generating the input,
#' \code{\link{plot_lr_distribution}} for visualization,
#' \code{\link{lr_combine}} for combining with other LR sources.
#'
#' @references
#' Marsico FL, Vigeland MD, Egeland T, Herrera Pinero F (2021). "Making
#' decisions in missing person identification cases with low statistical
#' power." \emph{Forensic Science International: Genetics}, 52, 102519.
#' \doi{10.1016/j.fsigen.2021.102519}
#'
#' @export
#' @examples
#' library(forrel)
#'
#' # Create pedigree and simulate
#' x <- linearPed(2)
#' x <- setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x <- profileSim(x, N = 1, ids = 2)
#'
#' # Simulate LRs
#' lr_sims <- sim_lr_genetic(x, missing = 5, numsims = 50, seed = 123)
#'
#' # Convert to dataframe
#' lr_df <- lr_to_dataframe(lr_sims)
#' head(lr_df)
#'
#' # Now can use with other functions
#' summary(log10(lr_df$Related))
#' plot_lr_distribution(lr_df)

lr_to_dataframe <- function(datasim) {

  # Input validation
  if (!is.list(datasim)) {
    stop("datasim must be a list (typically from sim_lr_genetic)")
  }

  if (!all(c("Related", "Unrelated") %in% names(datasim))) {
    stop("datasim must have 'Related' and 'Unrelated' components")
  }

  if (length(datasim$Related) == 0 || length(datasim$Unrelated) == 0) {
    stop("datasim$Related and datasim$Unrelated must not be empty")
  }

  # Helper function to safely extract LR value
  .extract_lr <- function(x, component_name) {
    if (!is.list(x)) {
      stop(sprintf("Invalid LR object in %s: expected list, got %s",
                   component_name, class(x)[1]))
    }
    if (!"LRtotal" %in% names(x)) {
      stop(sprintf("LR object in %s missing 'LRtotal' component", component_name))
    }
    if (!"H1:H2" %in% names(x[["LRtotal"]])) {
      stop(sprintf("LRtotal in %s missing 'H1:H2' value", component_name))
    }
    x[["LRtotal"]][["H1:H2"]]
  }

  # Extract LR values from Related simulations (H1)
  Related <- tryCatch(
    sapply(datasim$Related, function(x) .extract_lr(x, "Related")),
    error = function(e) {
      stop("Error extracting Related LR values: ", conditionMessage(e))
    }
  )

  # Extract LR values from Unrelated simulations (H2)
  Unrelated <- tryCatch(
    sapply(datasim$Unrelated, function(x) .extract_lr(x, "Unrelated")),
    error = function(e) {
      stop("Error extracting Unrelated LR values: ", conditionMessage(e))
    }
  )

  # Validate extracted values
  if (any(is.na(Related))) {
    warning(sprintf("%d NA values in Related LRs", sum(is.na(Related))))
  }
  if (any(is.na(Unrelated))) {
    warning(sprintf("%d NA values in Unrelated LRs", sum(is.na(Unrelated))))
  }

  # Combine into dataframe
  LRsimulated <- cbind(Unrelated, Related)
  colnames(LRsimulated) <- c("Unrelated", "Related")

  as.data.frame(LRsimulated)
}
