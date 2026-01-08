#' Compute Likelihood Ratios for Pigmentation Traits
#'
#' @description
#' Computes likelihood ratios (LRs) for each unique combination of hair color,
#' skin color, and eye color by dividing conditioned proportions (numerators)
#' by reference proportions (denominators).
#'
#' @param conditioned A data.frame with columns \code{hair_colour},
#'   \code{skin_colour}, \code{eye_colour}, and \code{numerators}.
#'   Typically output from \code{\link{compute_conditioned_prop}}.
#' @param unconditioned A data.frame with columns \code{hair_colour},
#'   \code{skin_colour}, \code{eye_colour}, and \code{f_h_s_y}.
#'   Typically output from \code{\link{compute_reference_prop}}.
#'
#' @return A data.frame with:
#'   \itemize{
#'     \item \code{hair_colour}, \code{skin_colour}, \code{eye_colour}: Trait combination
#'     \item \code{f_h_s_y}: Population frequency (denominator)
#'     \item \code{numerators}: Conditioned probability (numerator)
#'     \item \code{LR}: Likelihood ratio = numerators / f_h_s_y
#'   }
#'   Combinations not present in both inputs are excluded.
#'
#' @seealso
#' \code{\link{compute_conditioned_prop}} for computing numerators,
#' \code{\link{compute_reference_prop}} for computing denominators,
#' \code{\link{lr_pigmentation}} for simulating LR distributions.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Generate population data
#' pop_data <- sim_reference_pop(n = 500, seed = 123)
#'
#' # Compute proportions
#' conditioned <- compute_conditioned_prop(pop_data, 1, 1, 1, 0.01, 0.01, 0.01)
#' unconditioned <- compute_reference_prop(pop_data)
#'
#' # Compute LRs
#' lrs <- lr_compute_pigmentation(conditioned, unconditioned)
#' head(lrs)
#'
#' # Highest LRs (most discriminating combinations)
#' lrs[order(-lrs$LR), ][1:5, ]

lr_compute_pigmentation <- function(conditioned, unconditioned) {

  merged_data <- merge(unconditioned, conditioned,
                       by = c("hair_colour", "skin_colour", "eye_colour"))

  merged_data$LR <- merged_data$numerators / merged_data$f_h_s_y

  return(merged_data)
}
