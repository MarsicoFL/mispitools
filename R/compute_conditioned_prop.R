#' Compute Conditioned Proportions for Pigmentation Traits
#'
#' @description
#' Calculates the conditioned proportions (numerator probabilities) for
#' pigmentation trait combinations given the missing person's characteristics.
#' These proportions represent P(observed traits | H1), accounting for
#' observation errors.
#'
#' @param data A data.frame with columns \code{hair_colour}, \code{skin_colour},
#'   and \code{eye_colour}, typically output from \code{\link{sim_reference_pop}}.
#' @param h Integer (1-5). Missing person's hair color category.
#' @param s Integer (1-5). Missing person's skin color category.
#' @param y Integer (1-5). Missing person's eye color category.
#' @param eh Numeric (0-1). Error rate for observing hair color.
#' @param es Numeric (0-1). Error rate for observing skin color.
#' @param ey Numeric (0-1). Error rate for observing eye color.
#'
#' @return A data.frame with the original trait columns plus:
#'   \itemize{
#'     \item \code{numerators}: Conditioned probability for each combination,
#'       normalized to sum to 1
#'   }
#'   Only unique combinations are returned.
#'
#' @details
#' The function calculates the probability of observing each trait combination
#' given the MP's true characteristics and the error rates. Higher probabilities
#' are assigned to combinations matching the MP's traits, while combinations
#' with mismatches have probabilities proportional to the error rates.
#'
#' @seealso
#' \code{\link{sim_reference_pop}} for generating the input data,
#' \code{\link{compute_reference_prop}} for reference proportions,
#' \code{\link{lr_compute_pigmentation}} for computing LRs.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Generate reference population
#' pop_data <- sim_reference_pop(n = 500, seed = 123)
#'
#' # Compute conditioned proportions for MP with traits (1,1,1)
#' cond_prop <- compute_conditioned_prop(
#'   pop_data,
#'   h = 1, s = 1, y = 1,
#'   eh = 0.01, es = 0.01, ey = 0.01
#' )
#' head(cond_prop)

compute_conditioned_prop <- function(data, h, s, y, eh, es, ey) {

  numerators <- numeric(nrow(data))

  for (i in 1:nrow(data)) {
    C_h <- as.integer(data$hair_colour[i] == h)
    C_s <- as.integer(data$skin_colour[i] == s)
    C_y <- as.integer(data$eye_colour[i] == y)

    if (C_h && C_s && C_y) {
      numerators[i] <- 1 - eh - es - ey
    } else if (C_h && C_s) {
      numerators[i] <- (1 - ey) * eh * es
    } else if (C_h && C_y) {
      numerators[i] <- (1 - es) * eh * ey
    } else if (C_s && C_y) {
      numerators[i] <- (1 - eh) * es * ey
    } else if (C_h) {
      numerators[i] <- (1 - es - ey) * eh
    } else if (C_s) {
      numerators[i] <- (1 - eh - ey) * es
    } else if (C_y) {
      numerators[i] <- (1 - eh - es) * ey
    } else {
      numerators[i] <- eh * es * ey
    }
  }

  probs <- as.data.frame(cbind(data, numerators))
  probs <- unique(probs)
  probs$numerators <- probs$numerators / sum(probs$numerators)

  return(probs)
}
