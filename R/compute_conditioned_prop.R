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

  # Input validation

  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  required_cols <- c("hair_colour", "skin_colour", "eye_colour")
  if (!all(required_cols %in% names(data))) {
    stop("data must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (!is.numeric(h) || h < 1 || h > 5) {
    stop("h (hair color) must be an integer between 1 and 5")
  }
  if (!is.numeric(s) || s < 1 || s > 5) {
    stop("s (skin color) must be an integer between 1 and 5")
  }
  if (!is.numeric(y) || y < 1 || y > 5) {
    stop("y (eye color) must be an integer between 1 and 5")
  }
  if (!is.numeric(eh) || eh < 0 || eh > 1) {
    stop("eh (hair error rate) must be between 0 and 1")
  }
  if (!is.numeric(es) || es < 0 || es > 1) {
    stop("es (skin error rate) must be between 0 and 1")
  }
  if (!is.numeric(ey) || ey < 0 || ey > 1) {
    stop("ey (eye error rate) must be between 0 and 1")
  }

  numerators <- numeric(nrow(data))

  # Probability of correct observation = (1 - error_rate)
  # Probability of incorrect observation = error_rate
  # Assuming independence across traits
  p_h_match <- 1 - eh
  p_h_mismatch <- eh
  p_s_match <- 1 - es
  p_s_mismatch <- es
  p_y_match <- 1 - ey
  p_y_mismatch <- ey

  for (i in 1:nrow(data)) {
    match_h <- data$hair_colour[i] == h
    match_s <- data$skin_colour[i] == s
    match_y <- data$eye_colour[i] == y

    # Calculate probability as product of independent match/mismatch probabilities
    p_h <- if (match_h) p_h_match else p_h_mismatch
    p_s <- if (match_s) p_s_match else p_s_mismatch
    p_y <- if (match_y) p_y_match else p_y_mismatch

    numerators[i] <- p_h * p_s * p_y
  }

  probs <- as.data.frame(cbind(data, numerators))
  probs <- unique(probs)

  # Normalize to sum to 1
  total <- sum(probs$numerators)
  if (total == 0) {
    warning("All numerators are zero; returning unnormalized values")
  } else {
    probs$numerators <- probs$numerators / total
  }

  return(probs)
}
