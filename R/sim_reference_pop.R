#' Simulate Reference Population with Pigmentation Traits
#'
#' @description
#' Generates a simulated population dataset with correlated pigmentation
#' characteristics (hair color, skin color, eye color). The traits are
#' simulated using conditional probability distributions that reflect
#' realistic correlations between these characteristics.
#'
#' @param n Integer. Number of individuals to simulate. Default: 1000.
#' @param seed Integer. Random seed for reproducibility. Default: 1234.
#'
#' @return A data.frame with three columns:
#'   \itemize{
#'     \item \code{hair_colour}: Hair color category (1-5)
#'     \item \code{skin_colour}: Skin color category (1-5)
#'     \item \code{eye_colour}: Eye color category (1-5)
#'   }
#'   Categories are numbered 1 (lightest) to 5 (darkest).
#'
#' @details
#' Hair color categories:
#' \enumerate{
#'   \item Blonde/Light
#'   \item Light brown
#'   \item Medium brown
#'   \item Dark brown
#'   \item Black
#' }
#'
#' The simulation uses conditional probability distributions where:
#' \itemize{
#'   \item Hair color is sampled first from population frequencies
#'   \item Skin color is sampled conditional on hair color
#'   \item Eye color is sampled conditional on both hair and skin color
#' }
#'
#' This captures realistic correlations (e.g., darker hair tends to
#' co-occur with darker skin and eyes).
#'
#' @seealso
#' \code{\link{compute_conditioned_prop}} for computing proportions,
#' \code{\link{compute_reference_prop}} for reference frequencies,
#' \code{\link{lr_pigmentation}} for pigmentation LR calculations.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Simulate a population of 500 individuals
#' pop_data <- sim_reference_pop(n = 500, seed = 123)
#' head(pop_data)
#'
#' # Check trait distributions
#' table(pop_data$hair_colour)
#' table(pop_data$skin_colour)
#' table(pop_data$eye_colour)
#'
#' # Use for LR calculations
#' conditioned <- compute_conditioned_prop(pop_data, h = 1, s = 1, y = 1,
#'                                         eh = 0.01, es = 0.01, ey = 0.01)
#' unconditioned <- compute_reference_prop(pop_data)

sim_reference_pop <- function(n = 1000, seed = 1234) {

  set.seed(seed)

  # Input validation
  if (!is.numeric(n) || n < 1) {
    stop("n must be a positive integer")
  }
  n <- as.integer(n)

  # ==========================================================================
  # Conditional probability matrices for pigmentation traits
  # These capture correlations between hair, skin, and eye color
  # Rows = conditioning variable, Columns = outcome probabilities (categories 1-5)
  # ==========================================================================

  # P(hair color) - marginal distribution
  prob_hair <- c(0.4, 0.3, 0.15, 0.1, 0.05)

  # P(skin | hair) - 5x5 matrix, rows = hair color, cols = skin color probs
  prob_skin_given_hair <- matrix(
    c(0.10, 0.15, 0.40, 0.25, 0.10,  # hair = 1 (blonde)
      0.10, 0.20, 0.40, 0.20, 0.10,  # hair = 2 (light brown)
      0.60, 0.30, 0.05, 0.03, 0.02,  # hair = 3 (medium brown)
      0.70, 0.20, 0.05, 0.03, 0.02,  # hair = 4 (dark brown)
      0.20, 0.30, 0.20, 0.20, 0.10), # hair = 5 (black)
    nrow = 5, ncol = 5, byrow = TRUE
  )

  # P(eye | hair, skin) - stored as list of 5x5 matrices
  # prob_eye_given_hair_skin[[hair]][skin, eye]
  prob_eye_given_hair_skin <- list(
    # hair = 1 (blonde)
    matrix(c(
      0.50, 0.20, 0.10, 0.10, 0.10,   # skin = 1
      0.60, 0.10, 0.10, 0.10, 0.10,   # skin = 2
      0.70, 0.10, 0.10, 0.05, 0.05,   # skin = 3
      0.80, 0.05, 0.05, 0.05, 0.05,   # skin = 4
      0.90, 0.025, 0.025, 0.025, 0.025 # skin = 5
    ), nrow = 5, ncol = 5, byrow = TRUE),
    # hair = 2 (light brown)
    matrix(c(
      0.30, 0.30, 0.20, 0.10, 0.10,
      0.40, 0.30, 0.10, 0.10, 0.10,
      0.50, 0.20, 0.10, 0.10, 0.10,
      0.60, 0.15, 0.10, 0.10, 0.05,
      0.70, 0.10, 0.10, 0.05, 0.05
    ), nrow = 5, ncol = 5, byrow = TRUE),
    # hair = 3 (medium brown)
    matrix(c(
      0.10, 0.60, 0.10, 0.10, 0.10,
      0.10, 0.70, 0.10, 0.05, 0.05,
      0.15, 0.60, 0.10, 0.10, 0.05,
      0.20, 0.50, 0.10, 0.10, 0.10,
      0.25, 0.50, 0.10, 0.10, 0.05
    ), nrow = 5, ncol = 5, byrow = TRUE),
    # hair = 4 (dark brown)
    matrix(c(
      0.10, 0.70, 0.05, 0.10, 0.05,
      0.10, 0.60, 0.10, 0.10, 0.10,
      0.10, 0.70, 0.05, 0.10, 0.05,
      0.15, 0.60, 0.10, 0.10, 0.05,
      0.20, 0.60, 0.05, 0.10, 0.05
    ), nrow = 5, ncol = 5, byrow = TRUE),
    # hair = 5 (black)
    matrix(c(
      0.40, 0.20, 0.15, 0.15, 0.10,
      0.50, 0.20, 0.10, 0.10, 0.10,
      0.50, 0.15, 0.15, 0.10, 0.10,
      0.60, 0.15, 0.10, 0.10, 0.05,
      0.40, 0.10, 0.20, 0.20, 0.10
    ), nrow = 5, ncol = 5, byrow = TRUE)
  )

  # ==========================================================================
  # Vectorized sampling
  # ==========================================================================

  # Sample hair color
  hair_colour <- sample(1:5, n, replace = TRUE, prob = prob_hair)

  # Sample skin color conditional on hair (vectorized by hair category)
  skin_colour <- integer(n)
  for (h in 1:5) {
    idx <- which(hair_colour == h)
    if (length(idx) > 0) {
      skin_colour[idx] <- sample(1:5, length(idx), replace = TRUE,
                                  prob = prob_skin_given_hair[h, ])
    }
  }

  # Sample eye color conditional on hair and skin (vectorized by combination)
  eye_colour <- integer(n)
  for (h in 1:5) {
    for (s in 1:5) {
      idx <- which(hair_colour == h & skin_colour == s)
      if (length(idx) > 0) {
        eye_colour[idx] <- sample(1:5, length(idx), replace = TRUE,
                                   prob = prob_eye_given_hair_skin[[h]][s, ])
      }
    }
  }

  data.frame(
    hair_colour = hair_colour,
    skin_colour = skin_colour,
    eye_colour = eye_colour
  )
}
