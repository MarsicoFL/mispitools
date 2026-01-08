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

  # Hair color probabilities
  prob_hair <- c(0.4, 0.3, 0.15, 0.1, 0.05)
  hair_colour <- sample(1:5, n, replace = TRUE, prob = prob_hair)

  skin_colour <- numeric(n)
  eye_colour <- numeric(n)

  # Skin color conditional on hair color
  for (i in 1:n) {
    if (hair_colour[i] == 1) {
      skin_colour[i] <- sample(1:5, 1, prob = c(0.1, 0.15, 0.4, 0.25, 0.1))
    } else if (hair_colour[i] == 2) {
      skin_colour[i] <- sample(1:5, 1, prob = c(0.1, 0.2, 0.4, 0.2, 0.1))
    } else if (hair_colour[i] == 3) {
      skin_colour[i] <- sample(1:5, 1, prob = c(0.6, 0.3, 0.05, 0.03, 0.02))
    } else if (hair_colour[i] == 4) {
      skin_colour[i] <- sample(1:5, 1, prob = c(0.7, 0.2, 0.05, 0.03, 0.02))
    } else {
      skin_colour[i] <- sample(1:5, 1, prob = c(0.2, 0.3, 0.2, 0.2, 0.1))
    }
  }

  # Eye color conditional on hair and skin color
  for (i in 1:n) {
    if (hair_colour[i] == 1 & skin_colour[i] == 1) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.5, 0.2, 0.1, 0.1, 0.1))
    } else if (hair_colour[i] == 1 & skin_colour[i] == 2) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.6, 0.1, 0.1, 0.1, 0.1))
    } else if (hair_colour[i] == 1 & skin_colour[i] == 3) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.7, 0.1, 0.1, 0.05, 0.05))
    } else if (hair_colour[i] == 1 & skin_colour[i] == 4) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.8, 0.05, 0.05, 0.05, 0.05))
    } else if (hair_colour[i] == 1 & skin_colour[i] == 5) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.9, 0.025, 0.025, 0.025, 0.025))
    } else if (hair_colour[i] == 2 & skin_colour[i] == 1) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.3, 0.3, 0.2, 0.1, 0.1))
    } else if (hair_colour[i] == 2 & skin_colour[i] == 2) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.4, 0.3, 0.1, 0.1, 0.1))
    } else if (hair_colour[i] == 2 & skin_colour[i] == 3) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.5, 0.2, 0.1, 0.1, 0.1))
    } else if (hair_colour[i] == 2 & skin_colour[i] == 4) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.6, 0.15, 0.1, 0.1, 0.05))
    } else if (hair_colour[i] == 2 & skin_colour[i] == 5) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.7, 0.1, 0.1, 0.05, 0.05))
    } else if (hair_colour[i] == 3 & skin_colour[i] == 1) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.1, 0.6, 0.1, 0.1, 0.1))
    } else if (hair_colour[i] == 3 & skin_colour[i] == 2) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.1, 0.7, 0.1, 0.05, 0.05))
    } else if (hair_colour[i] == 3 & skin_colour[i] == 3) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.15, 0.6, 0.1, 0.1, 0.05))
    } else if (hair_colour[i] == 3 & skin_colour[i] == 4) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.2, 0.5, 0.1, 0.1, 0.1))
    } else if (hair_colour[i] == 3 & skin_colour[i] == 5) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.25, 0.5, 0.1, 0.1, 0.05))
    } else if (hair_colour[i] == 4 & skin_colour[i] == 1) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.1, 0.7, 0.05, 0.1, 0.05))
    } else if (hair_colour[i] == 4 & skin_colour[i] == 2) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.1, 0.6, 0.1, 0.1, 0.1))
    } else if (hair_colour[i] == 4 & skin_colour[i] == 3) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.1, 0.7, 0.05, 0.1, 0.05))
    } else if (hair_colour[i] == 4 & skin_colour[i] == 4) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.15, 0.6, 0.1, 0.1, 0.05))
    } else if (hair_colour[i] == 4 & skin_colour[i] == 5) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.2, 0.6, 0.05, 0.1, 0.05))
    } else if (hair_colour[i] == 5 & skin_colour[i] == 1) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.4, 0.2, 0.15, 0.15, 0.1))
    } else if (hair_colour[i] == 5 & skin_colour[i] == 2) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.5, 0.2, 0.1, 0.1, 0.1))
    } else if (hair_colour[i] == 5 & skin_colour[i] == 3) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.5, 0.15, 0.15, 0.1, 0.1))
    } else if (hair_colour[i] == 5 & skin_colour[i] == 4) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.6, 0.15, 0.1, 0.1, 0.05))
    } else if (hair_colour[i] == 5 & skin_colour[i] == 5) {
      eye_colour[i] <- sample(1:5, 1, prob = c(0.4, 0.1, 0.2, 0.2, 0.1))
    }
  }

  data <- data.frame(
    hair_colour = hair_colour,
    skin_colour = skin_colour,
    eye_colour = eye_colour
  )

  return(data)
}
