#' Generate Reference Properties for a Hypothetical Population
#'
#' This function simulates a dataset representing physical characteristics
#' (hair color, skin color, eye color) of a hypothetical population,
#' based on conditional probability distributions. The size of the simulated population
#' can be adjusted by the user.
#'
#' @param n The number of individuals in the simulated population.
#' @param seed Selected seed for simulations.
#' @return A \code{data.frame} with three columns: hair_colour, skin_colour, and eye_colour,
#' each representing the respective characteristics of each individual in the sample population.
#' The hair color is simulated based on predefined probabilities, and skin and eye colors
#' are generated conditionally based on the hair color.
#'
#' @examples
#' simRef(1000) # Generates a data frame with 1000 entries based on the defined distributions.
#' @export
simRef <- function(n = 1000, seed = 1234) {
  set.seed(seed)

  prob_hair <- c(0.4, 0.3, 0.15, 0.1, 0.05)

  hair_colour <- sample(1:5, n, replace = TRUE, prob = prob_hair)

  skin_colour <- numeric(n)
  eye_colour <- numeric(n)

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
  eye_colour = eye_colour)

  return(data)

 }
