#' Make preliminary investigation data simulations: a function for obtaining a database of preliminary investigation data for a missing child search.
#'
#' @param dateinit Minimun birth date of simulated persons of interest.
#' @param scenario Birth date distribution scenarios: (1) non-uniform, (2) uniform.
#' @param seed Select a seed for simulations. If it is defined, results will be reproducible. Suggested, seed = 123
#' @param ext Time extension for minimun birth date, rang in scenario 1 and days in scenario 2.
#' @param femaleprop Proportion of females.
#' @param numsims Number of simulated POIs.
#' @param region Birth region or place, argentine regions are listed as default. 
#' @param regionprob Birth region or place proportions.
#' @param birthprob Birth type probabilities: home birth, hospital birth and unknown-adoption.
#'
#' @return An object of class data.frame with preliminary investigation data.
#' @export
#'
#' @examples
#' makePrelim(
#'   dateinit = "1975/01/01",
#'   scenario = 1,
#'   femaleprop = 0.5,
#'   ext = 100,
#'   numsims = 10000,
#'   seed = 123,
#'   birthprob = c(0.09, 0.9, 0.01),
#'   region = c("Cuyo", "Patagonia", "Central region", "North west region", "Litoral",
#'             "Buenos Aires"),
#'   regionprob = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.5)
#'   )





makePrelim = function(dateinit = "1975/01/01", scenario = 1, femaleprop = 0.5, ext = 100, numsims = 10000, seed = 123, birthprob = c(0.09, 0.9, 0.01), region = c("Cuyo", "Patagonia", "Central region", "North west region", "Litoral", "Buenos Aires"), regionprob = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.5)) {

gender = c("female","male")
maleprop = 1 - femaleprop
birth = c("home birth", "hospital birth", "unknown-adoption")
set.seed(seed)

a <- seq(1, numsims, by = 1)

if(scenario == 1) {
  b <- as.data.frame(as.Date('1975/01/01') + ext*rgamma(numsims, 12))
}

else if (scenario == 2) {
  b <- as.data.frame(as.Date('1975/01/01') + runif(numsims, min = 0, max = ext))
}

  c <- sample(gender, numsims, replace = TRUE, prob = c(femaleprop, maleprop))
  d <- sample(birth, numsims, replace = TRUE, prob = birthprob)
  e <- sample(region, numsims, replace = TRUE, prob = regionprob)

PrelimDatasim <- cbind(a, b, c, d, e)
base::colnames(PrelimDatasim) <- c("POI-ID", "DBD", "Gender", "Birth-type", "Birth place")
base::structure(base::as.data.frame(PrelimDatasim))
}
