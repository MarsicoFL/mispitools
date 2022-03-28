#' Make preliminary investigation POI/UHR data simulations: a function for obtaining a database of preliminary investigation data for a missing person search.
#'
#' @param casetype Type of missing person search case. Two options are available: "migrants" or "children".
#' @param dateinit Minimun birth date of simulated persons of interest. Casetype: Children.
#' @param scenario Birth date distribution scenarios: (1) non-uniform, (2) uniform. Casetype: Children.
#' @param seed Select a seed for simulations. If it is defined, results will be reproducible. Casetype: All.
#' @param ext Time extension for minimun birth date, range in scenario 1 and days in scenario 2. Casetype: Children.
#' @param femaleprop Proportion of females. Casetype: All.
#' @param numsims Number of simulated POIs/UHRs. Casetype: All.
#' @param region Birth region or place in missing children case or place of discovery of the human remain in missing migrant case. Casetype: All.
#' @param regionprob Region proportions. Casetype: All.
#' @param birthprob Birth type probabilities: home birth, hospital birth and unknown-adoption. Casetype: Children.
#'
#' @return An object of class data.frame with preliminary investigation data.
#' @export
#' @importFrom stats rgamma
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @examples
#' makePOIprelim(
#'   dateinit = "1975/01/01",
#'   scenario = 1,
#'   femaleprop = 0.5,
#'   ext = 100,
#'   numsims = 10000,
#'   seed = 123,
#'   birthprob = c(0.09, 0.9, 0.01),
#'   region = c("North America", "South America", "Africa", "Asia", "Europe", "Oceania"),
#'   regionprob = c(0.2, 0.2, 0.2, 0.1, 0.2, 0.1))
#'   


makePOIprelim = function(casetype = "children", dateinit = "1975/01/01", scenario = 1, femaleprop = 0.5, ext = 100, numsims = 10000, seed = 123, birthprob = c(0.09, 0.9, 0.01), region = c("North America", "South America", "Africa", "Asia", "Europe", "Oceania"), regionprob = c(0.2, 0.2, 0.2, 0.1, 0.2, 0.1)) {

if(casetype == "children") {
sex = c("female","male")
maleprop = 1 - femaleprop
birth = c("home birth", "hospital birth", "unknown-adoption")
set.seed(seed)

a <- seq(1, numsims, by = 1)

if(scenario == 1) {
  b <- as.data.frame(as.Date(dateinit) + ext*rgamma(numsims, 12))
}

else if (scenario == 2) {
  b <- as.data.frame(as.Date(dateinit) + runif(numsims, min = 0, max = ext))
}

  c <- sample(sex, numsims, replace = TRUE, prob = c(femaleprop, maleprop))
  d <- sample(birth, numsims, replace = TRUE, prob = birthprob)
  e <- sample(region, numsims, replace = TRUE, prob = regionprob)

PrelimDatasim <- cbind(a, b, c, d, e)
base::colnames(PrelimDatasim) <- c("POI-ID", "DBD", "Sex", "Birth-type", "Region")
base::structure(base::as.data.frame(PrelimDatasim))}

else if (casetype == "migrants") { 
sex = c("female","male")
maleprop = 1 - femaleprop
age <- sample(seq.int(from = 4, to = 70, by = 1), size = numsims, replace = TRUE)
height <- rnorm(numsims, mean=170, sd=15)

set.seed(seed)

a <- seq(1, numsims, by = 1)
c <- sample(sex, numsims, replace = TRUE, prob = c(femaleprop, maleprop))
e <- sample(region, numsims, replace = TRUE, prob = regionprob)

PrelimDatasim <- cbind(a, age, c, height, e)
base::colnames(PrelimDatasim) <- c("UHR-ID", "Age", "Sex", "Height", "Region")
base::structure(base::as.data.frame(PrelimDatasim))
}}
