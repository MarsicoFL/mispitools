#' Simulate Preliminary Investigation Data for Missing Persons
#'
#' @description
#' Generates a simulated database of preliminary investigation data for
#' missing persons (MPs). This complements \code{\link{sim_poi_prelim}} which
#' generates data for persons of interest. Supports two case types: missing
#' children and missing migrants.
#'
#' @param casetype Character. Type of missing person case:
#'   \itemize{
#'     \item "children": Generates birth date, sex, birth month, and birth place
#'     \item "migrants": Generates age, sex, height, and region
#'   }
#'   Default: "children".
#' @param dateinit Character. Minimum birth date for simulated MPs in
#'   "YYYY/MM/DD" format. Only for casetype = "children". Default: "1975/01/01".
#' @param scenario Integer (1 or 2). Birth date distribution scenario:
#'   \itemize{
#'     \item 1: Non-uniform (gamma distribution)
#'     \item 2: Uniform distribution
#'   }
#'   Only for casetype = "children". Default: 1.
#' @param femaleprop Numeric (0-1). Proportion of females. Default: 0.5.
#' @param ext Numeric. Extension parameter for date simulation. Default: 100.
#' @param numsims Integer. Number of MPs to simulate. Default: 10000.
#' @param seed Integer. Random seed for reproducibility. Default: 123.
#' @param region Character vector. Names of regions/locations.
#'   Default: c("North America", "South America", "Africa", "Asia", "Europe", "Oceania").
#' @param regionprob Numeric vector. Probabilities for each region.
#'   Default: c(0.2, 0.2, 0.2, 0.1, 0.2, 0.1).
#'
#' @return A data.frame with columns depending on casetype:
#'   \itemize{
#'     \item \strong{children}: POI-ID, DBD, Sex, Month, Birth place
#'     \item \strong{migrants}: UHR-ID, Age, Sex, Height, Region
#'   }
#'
#' @details
#' This function generates the "ground truth" characteristics of missing
#' persons, while \code{\link{sim_poi_prelim}} generates the observed/recorded
#' characteristics of persons of interest (which may include observation
#' errors or falsified data).
#'
#' @seealso
#' \code{\link{sim_poi_prelim}} for simulating POI data,
#' \code{\link{sim_lr_prelim}} for using this data in LR calculations.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @importFrom stats rgamma runif rnorm
#' @examples
#' # Simulate missing children data
#' mp_children <- sim_mp_prelim(casetype = "children", numsims = 100, seed = 123)
#' head(mp_children)
#'
#' # Simulate missing migrants data
#' mp_migrants <- sim_mp_prelim(casetype = "migrants", numsims = 100, seed = 456)
#' head(mp_migrants)

sim_mp_prelim <- function(casetype = "children",
                          dateinit = "1975/01/01",
                          scenario = 1,
                          femaleprop = 0.5,
                          ext = 100,
                          numsims = 10000,
                          seed = 123,
                          region = c("North America", "South America", "Africa",
                                     "Asia", "Europe", "Oceania"),
                          regionprob = c(0.2, 0.2, 0.2, 0.1, 0.2, 0.1)) {

  if (casetype == "children") {
    sex <- c("female", "male")
    maleprop <- 1 - femaleprop
    set.seed(seed)

    a <- seq(1, numsims, by = 1)

    if (scenario == 1) {
      b <- as.data.frame(as.Date(dateinit) + ext * stats::rgamma(numsims, 12))
    } else if (scenario == 2) {
      b <- as.data.frame(as.Date(dateinit) + stats::runif(numsims, min = 0, max = ext))
    }

    c <- sample(sex, numsims, replace = TRUE, prob = c(femaleprop, maleprop))
    d <- sample(seq.int(from = 1, to = 9, by = 1), size = numsims, replace = TRUE)
    e <- sample(region, numsims, replace = TRUE, prob = regionprob)

    PrelimDatasim <- cbind(a, b, c, d, e)
    colnames(PrelimDatasim) <- c("POI-ID", "DBD", "Sex", "Month", "Birth place")
    as.data.frame(PrelimDatasim)

  } else if (casetype == "migrants") {
    sex <- c("female", "male")
    maleprop <- 1 - femaleprop
    set.seed(seed)

    age <- sample(seq.int(from = 1, to = 9, by = 1), size = numsims, replace = TRUE)
    height <- stats::rnorm(numsims, mean = 170, sd = 15)

    a <- seq(1, numsims, by = 1)
    c <- sample(sex, numsims, replace = TRUE, prob = c(femaleprop, maleprop))
    e <- sample(region, numsims, replace = TRUE, prob = regionprob)

    PrelimDatasim <- cbind(a, age, c, height, e)
    colnames(PrelimDatasim) <- c("UHR-ID", "Age", "Sex", "Height", "Region")
    as.data.frame(PrelimDatasim)
  }
}
