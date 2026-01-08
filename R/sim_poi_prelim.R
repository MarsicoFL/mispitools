#' Simulate Preliminary Investigation Data for Persons of Interest
#'
#' @description
#' Generates a simulated database of preliminary investigation data for
#' persons of interest (POIs) or unidentified human remains (UHRs). Supports
#' two case types: missing children and missing migrants.
#'
#' @param casetype Character. Type of missing person case:
#'   \itemize{
#'     \item "children": Generates birth date, sex, birth type, and region
#'     \item "migrants": Generates age, sex, height, and region
#'   }
#'   Default: "children".
#' @param dateinit Character. Minimum birth date for simulated POIs in
#'   "YYYY/MM/DD" format. Only used for casetype = "children".
#'   Default: "1975/01/01".
#' @param scenario Integer (1 or 2). Birth date distribution scenario:
#'   \itemize{
#'     \item 1: Non-uniform (gamma distribution, more realistic)
#'     \item 2: Uniform distribution
#'   }
#'   Only used for casetype = "children". Default: 1.
#' @param femaleprop Numeric (0-1). Proportion of females in the simulated
#'   population. Default: 0.5.
#' @param ext Numeric. Extension parameter:
#'   \itemize{
#'     \item Scenario 1: Scale factor for gamma distribution
#'     \item Scenario 2: Number of days range
#'   }
#'   Default: 100.
#' @param numsims Integer. Number of POIs/UHRs to simulate. Default: 10000.
#' @param seed Integer. Random seed for reproducibility. Default: 123.
#' @param birthprob Numeric vector of length 3. Probabilities for birth type:
#'   c(home_birth, hospital_birth, unknown/adoption). Only for "children".
#'   Default: c(0.09, 0.9, 0.01).
#' @param region Character vector. Names of regions/locations. Default:
#'   c("North America", "South America", "Africa", "Asia", "Europe", "Oceania").
#' @param regionprob Numeric vector. Probabilities for each region.
#'   Must sum to 1 and have same length as \code{region}.
#'   Default: c(0.2, 0.2, 0.2, 0.1, 0.2, 0.1).
#'
#' @return A data.frame with columns depending on casetype:
#'   \itemize{
#'     \item \strong{children}: POI-ID, DBD (declared birth date), Sex,
#'       Birth-type, Region
#'     \item \strong{migrants}: UHR-ID, Age, Sex, Height, Region
#'   }
#'
#' @details
#' For missing children cases, this function simulates characteristics of
#' children who may have been taken during periods of conflict or human rights
#' violations, with their identity documents potentially falsified.
#'
#' For missing migrants cases, this simulates characteristics of unidentified
#' human remains that may correspond to missing migrants.
#'
#' The birth date distribution in scenario 1 uses a gamma distribution with
#' shape=12, which creates a more realistic non-uniform pattern of births.
#'
#' @seealso
#' \code{\link{sim_mp_prelim}} for simulating missing person data,
#' \code{\link{sim_lr_prelim}} for using this data in LR calculations.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @importFrom stats rgamma
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @examples
#' # Simulate children case database
#' db_children <- sim_poi_prelim(
#'   casetype = "children",
#'   dateinit = "1975/01/01",
#'   scenario = 1,
#'   numsims = 100,
#'   seed = 123
#' )
#' head(db_children)
#'
#' # Simulate migrants case database
#' db_migrants <- sim_poi_prelim(
#'   casetype = "migrants",
#'   numsims = 100,
#'   seed = 456
#' )
#' head(db_migrants)
#' summary(db_migrants$Age)

sim_poi_prelim <- function(casetype = "children",
                           dateinit = "1975/01/01",
                           scenario = 1,
                           femaleprop = 0.5,
                           ext = 100,
                           numsims = 10000,
                           seed = 123,
                           birthprob = c(0.09, 0.9, 0.01),
                           region = c("North America", "South America", "Africa",
                                      "Asia", "Europe", "Oceania"),
                           regionprob = c(0.2, 0.2, 0.2, 0.1, 0.2, 0.1)) {

  if (casetype == "children") {
    sex <- c("female", "male")
    maleprop <- 1 - femaleprop
    birth <- c("home birth", "hospital birth", "unknown-adoption")
    set.seed(seed)

    a <- seq(1, numsims, by = 1)

    if (scenario == 1) {
      b <- as.data.frame(as.Date(dateinit) + ext * stats::rgamma(numsims, 12))
    } else if (scenario == 2) {
      b <- as.data.frame(as.Date(dateinit) + stats::runif(numsims, min = 0, max = ext))
    }

    c <- sample(sex, numsims, replace = TRUE, prob = c(femaleprop, maleprop))
    d <- sample(birth, numsims, replace = TRUE, prob = birthprob)
    e <- sample(region, numsims, replace = TRUE, prob = regionprob)

    PrelimDatasim <- cbind(a, b, c, d, e)
    base::colnames(PrelimDatasim) <- c("POI-ID", "DBD", "Sex", "Birth-type", "Region")
    base::structure(base::as.data.frame(PrelimDatasim))

  } else if (casetype == "migrants") {
    sex <- c("female", "male")
    maleprop <- 1 - femaleprop
    set.seed(seed)

    age <- sample(seq.int(from = 4, to = 70, by = 1), size = numsims, replace = TRUE)
    height <- stats::rnorm(numsims, mean = 170, sd = 15)

    a <- seq(1, numsims, by = 1)
    c <- sample(sex, numsims, replace = TRUE, prob = c(femaleprop, maleprop))
    e <- sample(region, numsims, replace = TRUE, prob = regionprob)

    PrelimDatasim <- cbind(a, age, c, height, e)
    base::colnames(PrelimDatasim) <- c("UHR-ID", "Age", "Sex", "Height", "Region")
    base::structure(base::as.data.frame(PrelimDatasim))
  }
}
