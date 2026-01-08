#' Likelihood Ratio for Age Variable
#'
#' @description
#' Simulates age observations and optionally computes likelihood ratios (LRs)
#' under either H1 (unidentified person is the missing person) or H2
#' (unidentified person is not the missing person).
#'
#' Ages are categorized into two groups based on whether they fall within
#' the missing person's estimated age range:
#' \itemize{
#'   \item T1: Age within range (MPa - MPr) to (MPa + MPr)
#'   \item T0: Age outside range
#' }
#'
#' @param MPa Numeric. Missing person's estimated age in years. Default: 40.
#' @param MPr Numeric. Age range tolerance (plus/minus years). The matching
#'   interval is (MPa - MPr) to (MPa + MPr). Default: 6.
#' @param UHRr Numeric. Additional uncertainty range for the unidentified
#'   person's age estimation. Default: 1.
#' @param gam Numeric. Gamma parameter for age uncertainty scaling.
#'   The uncertainty interval is age +/- (gam * age + UHRr). Default: 0.07.
#' @param nsims Integer. Number of simulations to perform. Default: 1000.
#' @param epa Numeric (0-1). Error rate for age categorization. Default: 0.05.
#' @param erRa Numeric (0-1). Error rate in the reference/database.
#'   Defaults to \code{epa}.
#' @param H Integer (1 or 2). Hypothesis to simulate under:
#'   \itemize{
#'     \item 1: H1 (Related) - Unidentified person IS the missing person
#'     \item 2: H2 (Unrelated) - Unidentified person is NOT the missing person
#'   }
#'   Default: 1.
#' @param modelA Character. Reference age distribution model:
#'   \itemize{
#'     \item "uniform": Assumes uniform age distribution (default)
#'     \item "custom": Uses empirical frequencies from simulations
#'   }
#' @param LR Logical. If TRUE, compute and return LR values. Default: FALSE.
#' @param seed Integer. Random seed for reproducibility. Default: 1234.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{group}: Age group classification ("T1" or "T0")
#'     \item \code{age}: Simulated age value
#'     \item \code{UHRmin}: Lower bound of uncertainty interval
#'     \item \code{UHRmax}: Upper bound of uncertainty interval
#'     \item \code{LRa}: Likelihood ratio (only if LR = TRUE)
#'   }
#'
#' @details
#' \strong{Under H1 (Related):}
#' Age is sampled to fall within the MP's range with probability (1 - erRa),
#' outside with probability erRa.
#'
#' \strong{Under H2 (Unrelated):}
#' Age is sampled uniformly from 1-80, then categorized.
#'
#' \strong{LR Calculation (uniform model):}
#' \itemize{
#'   \item LR(T1) = (1 - epa) / P(T1), where P(T1) = 2*MPr/80
#'   \item LR(T0) = epa / P(T0), where P(T0) = 1 - P(T1)
#' }
#'
#' @seealso
#' \code{\link{sim_lr_prelim}} for unified preliminary LR simulations,
#' \code{\link{lr_sex}}, \code{\link{lr_hair_color}} for other variables.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @import dplyr
#' @examples
#' # Simulate under H1 (related)
#' sim_h1 <- lr_age(MPa = 40, MPr = 6, H = 1, nsims = 100)
#' table(sim_h1$group)
#'
#' # Simulate under H2 with LR values
#' sim_h2 <- lr_age(MPa = 40, MPr = 6, H = 2, nsims = 100, LR = TRUE)
#' head(sim_h2)
#'
#' # Narrower age range (more discriminating)
#' sim_narrow <- lr_age(MPa = 35, MPr = 3, nsims = 500, LR = TRUE)
#' summary(sim_narrow$LRa)

lr_age <- function(MPa = 40,
                   MPr = 6,
                   UHRr = 1,
                   gam = 0.07,
                   nsims = 1000,
                   epa = 0.05,
                   erRa = epa,
                   H = 1,
                   modelA = c("uniform", "custom")[1],
                   LR = FALSE,
                   seed = 1234) {

  set.seed(seed)
  sims <- list()
  Age <- seq(1, 80)
  MPmin <- MPa - MPr
  MPmax <- MPa + MPr

  # Calculate probabilities under uniform model
  if (modelA == "uniform") {
    T1p <- (MPmax - MPmin) / length(Age)
    T0p <- 1 - T1p
    LR1 <- (1 - epa) / T1p
    LR0 <- epa / T0p
  }

  # Define age groups
  T1a <- Age[Age < MPmax & Age > MPmin]
  T0a <- Age[-T1a]

  if (H == 1) {
    # H1: Sample with high probability of being in range
    group <- unlist(sample(c("T1", "T0"), size = nsims, prob = c(1 - erRa, erRa), replace = TRUE))
    ages <- unlist(lapply(group, function(x) ifelse(x == "T1", sample(T1a, 1), sample(T0a, 1))))

    sims <- as.data.frame(cbind(group, ages))
    names(sims) <- c("group", "age")

    sims <- dplyr::mutate(sims, UHRmin = as.numeric(ages) - gam * as.numeric(ages) - UHRr)
    sims <- dplyr::mutate(sims, UHRmax = as.numeric(ages) + gam * as.numeric(ages) + UHRr)
  }
  else if (H == 2) {
    # H2: Sample uniformly from population
    ages <- unlist(sample(Age, nsims, replace = TRUE))
    group <- unlist(lapply(ages, function(x) ifelse(x > MPmin & x < MPmax, "T1", "T0")))

    sims <- as.data.frame(cbind(group, ages))
    names(sims) <- c("group", "age")
    sims <- dplyr::mutate(sims, UHRmin = as.numeric(ages) - gam * as.numeric(ages) - UHRr)
    sims <- dplyr::mutate(sims, UHRmax = as.numeric(ages) + gam * as.numeric(ages) + UHRr)
  }

  # Custom model: calculate probabilities from empirical data
  if (modelA == "custom") {
    T1p <- length(subset(sims$group, sims$group == "T1")) / length(sims$group)
    T0p <- 1 - T1p
    LR1 <- (1 - epa) / T1p
    LR0 <- epa / T0p
  }

  if (LR == TRUE) {
    sims <- dplyr::mutate(sims, LRa = ifelse(group == "T1", LR1, LR0))
    names(sims) <- c("group", "Age", "UHRmin", "UHRmax", "LRa")
    return(sims)
  }

  return(as.data.frame(sims))
}
