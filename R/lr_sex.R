#' Likelihood Ratio for Biological Sex
#'
#' @description
#' Simulates observations of biological sex and optionally computes likelihood
#' ratios (LRs) under either H1 (unidentified person is the missing person) or
#' H2 (unidentified person is not the missing person).
#'
#' @param MPs Character. Missing person's biological sex: "F" for female,
#'   "M" for male. Default: "F".
#' @param eps Numeric (0-1). Error rate (epsilon) for sex observation.
#'   Probability of misclassifying sex when recording. Default: 0.05.
#' @param erRs Numeric (0-1). Error rate in the database/reference.
#'   Defaults to \code{eps} if not specified.
#' @param nsims Integer. Number of simulations to perform. Default: 1000.
#' @param Ps Numeric vector of length 2. Sex proportions in the population,
#'   c(proportion_female, proportion_male). Must sum to 1. Default: c(0.5, 0.5).
#' @param H Integer (1 or 2). Hypothesis to simulate under:
#'   \itemize{
#'     \item 1: H1 (Related) - Unidentified person IS the missing person
#'     \item 2: H2 (Unrelated) - Unidentified person is NOT the missing person
#'   }
#'   Default: 1.
#' @param LR Logical. If TRUE, compute and return LR values for each
#'   simulated observation. Default: FALSE.
#' @param seed Integer. Random seed for reproducibility. Default: 1234.
#'
#' @return A data.frame with column \code{Sexo} containing simulated sex
#'   observations ("F" or "M"). If \code{LR = TRUE}, also includes column
#'   \code{LRs} with the likelihood ratio for each observation.
#'
#' @details
#' \strong{Under H1 (Related):}
#' The observed sex matches the MP's true sex with probability (1 - erRs),
#' and is incorrectly recorded with probability erRs.
#'
#' \strong{Under H2 (Unrelated):}
#' Sex is sampled from the population proportions Ps.
#'
#' \strong{LR Calculation:}
#' For a matching observation: \code{LR = (1 - eps) / Ps_MP}
#' For a non-matching observation: \code{LR = eps / Ps_other}
#'
#' @seealso
#' \code{\link{sim_lr_prelim}} for unified preliminary LR simulations,
#' \code{\link{lr_age}}, \code{\link{lr_hair_color}} for other variables.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Simulate under H1 (related)
#' sim_h1 <- lr_sex(MPs = "F", H = 1, nsims = 100)
#' table(sim_h1$Sexo)
#'
#' # Simulate under H2 (unrelated) with LR values
#' sim_h2 <- lr_sex(MPs = "F", H = 2, nsims = 100, LR = TRUE)
#' head(sim_h2)
#'
#' # Different population proportions
#' sim_custom <- lr_sex(
#'   MPs = "M",
#'   Ps = c(0.52, 0.48),  # 52% female population
#'   nsims = 500,
#'   LR = TRUE
#' )
#' summary(sim_custom$LRs)

lr_sex <- function(MPs = "F",
                   eps = 0.05,
                   erRs = eps,
                   nsims = 1000,
                   Ps = c(0.5, 0.5),
                   H = 1,
                   LR = FALSE,
                   seed = 1234) {

  set.seed(seed)
  sims <- list()
  S <- c("F", "M")

  # Identify MP's sex position and opposite
  MPss <- which(S == MPs)
  NoMPsn <- S[-MPss]
  noMPs <- which(S == NoMPsn)

  if (H == 1) {
    # H1: Sample with high probability of matching MP's sex
    x <- c(S[MPss], S[noMPs])
    sims <- as.data.frame(sample(x, size = nsims, prob = c(1 - erRs, erRs), replace = TRUE))
    names(sims) <- "Sexo"
  }
  else if (H == 2) {
    # H2: Sample from population proportions
    x <- c(S[MPss], S[noMPs])
    sims <- as.data.frame(sample(x, size = nsims, prob = c(Ps[1], Ps[2]), replace = TRUE))
    names(sims) <- "Sexo"
  }

  if (LR == TRUE) {
    # Calculate LR for each observation
    LRmatch <- (1 - eps) / Ps[MPss]
    LRnomatch <- eps / Ps[noMPs]
    LRs <- lapply(sims, function(x) ifelse(x == MPs, LRmatch, LRnomatch))
    sims <- cbind(sims, LRs)
    names(sims) <- c("Sexo", "LRs")
    return(sims)
  }

  return(sims)
}
