#' Simulate Posterior Odds Combining Genetic and Non-Genetic Evidence
#'
#' @description
#' Simulates posterior odds distributions by combining prior probabilities
#' with likelihood ratios from both genetic and non-genetic (preliminary
#' investigation) evidence. This implements a full Bayesian integration
#' of multiple evidence sources.
#'
#' @param datasim Output from \code{\link{sim_lr_genetic}} containing
#'   genetic LR simulations.
#' @param Prior Numeric (0-1). Prior probability for H1 (that POI is MP).
#'   Default: 0.01.
#' @param PriorModel Character. How to incorporate preliminary evidence:
#'   \itemize{
#'     \item "prelim": Combines prior with preliminary data LRs (sex, age, color)
#'     \item "uniform": Uses only the prior probability without preliminary LRs
#'   }
#'   Default: "prelim".
#' @param eps Numeric (0-1). Error rate for sex observation. Default: 0.05.
#' @param erRs Numeric (0-1). Error rate for sex in database. Default: 0.01.
#' @param epc Matrix. Hair color error matrix from \code{\link{error_matrix_hair}}.
#' @param erRc Matrix. Hair color error matrix for database.
#' @param MPc Integer (1-5). Missing person's hair color. Default: 1.
#' @param epa Numeric (0-1). Error rate for age. Default: 0.05.
#' @param erRa Numeric (0-1). Error rate for age in database. Default: 0.01.
#' @param MPa Numeric. Missing person's age. Default: 10.
#' @param MPr Numeric. Age range tolerance. Default: 2.
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{Unrelated}: Posterior odds under H2
#'     \item \code{Related}: Posterior odds under H1
#'   }
#'
#' @details
#' Posterior odds are calculated as:
#' \deqn{Posterior = Prior \times LR_{prelim} \times LR_{genetic}}
#'
#' Where:
#' \itemize{
#'   \item Prior = P(H1) / P(H2) = Prior / (1 - Prior)
#'   \item LR_prelim = LR_sex * LR_age * LR_color
#'   \item LR_genetic = from genetic simulation
#' }
#'
#' @seealso
#' \code{\link{sim_lr_genetic}} for genetic simulations,
#' \code{\link{lr_sex}}, \code{\link{lr_age}}, \code{\link{lr_hair_color}}
#' for preliminary evidence LRs.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' library(forrel)
#'
#' # Setup pedigree
#' x <- linearPed(2)
#' x <- setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x <- profileSim(x, N = 1, ids = 2)
#'
#' # Simulate genetic LRs
#' datasim <- sim_lr_genetic(x, missing = 5, numsims = 50, seed = 123)
#'
#' # Compute posterior odds with preliminary evidence
#' post <- sim_posterior(datasim, Prior = 0.01, PriorModel = "prelim")
#' head(post)
#'
#' # Visualize
#' plot_lr_distribution(post)

sim_posterior <- function(datasim,
                          Prior = 0.01,
                          PriorModel = c("prelim", "uniform")[1],
                          eps = 0.05,
                          erRs = 0.01,
                          epc = error_matrix_hair(),
                          erRc = error_matrix_hair(),
                          MPc = 1,
                          epa = 0.05,
                          erRa = 0.01,
                          MPa = 10,
                          MPr = 2) {

  # Extract LR values from genetic simulations
  unrelated_values <- vector()
  related_values <- vector()
  list_length <- length(datasim[["Unrelated"]])

  for (i in 1:list_length) {
    unrelated_value <- datasim[["Unrelated"]][[i]][["LRtotal"]][["H1:H2"]]
    related_value <- datasim[["Related"]][[i]][["LRtotal"]][["H1:H2"]]

    unrelated_values <- c(unrelated_values, unrelated_value)
    related_values <- c(related_values, related_value)
  }

  results_df <- data.frame(Unrelated = unrelated_values, Related = related_values)
  datasim <- results_df
  nsims <- nrow(datasim)

  if (PriorModel == "prelim") {
    LRgenH2 <- datasim$Unrelated

    # Simulations for H2
    LRs2 <- lr_sex(LR = TRUE, H = 2, nsims = nsims, eps = eps, erRs = erRs)
    LRc2 <- lr_hair_color(LR = TRUE, H = 2, nsims = nsims, epc = epc, erRc = erRc, MPc = MPc)
    LRa2 <- lr_age(LR = TRUE, H = 2, nsims = nsims, epa = epa, erRa = erRa, MPa = MPa, MPr = MPr)

    datasim2 <- data.frame(
      LRs = LRs2$LRs,
      LRc = LRc2$LRc,
      LRa = LRa2$LRa
    )

    datasim2$LRtot <- datasim2$LRa * datasim2$LRc * datasim2$LRs
    datasim2$Label <- "H2"
    datasim2$LRgenH2 <- LRgenH2
    datasim2$PriorOdds <- Prior / (1 - Prior)
    datasim2$PostOdds <- datasim2$PriorOdds * datasim2$LRtot
    datasim2$PostOddsGen <- datasim2$PostOdds * datasim2$LRgenH2

    # Simulations for H1
    LRgenH1 <- datasim$Related

    LRs1 <- lr_sex(LR = TRUE, H = 1, nsims = nsims, eps = eps, erRs = erRs)
    LRc1 <- lr_hair_color(LR = TRUE, H = 1, nsims = nsims, epc = epc, erRc = erRc, MPc = MPc)
    LRa1 <- lr_age(LR = TRUE, H = 1, nsims = nsims, epa = epa, erRa = erRa, MPa = MPa, MPr = MPr)

    datasim1 <- data.frame(
      LRs = LRs1$LRs,
      LRc = LRc1$LRc,
      LRa = LRa1$LRa
    )

    datasim1$LRtot <- datasim1$LRa * datasim1$LRc * datasim1$LRs
    datasim1$Label <- "H1"
    datasim1$LRgenH1 <- LRgenH1
    datasim1$PriorOdds <- Prior / (1 - Prior)
    datasim1$PostOdds <- datasim1$PriorOdds * datasim1$LRtot
    datasim1$PostOddsGen <- datasim1$PostOdds * datasim1$LRgenH1

  } else {
    LRgenH2 <- datasim$Unrelated
    datasim2 <- data.frame(Unrelated = LRgenH2)
    datasim2$PriorOdds <- Prior / (1 - Prior)
    datasim2$PostOdds <- datasim2$PriorOdds * datasim2$Unrelated

    LRgenH1 <- datasim$Related
    datasim1 <- data.frame(Related = LRgenH1)
    datasim1$PriorOdds <- Prior / (1 - Prior)
    datasim1$PostOdds <- datasim1$PriorOdds * datasim1$Related
  }

  data <- data.frame(
    Unrelated = datasim2$PostOddsGen,
    Related = datasim1$PostOddsGen
  )

  return(data)
}
