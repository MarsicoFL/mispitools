#' Likelihood Ratio for Hair Color
#'
#' @description
#' Simulates hair color observations and optionally computes likelihood ratios
#' (LRs) under either H1 (unidentified person is the missing person) or H2
#' (unidentified person is not the missing person).
#'
#' Hair color is categorized into 5 groups:
#' 1=Black, 2=Brown, 3=Blonde, 4=Red, 5=Gray/White
#'
#' @param MPc Integer (1-5). Missing person's hair color category. Default: 1.
#' @param epc Matrix. Hair color error/confusion matrix, typically created
#'   with \code{\link{error_matrix_hair}}. Rows represent true colors,
#'   columns represent observed colors. Default: \code{error_matrix_hair()}.
#' @param erRc Matrix. Error matrix for the reference/database.
#'   Defaults to \code{epc}.
#' @param numsims Integer. Number of simulations to perform. Default: 1000.
#' @param nsims Deprecated. Use \code{numsims} instead.
#' @param Pc Numeric vector of length 5. Hair color proportions in the
#'   population. Must sum to 1. Default: c(0.3, 0.2, 0.25, 0.15, 0.1).
#' @param H Integer (1 or 2). Hypothesis to simulate under:
#'   \itemize{
#'     \item 1: H1 (Related) - Unidentified person IS the missing person
#'     \item 2: H2 (Unrelated) - Unidentified person is NOT the missing person
#'   }
#'   Default: 1.
#' @param Qprop Integer. Query color for testing. Defaults to \code{MPc}.
#' @param LR Logical. If TRUE, compute and return LR values. Default: FALSE.
#' @param seed Integer. Random seed for reproducibility. Default: 1234.
#'
#' @return A data.frame with column \code{Col} containing simulated color
#'   observations (1-5). If \code{LR = TRUE}, also includes column \code{LRc}
#'   with the likelihood ratio for each observation.
#'
#' @details
#' \strong{Under H1 (Related):}
#' Observed color is sampled using the row of the error matrix corresponding
#' to the MP's true hair color. This accounts for observation errors.
#'
#' \strong{Under H2 (Unrelated):}
#' Color is sampled from the population proportions Pc.
#'
#' \strong{LR Calculation:}
#' LR = P(observed color | true color is MPc) / P(observed color in population)
#' \code{LR = epc(MPc, observed) / Pc(observed)}
#'
#' @seealso
#' \code{\link{error_matrix_hair}} for creating the error matrix,
#' \code{\link{lr_pigmentation}} for combined pigmentation traits,
#' \code{\link{sim_lr_prelim}} for unified preliminary LR simulations.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Simulate under H1 (related) with black hair MP
#' sim_h1 <- lr_hair_color(MPc = 1, H = 1, numsims = 100)
#' table(sim_h1$Col)
#'
#' # Simulate under H2 with LR values
#' sim_h2 <- lr_hair_color(MPc = 2, H = 2, numsims = 100, LR = TRUE)
#' head(sim_h2)
#' summary(sim_h2$LRc)
#'
#' # Custom population proportions
#' sim_custom <- lr_hair_color(
#'   MPc = 3,  # Blonde
#'   Pc = c(0.1, 0.4, 0.3, 0.1, 0.1),  # Different population
#'   numsims = 500,
#'   LR = TRUE
#' )

lr_hair_color <- function(MPc = 1,
                          epc = error_matrix_hair(),
                          erRc = epc,
                          numsims = 1000,
                          Pc = c(0.3, 0.2, 0.25, 0.15, 0.1),
                          H = 1,
                          Qprop = MPc,
                          LR = FALSE,
                          seed = 1234,
                          nsims = NULL) {

  # Handle deprecated nsims parameter
  if (!is.null(nsims)) {
    warning("Parameter 'nsims' is deprecated. Use 'numsims' instead.",
            call. = FALSE)
    numsims <- nsims
  }

  # Input validation
  if (!is.numeric(MPc) || length(MPc) != 1 || MPc < 1 || MPc > 5 || MPc != round(MPc)) {
    stop("MPc must be an integer between 1 and 5")
  }
  if (!is.matrix(epc) || nrow(epc) != 5 || ncol(epc) != 5) {
    stop("epc must be a 5x5 matrix")
  }
  if (!is.matrix(erRc) || nrow(erRc) != 5 || ncol(erRc) != 5) {
    stop("erRc must be a 5x5 matrix")
  }
  if (!is.numeric(numsims) || length(numsims) != 1 || numsims < 1) {
    stop("numsims must be a positive integer")
  }
  if (!is.numeric(Pc) || length(Pc) != 5) {
    stop("Pc must be a numeric vector of length 5")
  }
  if (any(Pc < 0)) {
    stop("Pc values must be non-negative")
  }
  if (abs(sum(Pc) - 1) > 1e-6) {
    warning("Pc does not sum to 1; normalizing")
    Pc <- Pc / sum(Pc)
  }
  if (any(Pc == 0)) {
    warning("Pc contains zero values; LR may be Inf for some observations")
  }
  if (!H %in% c(1, 2)) {
    stop("H must be 1 or 2")
  }

  sims <- list()
  Col <- c(1, 2, 3, 4, 5)

  set.seed(seed)

  if (H == 1) {
    # H1: Sample using error matrix row for MP's true color
    x <- Col
    sims <- as.data.frame(sample(x, size = numsims, prob = erRc[MPc, ], replace = TRUE))
    names(sims) <- "Col"
  }
  else if (H == 2) {
    # H2: Sample from population proportions
    x <- Col
    sims <- as.data.frame(sample(x, size = numsims, prob = Pc, replace = TRUE))
    names(sims) <- "Col"
  }

  if (LR == TRUE) {
    # Calculate LR for each observation
    # LR = P(observed | true=MPc) / P(observed in population)
    LRs <- lapply(sims, function(x) epc[MPc, x] / Pc[x])
    sims <- cbind(sims, LRs)
    names(sims) <- c("Col", "LRc")
    return(sims)
  }

  return(sims)
}
