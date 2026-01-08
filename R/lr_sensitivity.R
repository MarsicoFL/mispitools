#' Sensitivity Analysis for Likelihood Ratios
#'
#' @description
#' Evaluates how the likelihood ratio changes when model parameters vary.
#' This is essential for understanding the robustness of forensic conclusions
#' and for communicating uncertainty to decision-makers.
#'
#' @param evidence_type Character. Type of evidence to analyze.
#'   Options: "sex", "age", "hair", "region".
#' @param param Character. Parameter to vary. Options depend on evidence_type:
#'   \itemize{
#'     \item \code{"sex"}: "eps" (error rate), "freq" (population frequency)
#'     \item \code{"age"}: "eps" (error rate), "range" (age interval)
#'     \item \code{"hair"}: "eps" (error rate), "freq" (population frequency)
#'     \item \code{"region"}: "eps" (error rate), "nreg" (number of regions)
#'   }
#' @param range Numeric vector of length 2. Range of parameter values to test.
#'   Default depends on param type.
#' @param steps Integer. Number of steps in the range. Default: 20.
#' @param match Logical. TRUE for matching evidence (same sex/age in range/etc),
#'   FALSE for mismatching. Default: TRUE.
#' @param baseline List. Baseline parameter values. If NULL, uses defaults.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{param_value}: Parameter value tested
#'     \item \code{LR}: Resulting likelihood ratio
#'     \item \code{log10_LR}: Log10 of LR (useful for plotting)
#'   }
#'
#' @details
#' Sensitivity analysis is critical in forensic science because:
#' \enumerate{
#'   \item Parameters (error rates, population frequencies) are often estimated
#'     with uncertainty
#'   \item Different reference populations may have different frequencies
#'   \item The analysis reveals which parameters most affect conclusions
#' }
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item Steep curves indicate high sensitivity (conclusions depend strongly
#'     on parameter choice)
#'   \item Flat curves indicate robustness (conclusions stable across
#'     reasonable parameter values)
#' }
#'
#' @seealso
#' \code{\link{lr_sex}}, \code{\link{lr_age}}, \code{\link{lr_hair_color}}
#' for individual LR calculations.
#'
#' @references
#' Kling D, Tillmar AO, Egeland T (2014). "Familias 3-Extensions and new
#' functionality." \emph{Forensic Science International: Genetics}, 13, 121-127.
#'
#' @export
#' @examples
#' # How does sex LR change with error rate?
#' sens_eps <- lr_sensitivity("sex", param = "eps", range = c(0.01, 0.20))
#' plot(sens_eps$param_value, sens_eps$log10_LR, type = "l",
#'      xlab = "Error rate", ylab = "log10(LR)",
#'      main = "Sex LR sensitivity to error rate")
#' abline(h = 0, lty = 2)
#'
#' # How does sex LR change with population frequency?
#' sens_freq <- lr_sensitivity("sex", param = "freq", range = c(0.3, 0.7))
#' plot(sens_freq$param_value, sens_freq$log10_LR, type = "l",
#'      xlab = "Female frequency", ylab = "log10(LR)",
#'      main = "Sex LR sensitivity to population frequency")
#'
#' # Age LR sensitivity to range parameter
#' sens_range <- lr_sensitivity("age", param = "range", range = c(2, 15))
#' plot(sens_range$param_value, sens_range$log10_LR, type = "l",
#'      xlab = "Age range (+/- years)", ylab = "log10(LR)",
#'      main = "Age LR sensitivity to range")

lr_sensitivity <- function(evidence_type,
                           param,
                           range = NULL,
                           steps = 20,
                           match = TRUE,
                           baseline = NULL) {

  # Validate evidence type
  valid_types <- c("sex", "age", "hair", "region")
  if (!evidence_type %in% valid_types) {
    stop("evidence_type must be one of: ", paste(valid_types, collapse = ", "))
  }

  # Validate param for each evidence type
  valid_params <- list(
    sex = c("eps", "freq"),
    age = c("eps", "range"),
    hair = c("eps", "freq"),
    region = c("eps", "nreg")
  )

  if (!param %in% valid_params[[evidence_type]]) {
    stop("For evidence_type '", evidence_type, "', param must be one of: ",
         paste(valid_params[[evidence_type]], collapse = ", "))
  }

  # Set default ranges
  if (is.null(range)) {
    range <- switch(param,
      eps = c(0.01, 0.20),
      freq = c(0.30, 0.70),
      range = c(2, 15),
      nreg = c(2, 10)
    )
  }

  # Generate parameter sequence
  param_values <- seq(range[1], range[2], length.out = steps)

  # Calculate LR for each parameter value
  LRs <- numeric(steps)

  for (i in seq_along(param_values)) {
    val <- param_values[i]

    LRs[i] <- switch(evidence_type,
      sex = calc_lr_sex_sensitivity(param, val, match, baseline),
      age = calc_lr_age_sensitivity(param, val, match, baseline),
      hair = calc_lr_hair_sensitivity(param, val, match, baseline),
      region = calc_lr_region_sensitivity(param, val, match, baseline)
    )
  }

  data.frame(
    param_value = param_values,
    LR = LRs,
    log10_LR = log10(LRs)
  )
}

# Internal helper functions for sensitivity calculations

calc_lr_sex_sensitivity <- function(param, val, match, baseline) {
  # Defaults
  eps <- if (!is.null(baseline$eps)) baseline$eps else 0.05
  freq <- if (!is.null(baseline$freq)) baseline$freq else 0.5

  if (param == "eps") eps <- val
  if (param == "freq") freq <- val

  # LR calculation
  # P(match | H1) = 1 - eps
  # P(mismatch | H1) = eps
  # P(match | H2) = freq (frequency of MP's sex in population)
  # P(mismatch | H2) = 1 - freq

  if (match) {
    LR <- (1 - eps) / freq
  } else {
    LR <- eps / (1 - freq)
  }
  LR
}

calc_lr_age_sensitivity <- function(param, val, match, baseline) {
  # Defaults
  eps <- if (!is.null(baseline$eps)) baseline$eps else 0.05
  age_range <- if (!is.null(baseline$range)) baseline$range else 6
  max_age <- 80

  if (param == "eps") eps <- val
  if (param == "range") age_range <- val

  # P(in range | H2) under uniform distribution
  p_range_h2 <- (2 * age_range) / max_age

  if (match) {
    LR <- (1 - eps) / p_range_h2
  } else {
    LR <- eps / (1 - p_range_h2)
  }
  LR
}

calc_lr_hair_sensitivity <- function(param, val, match, baseline) {
  # Defaults
  eps <- if (!is.null(baseline$eps)) baseline$eps else 0.05
  freq <- if (!is.null(baseline$freq)) baseline$freq else 0.3  # black hair

  if (param == "eps") eps <- val
  if (param == "freq") freq <- val

  # Simplified model: match/mismatch with single error rate
  if (match) {
    LR <- (1 - eps) / freq
  } else {
    LR <- eps / (1 - freq)
  }
  LR
}

calc_lr_region_sensitivity <- function(param, val, match, baseline) {
  # Defaults
  eps <- if (!is.null(baseline$eps)) baseline$eps else 0.05
  nreg <- if (!is.null(baseline$nreg)) baseline$nreg else 6

  if (param == "eps") eps <- val
  if (param == "nreg") nreg <- round(val)

  # P(match | H2) = 1/nreg (uniform)
  p_match_h2 <- 1 / nreg
  p_mismatch_h2 <- 1 - p_match_h2

  if (match) {
    LR <- (1 - eps) / p_match_h2
  } else {
    LR <- eps / p_mismatch_h2
  }
  LR
}
