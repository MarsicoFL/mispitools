# Computations used by mispitools_app(), kept outside the app so they can
# be called and tested without starting Shiny. They used to be defined
# inside mispitools_app() itself. The bodies are unchanged from that
# version; only the names gained an app_ prefix so they don't get mixed
# up with the exported lr_*() functions, which use different models.

#' LR for sex evidence, as shown in the app
#'
#' @param MPs Sex of the missing person, `"F"` or `"M"`.
#' @param obsSex Observed sex of the unidentified person, `"F"` or `"M"`.
#' @param eps Probability that the recorded sex is wrong.
#' @param propF Proportion of females in the population.
#' @return The LR, a length-one numeric named after `obsSex`.
#' @noRd
app_lr_sex <- function(MPs, obsSex, eps, propF) {
  propS <- c(F = propF, M = 1 - propF)
  if (MPs == obsSex) {
    LR <- (1 - eps) / propS[obsSex]
  } else {
    LR <- eps / propS[obsSex]
  }
  LR
}

#' LR for age evidence, as shown in the app
#'
#' Population ages are taken as uniform on 80 years. The interval
#' `MPa +/- MPr` is clamped to `[1, 80]`, and P(in range | H2) is kept
#' inside `[0.001, 0.999]`.
#'
#' @param MPa Missing person's age.
#' @param MPr Half-width of the age interval.
#' @param obsAge Observed (estimated) age of the unidentified person.
#' @param epa Probability that the age falls outside the interval under H1.
#' @return The LR as a single number.
#' @noRd
app_lr_age <- function(MPa, MPr, obsAge, epa) {
  MPmin <- max(1, MPa - MPr)
  MPmax <- min(80, MPa + MPr)
  T1p <- (MPmax - MPmin) / 80
  T1p <- max(0.001, min(0.999, T1p))
  T0p <- 1 - T1p

  if (obsAge >= MPmin && obsAge <= MPmax) {
    LR <- (1 - epa) / T1p
  } else {
    LR <- epa / T0p
  }
  LR
}

#' LR for hair colour evidence, as shown in the app
#'
#' Uniform error model: the observed colour equals the missing person's
#' with probability `1 - epc`, otherwise it is any of the other colours
#' with equal probability.
#'
#' @param MPc Missing person's hair colour, an index into `propC`.
#' @param obsC Observed hair colour, an index into `propC`.
#' @param epc Probability that the observed colour is wrong.
#' @param propC Population frequency of each colour.
#' @return The LR, a length-one numeric (keeps the name of `propC[obsC]`).
#' @noRd
app_lr_hair <- function(MPc, obsC, epc, propC) {
  # Uniform error model
  n_colors <- length(propC)
  if (MPc == obsC) {
    p_obs_h1 <- 1 - epc
  } else {
    p_obs_h1 <- epc / (n_colors - 1)
  }
  LR <- p_obs_h1 / propC[obsC]
  LR
}
