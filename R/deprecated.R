#' Deprecated functions in mispitools
#'
#' @description
#' These functions are deprecated and will be removed in a future version
#' of mispitools. Please use the new function names instead.
#'
#' The function names have been updated to follow a consistent snake_case
#' naming convention with semantic prefixes:
#' \itemize{
#'   \item \code{sim_*}: Simulation functions
#'   \item \code{lr_*}: Likelihood ratio functions
#'   \item \code{plot_*}: Visualization functions
#'   \item \code{cpt_*}: Conditional probability table functions
#'   \item \code{kl_*}: Kullback-Leibler divergence functions
#' }
#'
#' @name mispitools-deprecated
#' @keywords internal
NULL

# ==============================================================================
# Simulation functions
# ==============================================================================

#' @rdname mispitools-deprecated
#' @export
simLRgen <- function(...) {
  .Deprecated("sim_lr_genetic", package = "mispitools",
              msg = "simLRgen() is deprecated. Use sim_lr_genetic() instead.")
  sim_lr_genetic(...)
}

#' @rdname mispitools-deprecated
#' @export
simLRprelim <- function(...) {
  .Deprecated("sim_lr_prelim", package = "mispitools",
              msg = "simLRprelim() is deprecated. Use sim_lr_prelim() instead.")
  sim_lr_prelim(...)
}

#' @rdname mispitools-deprecated
#' @export
simLR2dataframe <- function(...) {
  .Deprecated("lr_to_dataframe", package = "mispitools",
              msg = "simLR2dataframe() is deprecated. Use lr_to_dataframe() instead.")
  lr_to_dataframe(...)
}

#' @rdname mispitools-deprecated
#' @export
simRef <- function(...) {
  .Deprecated("sim_reference_pop", package = "mispitools",
              msg = "simRef() is deprecated. Use sim_reference_pop() instead.")
  sim_reference_pop(...)
}

#' @rdname mispitools-deprecated
#' @export
makePOIgen <- function(...) {
  .Deprecated("sim_poi_genetic", package = "mispitools",
              msg = "makePOIgen() is deprecated. Use sim_poi_genetic() instead.")
  sim_poi_genetic(...)
}

#' @rdname mispitools-deprecated
#' @export
makePOIprelim <- function(...) {
  .Deprecated("sim_poi_prelim", package = "mispitools",
              msg = "makePOIprelim() is deprecated. Use sim_poi_prelim() instead.")
  sim_poi_prelim(...)
}

#' @rdname mispitools-deprecated
#' @export
makeMPprelim <- function(...) {
  .Deprecated("sim_mp_prelim", package = "mispitools",
              msg = "makeMPprelim() is deprecated. Use sim_mp_prelim() instead.")
  sim_mp_prelim(...)
}

#' @rdname mispitools-deprecated
#' @export
postSim <- function(...) {
  .Deprecated("sim_posterior", package = "mispitools",
              msg = "postSim() is deprecated. Use sim_posterior() instead.")
  sim_posterior(...)
}

# ==============================================================================
# LR calculation functions
# ==============================================================================

#' @rdname mispitools-deprecated
#' @export
LRsex <- function(...) {
  .Deprecated("lr_sex", package = "mispitools",
              msg = "LRsex() is deprecated. Use lr_sex() instead.")
  lr_sex(...)
}

#' @rdname mispitools-deprecated
#' @export
LRage <- function(...) {
  .Deprecated("lr_age", package = "mispitools",
              msg = "LRage() is deprecated. Use lr_age() instead.")
  lr_age(...)
}

#' @rdname mispitools-deprecated
#' @export
LRcol <- function(...) {
  .Deprecated("lr_hair_color", package = "mispitools",
              msg = "LRcol() is deprecated. Use lr_hair_color() instead.")
  lr_hair_color(...)
}

#' @rdname mispitools-deprecated
#' @export
LRcolors <- function(...) {
  .Deprecated("lr_pigmentation", package = "mispitools",
              msg = "LRcolors() is deprecated. Use lr_pigmentation() instead.")
  lr_pigmentation(...)
}

#' @rdname mispitools-deprecated
#' @export
LRdate <- function(...) {
  .Deprecated("lr_birthdate", package = "mispitools",
              msg = "LRdate() is deprecated. Use lr_birthdate() instead.")
  lr_birthdate(...)
}

#' @rdname mispitools-deprecated
#' @export
combLR <- function(...) {
  .Deprecated("lr_combine", package = "mispitools",
              msg = "combLR() is deprecated. Use lr_combine() instead.")
  lr_combine(...)
}

#' @rdname mispitools-deprecated
#' @export
compute_LRs_colors <- function(...) {
  .Deprecated("lr_compute_pigmentation", package = "mispitools",
              msg = "compute_LRs_colors() is deprecated. Use lr_compute_pigmentation() instead.")
  lr_compute_pigmentation(...)
}

# ==============================================================================
# CPT functions
# ==============================================================================

#' @rdname mispitools-deprecated
#' @export
CPT_POP <- function(...) {
  .Deprecated("cpt_population", package = "mispitools",
              msg = "CPT_POP() is deprecated. Use cpt_population() instead.")
  cpt_population(...)
}

#' @rdname mispitools-deprecated
#' @export
CPT_MP <- function(...) {
  .Deprecated("cpt_missing_person", package = "mispitools",
              msg = "CPT_MP() is deprecated. Use cpt_missing_person() instead.")
  cpt_missing_person(...)
}

#' @rdname mispitools-deprecated
#' @export
Cmodel <- function(...) {
  .Deprecated("error_matrix_hair", package = "mispitools",
              msg = "Cmodel() is deprecated. Use error_matrix_hair() instead.")
  error_matrix_hair(...)
}

# ==============================================================================
# Visualization functions
# ==============================================================================

#' @rdname mispitools-deprecated
#' @export
LRdist <- function(...) {
  .Deprecated("plot_lr_distribution", package = "mispitools",
              msg = "LRdist() is deprecated. Use plot_lr_distribution() instead.")
  plot_lr_distribution(...)
}

#' @rdname mispitools-deprecated
#' @export
deplot <- function(...) {
  .Deprecated("plot_decision_curve", package = "mispitools",
              msg = "deplot() is deprecated. Use plot_decision_curve() instead.")
  plot_decision_curve(...)
}

#' @rdname mispitools-deprecated
#' @export
CondPlot <- function(...) {
  .Deprecated("plot_cpt", package = "mispitools",
              msg = "CondPlot() is deprecated. Use plot_cpt() instead.")
  plot_cpt(...)
}

# ==============================================================================
# Decision functions
# ==============================================================================

#' @rdname mispitools-deprecated
#' @export
DeT <- function(...) {
  .Deprecated("decision_threshold", package = "mispitools",
              msg = "DeT() is deprecated. Use decision_threshold() instead.")
  decision_threshold(...)
}

#' @rdname mispitools-deprecated
#' @export
Trates <- function(...) {
  .Deprecated("threshold_rates", package = "mispitools",
              msg = "Trates() is deprecated. Use threshold_rates() instead.")
  threshold_rates(...)
}

# ==============================================================================
# Utility functions
# ==============================================================================

#' @rdname mispitools-deprecated
#' @export
getfreqs <- function(...) {
  .Deprecated("get_allele_freqs", package = "mispitools",
              msg = "getfreqs() is deprecated. Use get_allele_freqs() instead.")
  get_allele_freqs(...)
}

#' @rdname mispitools-deprecated
#' @export
conditionedProp <- function(...) {
  .Deprecated("compute_conditioned_prop", package = "mispitools",
              msg = "conditionedProp() is deprecated. Use compute_conditioned_prop() instead.")
  compute_conditioned_prop(...)
}

#' @rdname mispitools-deprecated
#' @export
refProp <- function(...) {
  .Deprecated("compute_reference_prop", package = "mispitools",
              msg = "refProp() is deprecated. Use compute_reference_prop() instead.")
  compute_reference_prop(...)
}

# ==============================================================================
# KL divergence functions
# ==============================================================================

#' @rdname mispitools-deprecated
#' @export
bidirectionalKL <- function(...) {
  .Deprecated("kl_bidirectional", package = "mispitools",
              msg = "bidirectionalKL() is deprecated. Use kl_bidirectional() instead.")
  kl_bidirectional(...)
}

#' @rdname mispitools-deprecated
#' @export
klPIE <- function(...) {
  .Deprecated("kl_pie", package = "mispitools",
              msg = "klPIE() is deprecated. Use kl_pie() instead.")
  kl_pie(...)
}

#' @rdname mispitools-deprecated
#' @export
multi_kl_divergence <- function(...) {
  .Deprecated("kl_multi", package = "mispitools",
              msg = "multi_kl_divergence() is deprecated. Use kl_multi() instead.")
  kl_multi(...)
}

# ==============================================================================
# Shiny apps
# ==============================================================================

#' @rdname mispitools-deprecated
#' @export
mispiApp <- function(...) {
  .Deprecated("mispitools_app", package = "mispitools",
              msg = "mispiApp() is deprecated. Use mispitools_app() instead.")
  mispitools_app(...)
}

#' @rdname mispitools-deprecated
#' @export
lrComparisonApp <- function(...) {
  .Deprecated("mispitools_app", package = "mispitools",
              msg = "lrComparisonApp() is deprecated. Use mispitools_app() instead.")
  mispitools_app(...)
}

#' @rdname mispitools-deprecated
#' @export
app_mispitools <- function(...) {
  .Deprecated("mispitools_app", package = "mispitools",
              msg = "app_mispitools() is deprecated. Use mispitools_app() instead.")
  mispitools_app(...)
}

#' @rdname mispitools-deprecated
#' @export
app_lr_comparison <- function(...) {
  .Deprecated("mispitools_app", package = "mispitools",
              msg = "app_lr_comparison() is deprecated. Use mispitools_app() instead.")
  mispitools_app(...)
}
