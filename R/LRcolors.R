#' Simulate LR distributions for pigmentation traits under H1 and H2
#'
#' This function samples LR values with probabilities proportional to the numerator
#' (related, H1) and the denominator (unrelated, H2) across hair/skin/eye colour combinations.
#' Use `compute_LRs_colors()` to prepare the input data with matching trait columns.
#'
#' @param df A data.frame with `numerators`, `f_h_s_y`, and `LR` columns, typically the output
#'   of `compute_LRs_colors()`.
#' @param seed For replication purposes.
#' @param nsim Number of LRs simulated per hypothesis.
#' @return A data.frame with columns `Unrelated` (H2) and `Related` (H1).
#' @export
#' @examples
#' data <- simRef()
#' conditioned <- conditionedProp(data, 1, 1, 1, 0.01, 0.01, 0.01)
#' unconditioned <- refProp(data)
#' lrs <- compute_LRs_colors(conditioned, unconditioned)
#' LRcolors(lrs, nsim = 100)

LRcolors <- function(df, seed = 1234, nsim = 500) {
  set.seed(seed)
  
  LR <- df$LR
  f_h_s_y <- df$f_h_s_y
  numerators <- df$numerators
  
  prob_unrelated <- f_h_s_y / sum(f_h_s_y)
  prob_related <- numerators / sum(numerators)
  
  Unrelated <- sample(LR, nsim, replace = TRUE, prob = prob_unrelated)
  Related <- sample(LR, nsim, replace = TRUE, prob = prob_related)
  
  result_df <- data.frame(Unrelated = Unrelated, Related = Related)
  
  return(result_df)
}
