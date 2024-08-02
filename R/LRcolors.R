#' Simulate LR values considering H1 and H2
#'
#'
#' @param df A data.frame containing the characteristics of individuals, numerator, f_h_s_y and LRs. Output from compute_LRs function.
#' @param seed For replication purposes.
#' @param nsim Number of LRs simulated.
#' @return LR distribution considering H1 (Related) and H2 (Unrelated).
#' 
#' @export

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
