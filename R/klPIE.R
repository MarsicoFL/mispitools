#' Calculate Kullback-Leibler Divergence with Base 10 Logarithm
#'
#' This function computes the Kullback-Leibler (KL) divergence between two 
#' probability distributions represented by matrices, using a base 10 logarithm. 
#' The function calculates KL divergence in both directions (P || Q and Q || P) 
#' and handles zero probabilities by replacing them with a minimum value to avoid 
#' undefined logarithms.
#'
#' @param P A numeric matrix representing the first probability distribution. 
#' The entire matrix should sum to 1.
#' @param Q A numeric matrix representing the second probability distribution.
#' The entire matrix should sum to 1.
#' @param min_value A numeric value representing the minimum value to replace 
#' any zero probabilities. Defaults to \code{1e-12}.
#'
#' @return A named numeric vector with two elements:
#' \describe{
#'   \item{"P || Q"}{The KL divergence from P to Q (P || Q).}
#'   \item{"Q || P"}{The KL divergence from Q to P (Q || P).}
#' }
#' 
#'
#' @export
klPIE <- function(P, Q, min_value = 1e-12) {
  
  if (!is.numeric(P) || !is.numeric(Q)) {
    stop("P and Q must be numeric matrices.")
  }
  
  if (any(dim(P) != dim(Q))) {
    stop("P and Q must have the same dimensions.")
  }
  
  
  P <- pmax(P, min_value)
  Q <- pmax(Q, min_value)
  
  kl_pq <- sum(P * log10(P / Q))
  kl_qp <- sum(Q * log10(Q / P))
  
  return(c("P || Q" = kl_pq, "Q || P" = kl_qp))
}
