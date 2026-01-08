#' Kullback-Leibler Divergence for Probability Matrices
#'
#' @description
#' Computes the Kullback-Leibler (KL) divergence between two probability
#' matrices using base-10 logarithm. Calculates divergence in both directions.
#'
#' This is useful for comparing conditional probability tables (CPTs) or
#' other matrix representations of probability distributions.
#'
#' @param P A numeric matrix representing the first probability distribution.
#'   Should sum to 1 (or will be treated as unnormalized).
#' @param Q A numeric matrix representing the second probability distribution.
#'   Must have the same dimensions as P.
#' @param min_value Numeric. Minimum value to replace zeros (to avoid
#'   undefined logarithms). Default: 1e-12.
#'
#' @return A named numeric vector with two elements:
#'   \itemize{
#'     \item \code{"P || Q"}: KL divergence from P to Q
#'     \item \code{"Q || P"}: KL divergence from Q to P
#'   }
#'
#' @details
#' The KL divergence is computed as:
#' \deqn{KL(P || Q) = \sum_i P_i \log_{10}(P_i / Q_i)}
#'
#' Zero values in P or Q are replaced with \code{min_value} to avoid
#' undefined operations.
#'
#' @seealso
#' \code{\link{kl_bidirectional}} for allele frequency comparisons,
#' \code{\link{cpt_population}}, \code{\link{cpt_missing_person}} for
#' creating probability matrices.
#'
#' @references
#' Kullback S, Leibler RA (1951). "On Information and Sufficiency."
#' \emph{The Annals of Mathematical Statistics}, 22(1), 79-86.
#'
#' @export
#' @examples
#' # Compare two CPTs
#' cpt1 <- cpt_population()
#' cpt2 <- cpt_population(propS = c(0.6, 0.4))
#'
#' kl_pie(cpt1, cpt2)

kl_pie <- function(P, Q, min_value = 1e-12) {

  if (!is.numeric(P) || !is.numeric(Q)) {
    stop("P and Q must be numeric matrices.")
  }

  if (any(dim(P) != dim(Q))) {
    stop("P and Q must have the same dimensions.")
  }

  # Replace zeros with minimum value
  P <- pmax(P, min_value)
  Q <- pmax(Q, min_value)

  # Calculate KL divergence in both directions
  kl_pq <- sum(P * log10(P / Q))
  kl_qp <- sum(Q * log10(Q / P))

  return(c("P || Q" = kl_pq, "Q || P" = kl_qp))
}
