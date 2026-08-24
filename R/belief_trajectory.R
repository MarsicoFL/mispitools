#' Compute a belief trajectory from sequential evidence
#'
#' Given a prior distribution over hypotheses and a sequence of likelihood
#' ratio vectors (one per evidence step), computes the full Bayesian belief
#' trajectory by sequential multiplicative updating.
#'
#' @param prior Numeric vector. Prior distribution over the \eqn{n} hypotheses
#'   (must sum to 1, entries in \eqn{[0, 1]}).
#' @param lr_list List of numeric vectors. Each element is the likelihood
#'   ratio vector for an evidence step, with one entry per hypothesis
#'   (same length as \code{prior}).
#'
#' @return A matrix with \eqn{(T+1)} rows and \eqn{n} columns, where row 1 is
#'   the prior and row \eqn{t+1} is the posterior distribution after the
#'   first \eqn{t} evidence steps. Each row is a probability distribution.
#'
#' @details
#' The update rule is the classical multiplicative Bayesian update:
#' \deqn{P_t(H_i) \propto P_{t-1}(H_i) \cdot \mathrm{LR}_t(H_i)}
#' with renormalization after each step. In the two-hypothesis case
#' (\eqn{n = 2}), the ratio of the components of
#' \code{lr_list[[t]]} recovers the classical forensic likelihood ratio.
#'
#' @examples
#' # Two hypotheses, three evidence steps
#' prior <- c(0.5, 0.5)
#' lr_list <- list(c(3, 1), c(1, 2), c(2, 1))
#' traj <- belief_trajectory(prior, lr_list)
#' print(traj)
#'
#' @seealso \code{\link{binary_belief_trajectory}} for the two-hypothesis
#'   special case, \code{\link{trajectory_metrics}} for summary statistics.
#'
#' @references
#' Marsico, F. L. & Egeland, T. (in preparation).
#' Belief dynamics during the investigative process.
#'
#' @export
belief_trajectory <- function(prior, lr_list) {
  if (!is.numeric(prior))
    stop("prior must be a numeric vector")
  if (abs(sum(prior) - 1) > 1e-10)
    stop("prior must sum to 1")
  if (!is.list(lr_list))
    stop("lr_list must be a list of numeric vectors")
  n <- length(prior)
  TT <- length(lr_list)
  traj <- matrix(0, nrow = TT + 1, ncol = n)
  traj[1, ] <- prior
  for (t in seq_len(TT)) {
    lrt <- lr_list[[t]]
    if (length(lrt) != n)
      stop(sprintf("LR vector at step %d has length %d, expected %d",
                   t, length(lrt), n))
    if (any(lrt < 0))
      stop(sprintf("LR vector at step %d contains negative values", t))
    unnorm <- traj[t, ] * lrt
    s <- sum(unnorm)
    if (s <= 0)
      stop(sprintf("Posterior is degenerate at step %d (sum is zero)", t))
    traj[t + 1, ] <- unnorm / s
  }
  traj
}


#' Compute a binary belief trajectory from per-marker likelihood ratios
#'
#' For the two-hypothesis forensic case, computes the trajectory of
#' \eqn{P(H_1)} as markers are added sequentially.
#'
#' @param per_marker_lrs Named numeric vector. Per-marker likelihood ratios
#'   \eqn{\mathrm{LR}_k = P(\mathrm{profile}_k | H_1) / P(\mathrm{profile}_k | H_2)}.
#' @param prior_odds Numeric scalar. Prior odds for \eqn{H_1} vs \eqn{H_2}
#'   (default 1, i.e., equal prior).
#'
#' @return A data frame with one row per step (including the prior at step 0)
#'   and the following columns:
#'   \describe{
#'     \item{\code{step}}{Integer, 0 for prior, 1 to \eqn{K} for markers.}
#'     \item{\code{marker}}{Character, "Prior" or the marker name.}
#'     \item{\code{log10_lr}}{Per-step \eqn{\log_{10} \mathrm{LR}_k}.}
#'     \item{\code{cum_log10_lr}}{Cumulative \eqn{\log_{10} \mathrm{LR}} up to
#'       and including step \code{step}.}
#'     \item{\code{posterior_h1}}{\eqn{P(H_1)} at step \code{step}.}
#'     \item{\code{posterior_h2}}{\eqn{P(H_2)} at step \code{step}.}
#'   }
#'
#' @examples
#' lrs <- c(D3S1358 = 5.2, TH01 = 1.8, D21S11 = 12.0, D18S51 = 3.1)
#' binary_belief_trajectory(lrs)
#'
#' @seealso \code{\link{belief_trajectory}} for the general multi-hypothesis
#'   case, \code{\link{concentration_index}} for fragility diagnostics.
#'
#' @references
#' Marsico, F. L. & Egeland, T. (in preparation).
#' Belief dynamics during the investigative process.
#'
#' @export
binary_belief_trajectory <- function(per_marker_lrs, prior_odds = 1) {
  if (!is.numeric(per_marker_lrs))
    stop("per_marker_lrs must be a numeric vector")
  if (any(per_marker_lrs <= 0))
    stop("per_marker_lrs must be positive")
  if (!is.numeric(prior_odds) || length(prior_odds) != 1 || prior_odds <= 0)
    stop("prior_odds must be a positive scalar")

  M <- length(per_marker_lrs)
  marker_names <- names(per_marker_lrs)
  if (is.null(marker_names)) marker_names <- paste0("M", seq_len(M))

  log_lrs <- log10(per_marker_lrs)
  cum_log_lr <- cumsum(log_lrs)
  cum_lr <- 10 ^ cum_log_lr

  p_h1 <- c(prior_odds / (1 + prior_odds),
            (prior_odds * cum_lr) / (1 + prior_odds * cum_lr))

  data.frame(
    step = 0:M,
    marker = c("Prior", marker_names),
    log10_lr = c(0, log_lrs),
    cum_log10_lr = c(0, cum_log_lr),
    posterior_h1 = p_h1,
    posterior_h2 = 1 - p_h1,
    stringsAsFactors = FALSE
  )
}


#' Concentration index of per-step evidence contributions
#'
#' Computes the fraction of the total absolute evidence weight carried by
#' the single most impactful step. Ranges from \eqn{1/T} (uniform
#' contributions) to \eqn{1} (single-step dominance). Under the
#' single-failure-mode axiomatization (Theorem T-001 in the supplementary
#' material of Marsico & Egeland, in prep.), this is the canonical
#' concentration measure for forensic fragility: it equals the fraction of
#' the total weight of evidence lost under worst-case removal of a single
#' piece of evidence.
#'
#' @param weights Numeric vector. Per-step contributions, e.g.,
#'   \eqn{\log_{10}(\mathrm{LR}_k)} for the two-hypothesis forensic case or
#'   per-step Kullback-Leibler divergences for the general case.
#'
#' @return Numeric scalar in \eqn{[1/T, 1]} where \eqn{T} is the length of
#'   \code{weights}. Returns 0 if all weights are zero.
#'
#' @details
#' The formula is
#' \deqn{C_W(w) = \frac{\max_k |w_k|}{\sum_k |w_k|}}
#' with absolute values taken to handle both supporting (\eqn{w_k > 0}) and
#' excluding (\eqn{w_k < 0}) evidence symmetrically. For inclusion-fragility
#' specifically (restricting to positive weights), use
#' \code{\link{concentration_index_positive}}. For robustness checks, compare
#' against \code{\link{herfindahl_index}} and \code{\link{shannon_concentration}}.
#'
#' @examples
#' # Balanced: all markers contribute equally
#' concentration_index(rep(0.7, 15))  # -> 1/15
#'
#' # Concentrated: one marker dominates
#' concentration_index(c(5, rep(0.1, 14)))  # -> ~0.78
#'
#' @seealso
#' \code{\link{concentration_index_positive}},
#' \code{\link{herfindahl_index}},
#' \code{\link{shannon_concentration}},
#' \code{\link{leave_one_out}}.
#'
#' @references
#' Marsico, F. L. & Egeland, T. (in preparation).
#' Belief dynamics during the investigative process.
#' Slooten, K. (2021). The analogy between DNA kinship and DNA mixture
#' evaluation. \emph{Forensic Science International: Genetics} 51, 102444.
#'
#' @export
concentration_index <- function(weights) {
  if (!is.numeric(weights))
    stop("weights must be a numeric vector")
  aw <- abs(weights)
  s <- sum(aw)
  if (s == 0) return(0)
  max(aw) / s
}


#' Inclusion-fragility concentration index (positive weights only)
#'
#' Signed variant of \code{\link{concentration_index}} that uses only the
#' supporting (\eqn{w_k > 0}) contributions. For forensic inclusion cases,
#' this measures the fragility of the \emph{support} for the prosecution
#' hypothesis specifically, ignoring any excluding markers.
#'
#' @param weights Numeric vector of per-step contributions
#'   (typically \eqn{\log_{10} \mathrm{LR}_k}).
#'
#' @return Numeric scalar in \eqn{[0, 1]}. Returns 0 if no weight is positive.
#'
#' @details
#' Defined as
#' \deqn{C_W^+(w) = \frac{\max_k w_k^+}{\sum_k w_k^+}}
#' with \eqn{w_k^+ = \max(w_k, 0)}. Under Theorem T-003.1 of the Belief
#' Dynamics framework, \eqn{C_W^+} equals the fraction of the total
#' positive weight of evidence lost under adversarial removal of the single
#' most impactful supporting marker.
#'
#' @examples
#' # Mostly supporting, one dominant
#' concentration_index_positive(c(2.5, 0.3, 0.4, 0.2, 0.1))
#'
#' # Mixed with one strong exclusion: the exclusion does NOT contribute
#' concentration_index_positive(c(1.0, 0.8, 0.9, -5.0, 1.1))
#'
#' @seealso \code{\link{concentration_index}},
#' \code{\link{leave_one_out}}.
#'
#' @references
#' Marsico, F. L. & Egeland, T. (in preparation).
#' Belief dynamics during the investigative process.
#'
#' @export
concentration_index_positive <- function(weights) {
  if (!is.numeric(weights))
    stop("weights must be a numeric vector")
  wp <- pmax(weights, 0)
  s <- sum(wp)
  if (s == 0) return(0)
  max(wp) / s
}


#' Herfindahl-Hirschman concentration of evidence contributions
#'
#' Computes the Herfindahl-Hirschman index on the normalized per-step
#' contributions. A robustness check / alternative concentration measure to
#' \code{\link{concentration_index}}. Ranges from \eqn{1/T} (uniform) to
#' \eqn{1} (single-step dominance).
#'
#' @param weights Numeric vector of per-step contributions.
#'
#' @return Numeric scalar in \eqn{[1/T, 1]}. Returns 0 if all weights
#'   are zero.
#'
#' @details
#' \deqn{H(w) = \sum_{k=1}^T p_k^2, \qquad p_k = \frac{|w_k|}{\sum_j |w_j|}}
#' Unlike \code{\link{concentration_index}} (which depends only on the
#' maximum and the sum), Herfindahl depends on the full distribution of
#' weights and is therefore more sensitive to "second-tier" contributors.
#'
#' @examples
#' herfindahl_index(rep(1, 10))           # 0.1
#' herfindahl_index(c(10, rep(0.1, 9)))  # ~0.92
#'
#' @seealso \code{\link{concentration_index}},
#' \code{\link{shannon_concentration}}.
#'
#' @export
herfindahl_index <- function(weights) {
  if (!is.numeric(weights))
    stop("weights must be a numeric vector")
  aw <- abs(weights)
  s <- sum(aw)
  if (s == 0) return(0)
  p <- aw / s
  sum(p ^ 2)
}


#' Shannon-entropy-based concentration of evidence contributions
#'
#' Computes a concentration measure based on the Shannon entropy of the
#' normalized per-step contributions. Robustness alternative to
#' \code{\link{concentration_index}} and \code{\link{herfindahl_index}}.
#'
#' @param weights Numeric vector of per-step contributions.
#'
#' @return Numeric scalar in \eqn{[0, 1]}. Returns 0 for uniform weights
#'   and approaches 1 for single-step dominance. Returns 0 if all weights
#'   are zero.
#'
#' @details
#' Defined as
#' \deqn{C_S(w) = 1 - \frac{H(p)}{\log T}, \qquad p_k = \frac{|w_k|}{\sum_j |w_j|},}
#' where \eqn{H(p) = -\sum_k p_k \log p_k} is the natural-log Shannon
#' entropy of the normalized contributions. This is the complement of the
#' normalized Shannon entropy: uniform weights give \eqn{C_S = 0},
#' single-step dominance gives \eqn{C_S \to 1}.
#'
#' @examples
#' shannon_concentration(rep(1, 10))           # 0
#' shannon_concentration(c(10, rep(0.01, 9)))  # ~0.95
#'
#' @seealso \code{\link{concentration_index}},
#' \code{\link{herfindahl_index}}.
#'
#' @export
shannon_concentration <- function(weights) {
  if (!is.numeric(weights))
    stop("weights must be a numeric vector")
  aw <- abs(weights)
  s <- sum(aw)
  if (s == 0) return(0)
  p <- aw / s
  p_pos <- p[p > 0]
  if (length(p_pos) <= 1) return(1)
  H <- -sum(p_pos * log(p_pos))
  Hmax <- log(length(p))
  1 - H / Hmax
}


#' Leave-one-out fragility analysis of evidence contributions
#'
#' For each step, computes the total weight of evidence with that step
#' removed. Identifies cases where a single step dominates: if
#' \code{fraction[k]} is close to 1, the total weight hinges on that single
#' step and the conclusion would be substantially weakened if it were
#' excluded or challenged.
#'
#' @param per_marker_lrs Named numeric vector. Per-marker likelihood ratios
#'   (positive values).
#'
#' @return A data frame with columns
#'   \describe{
#'     \item{\code{marker}}{Marker name.}
#'     \item{\code{log10_lr}}{\eqn{\log_{10} \mathrm{LR}_k}, this marker's
#'       contribution to the total.}
#'     \item{\code{total_without}}{Total \eqn{\log_{10} \mathrm{LR}} with
#'       this marker removed.}
#'     \item{\code{fraction}}{\eqn{|\log_{10} \mathrm{LR}_k| /
#'       \sum_j |\log_{10} \mathrm{LR}_j|}, this marker's share of the
#'       absolute total.}
#'   }
#'
#' @examples
#' lrs <- c(D3S1358 = 5.2, TH01 = 1.8, D21S11 = 120.0, D18S51 = 3.1)
#' leave_one_out(lrs)
#'
#' @seealso \code{\link{concentration_index}},
#' \code{\link{concentration_index_positive}}.
#'
#' @references
#' Marsico, F. L. & Egeland, T. (in preparation).
#' Belief dynamics during the investigative process.
#'
#' @export
leave_one_out <- function(per_marker_lrs) {
  if (!is.numeric(per_marker_lrs))
    stop("per_marker_lrs must be a numeric vector")
  if (any(per_marker_lrs <= 0))
    stop("per_marker_lrs must be positive")

  M <- length(per_marker_lrs)
  marker_names <- names(per_marker_lrs)
  if (is.null(marker_names)) marker_names <- paste0("M", seq_len(M))

  log_lrs <- log10(per_marker_lrs)
  total <- sum(log_lrs)
  abs_total <- sum(abs(log_lrs))

  data.frame(
    marker = marker_names,
    log10_lr = log_lrs,
    total_without = total - log_lrs,
    fraction = if (abs_total > 0) abs(log_lrs) / abs_total else rep(0, M),
    stringsAsFactors = FALSE
  )
}


#' Shannon entropy in base-10 (bans)
#'
#' Computes the Shannon entropy of a probability distribution, using
#' base-10 logarithm so the result is expressed in \emph{bans}
#' (the forensic convention for weight-of-evidence units).
#'
#' @param p Numeric vector. Probability distribution (must sum to 1, entries
#'   in \eqn{[0, 1]}).
#' @param epsilon Numeric scalar. Small positive constant to avoid
#'   \eqn{\log(0)} for zero-probability outcomes. Default \eqn{10^{-20}}.
#'
#' @return Numeric scalar. Shannon entropy \eqn{H(p) = -\sum_i p_i
#'   \log_{10} p_i} in bans.
#'
#' @examples
#' entropy_log10(c(0.5, 0.5))          # log10(2)
#' entropy_log10(rep(1/12, 12))        # log10(12)
#' entropy_log10(c(1, 0, 0, 0))        # ~0
#'
#' @seealso \code{\link{kl_divergence_log10}},
#' \code{\link{trajectory_metrics}}.
#'
#' @export
entropy_log10 <- function(p, epsilon = 1e-20) {
  if (!is.numeric(p))
    stop("p must be a numeric vector")
  p <- pmax(p, epsilon)
  -sum(p * log10(p))
}


#' Kullback-Leibler divergence in base-10 (bans)
#'
#' Computes the Kullback-Leibler divergence \eqn{D_{\mathrm{KL}}(p \| q)}
#' between two distributions, in base-10 logarithm so the result is
#' expressed in \emph{bans}.
#'
#' @param p Numeric vector. Target distribution.
#' @param q Numeric vector. Reference distribution (same length as \code{p}).
#' @param epsilon Numeric scalar. Small positive constant to avoid
#'   \eqn{\log(0)}. Default \eqn{10^{-20}}.
#'
#' @return Numeric scalar. \eqn{D_{\mathrm{KL}}(p \| q) = \sum_i p_i
#'   \log_{10}(p_i / q_i)} in bans.
#'
#' @examples
#' p <- c(0.7, 0.3)
#' q <- c(0.5, 0.5)
#' kl_divergence_log10(p, q)          # positive
#' kl_divergence_log10(p, p)          # 0
#'
#' @seealso \code{\link{entropy_log10}},
#' \code{\link{trajectory_metrics}}.
#'
#' @export
kl_divergence_log10 <- function(p, q, epsilon = 1e-20) {
  if (!is.numeric(p) || !is.numeric(q))
    stop("p and q must be numeric vectors")
  if (length(p) != length(q))
    stop("p and q must have the same length")
  p <- pmax(p, epsilon)
  q <- pmax(q, epsilon)
  sum(p * log10(p / q))
}


#' Trajectory metrics from a belief trajectory matrix
#'
#' Computes all trajectory-level metrics for a belief trajectory: entropy
#' at each step, per-step Kullback-Leibler divergence (Bayesian surprise),
#' cumulative divergence from the prior, per-step total-variation distance,
#' total path length, and a family of concentration indices (max/sum,
#' Herfindahl, Shannon-based).
#'
#' @param traj_matrix Matrix. Output from \code{\link{belief_trajectory}}
#'   with rows indexing time steps (row 1 = prior, row \eqn{T+1} = final
#'   posterior) and columns indexing hypotheses.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{entropy}}{Numeric vector of length \eqn{T+1}: Shannon
#'       entropy at each step.}
#'     \item{\code{kl_step}}{Numeric vector of length \eqn{T}: per-step
#'       Kullback-Leibler divergence \eqn{D_{\mathrm{KL}}(P_t \| P_{t-1})},
#'       equivalently Bayesian surprise (Itti & Baldi, 2009).}
#'     \item{\code{kl_from_prior}}{Numeric vector of length \eqn{T+1}:
#'       cumulative Kullback-Leibler divergence from the prior,
#'       \eqn{D_{\mathrm{KL}}(P_t \| P_0)}.}
#'     \item{\code{tv_step}}{Numeric vector of length \eqn{T}: per-step
#'       total-variation distance \eqn{\mathrm{TV}(P_t, P_{t-1})}.}
#'     \item{\code{path_length}}{Numeric scalar. Total path length in
#'       total-variation distance, \eqn{\sum_t \mathrm{TV}(P_t, P_{t-1})}.}
#'     \item{\code{concentration}}{Numeric scalar. Max/sum concentration
#'       index based on per-step information gains (see
#'       \code{\link{concentration_index}}).}
#'     \item{\code{concentration_herfindahl}}{Numeric scalar. Herfindahl
#'       concentration (see \code{\link{herfindahl_index}}).}
#'     \item{\code{concentration_shannon}}{Numeric scalar. Shannon-based
#'       concentration (see \code{\link{shannon_concentration}}).}
#'   }
#'
#' @examples
#' prior <- c(0.5, 0.5)
#' lrs <- list(c(3, 1), c(1, 2), c(2, 1))
#' traj <- belief_trajectory(prior, lrs)
#' metrics <- trajectory_metrics(traj)
#' str(metrics)
#'
#' @seealso \code{\link{belief_trajectory}},
#' \code{\link{entropy_log10}},
#' \code{\link{kl_divergence_log10}},
#' \code{\link{concentration_index}}.
#'
#' @references
#' Marsico, F. L. & Egeland, T. (in preparation).
#' Belief dynamics during the investigative process.
#' Itti, L. & Baldi, P. (2009). Bayesian surprise attracts human attention.
#' \emph{Vision Research} 49, 1295-1306.
#'
#' @export
trajectory_metrics <- function(traj_matrix) {
  if (!is.matrix(traj_matrix) || !is.numeric(traj_matrix))
    stop("traj_matrix must be a numeric matrix")
  if (nrow(traj_matrix) < 2)
    stop("traj_matrix must have at least 2 rows (prior + one step)")

  TT <- nrow(traj_matrix) - 1
  ent <- apply(traj_matrix, 1, entropy_log10)

  kl_s <- vapply(seq_len(TT), function(t)
    kl_divergence_log10(traj_matrix[t + 1, ], traj_matrix[t, ]),
    numeric(1))

  kl_p <- vapply(seq_len(nrow(traj_matrix)), function(t)
    kl_divergence_log10(traj_matrix[t, ], traj_matrix[1, ]),
    numeric(1))

  tv_s <- vapply(seq_len(TT), function(t)
    0.5 * sum(abs(traj_matrix[t + 1, ] - traj_matrix[t, ])),
    numeric(1))

  list(
    entropy = ent,
    kl_step = kl_s,
    kl_from_prior = kl_p,
    tv_step = tv_s,
    path_length = sum(tv_s),
    concentration = concentration_index(kl_s),
    concentration_herfindahl = herfindahl_index(kl_s),
    concentration_shannon = shannon_concentration(kl_s)
  )
}
