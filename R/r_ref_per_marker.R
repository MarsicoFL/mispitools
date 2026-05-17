## R reference engine — per-marker KL and LR distribution (F1.6).
##
## Internal-only helpers (no @export) that consume `cpt_marker_joint_R`.
## The C++ engine (F2.x and F3.x) will be validated against these as
## oracles in F1.7 and onward. Slow on large pedigrees by design.
##
## Both quantities are derived from the joint (P_H1, P_H2) returned by
## `cpt_marker_joint_R(model, poi)`, evaluated on the joint support
## (union of supports of P_H1 and P_H2 — states with both = 0 are
## already filtered upstream).
##
## Conventions for boundary states (the joint reports both):
##   P1 > 0, P2 = 0  →  log10 LR = +Inf  →  contributes +Inf to
##                       KL(H1 || H2) and 0 to KL(H2 || H1).
##   P1 = 0, P2 > 0  →  log10 LR = -Inf  →  contributes 0 to
##                       KL(H1 || H2) and +Inf to KL(H2 || H1).
##
## Under `mutation = "none"`, Mendelian-incompatible states (e.g. AA x AA
## with an AB child) yield P_H1 = 0 with P_H2 > 0, so
## `kl_h2_to_h1 = +Inf` and `e_log10_lr_h2 = -Inf`. Under
## `mutation = "equal"` or `"stepwise"` with rate > 0, P_H1 > 0
## everywhere P_H2 > 0 and both KLs are finite.
##
## Sums are taken with the standard limit `0 * log(0/x) = 0`.

#' @noRd
per_marker_lr_dist_R <- function(model, poi = NULL, aggregate = TRUE) {
  if (!inherits(model, "marker_model")) {
    stop("`model` must be a 'marker_model' object.", call. = FALSE)
  }
  joint <- cpt_marker_joint_R(model, poi = poi)
  log10_lr <- log10_lr_from_probs(joint$P_H1, joint$P_H2)

  out <- data.frame(
    log10_lr = log10_lr,
    p_h1 = joint$P_H1,
    p_h2 = joint$P_H2
  )
  if (aggregate) {
    out <- aggregate_lr_dist(out)
  }
  attr(out, "marker_id") <- attr(joint, "marker_id")
  attr(out, "poi") <- attr(joint, "poi")
  out
}

#' @noRd
per_marker_kl_R <- function(model, poi = NULL) {
  if (!inherits(model, "marker_model")) {
    stop("`model` must be a 'marker_model' object.", call. = FALSE)
  }
  joint <- cpt_marker_joint_R(model, poi = poi)
  P1 <- joint$P_H1
  P2 <- joint$P_H2
  log10_lr <- log10_lr_from_probs(P1, P2)
  ln10 <- log(10)

  e_h1 <- weighted_log10_lr_sum(P1, log10_lr)
  e_h2 <- weighted_log10_lr_sum(P2, log10_lr)

  data.frame(
    marker = attr(joint, "marker_id"),
    e_log10_lr_h1 = e_h1,
    e_log10_lr_h2 = e_h2,
    kl_h1_to_h2 = e_h1 * ln10,
    kl_h2_to_h1 = -e_h2 * ln10,
    stringsAsFactors = FALSE
  )
}

## R reference for the F4.2 composition step. Exact sequential
## convolution of independent per-feature LR distributions: the total
## log10 LR is the sum of the per-feature log10 LRs (conditional
## independence), so the composed distribution is their convolution.
## The C++ kernel `core::lr_dist_compose` (exact mode) reproduces this
## bit-for-bit (1e-12 tol) — the cartesian product is iterated in the
## same order (feature outermost, accumulator innermost) so the
## post-sort summation in `aggregate_lr_dist` matches.

#' @noRd
lr_dist_compose_R <- function(dists) {
  acc <- data.frame(log10_lr = 0, p_h1 = 1, p_h2 = 1)
  for (d in dists) {
    n_a <- nrow(acc)
    n_d <- nrow(d)
    ia <- rep(seq_len(n_a), times = n_d)
    id <- rep(seq_len(n_d), each = n_a)
    new_lr <- acc$log10_lr[ia] + d$log10_lr[id]
    new_p1 <- acc$p_h1[ia] * d$p_h1[id]
    new_p2 <- acc$p_h2[ia] * d$p_h2[id]
    keep <- new_p1 > 0 | new_p2 > 0
    nd <- data.frame(
      log10_lr = new_lr[keep],
      p_h1 = new_p1[keep],
      p_h2 = new_p2[keep]
    )
    if (nrow(nd) == 0L) {
      return(nd)
    }
    acc <- aggregate_lr_dist(nd)
  }
  acc
}

#' @noRd
log10_lr_from_probs <- function(P1, P2) {
  out <- rep(NA_real_, length(P1))
  pos1 <- P1 > 0
  pos2 <- P2 > 0
  both <- pos1 & pos2
  out[both] <- log10(P1[both]) - log10(P2[both])
  out[pos1 & !pos2] <- Inf
  out[!pos1 & pos2] <- -Inf
  out
}

#' @noRd
weighted_log10_lr_sum <- function(weight, log10_lr) {
  active <- weight > 0
  if (!any(active)) return(0)
  parts <- weight[active] * log10_lr[active]
  sum(parts)
}

## R reference for the F4.3 decision-theoretic primitives. Every
## quantity is derived from an (log10_lr, p_h1, p_h2) sparse LR
## distribution `d` (the output of `per_marker_lr_dist_R` /
## `lr_dist_compose_R`). The C++ kernel in `src/core/decision.cpp`
## reproduces these to 1e-12: same active-weight filter, same
## ascending summation order, same tie / ±Inf conventions.

#' @noRd
weighted_centered_sq <- function(weight, log10_lr, m) {
  active <- weight > 0
  if (!any(active)) return(0)
  if (!is.finite(m)) return(Inf)
  parts <- weight[active] * (log10_lr[active] - m)^2
  sum(parts)
}

#' @noRd
lr_dist_summary_R <- function(d) {
  lr <- d$log10_lr
  p1 <- d$p_h1
  p2 <- d$p_h2
  m1 <- weighted_log10_lr_sum(p1, lr)
  m2 <- weighted_log10_lr_sum(p2, lr)
  v1 <- weighted_centered_sq(p1, lr, m1)
  v2 <- weighted_centered_sq(p2, lr, m2)
  list(
    mean_h1 = m1,
    mean_h2 = m2,
    var_h1 = v1,
    var_h2 = v2,
    sd_h1 = sqrt(v1),
    sd_h2 = sqrt(v2),
    mass_h1 = sum(p1),
    mass_h2 = sum(p2),
    has_pos_inf = any(is.infinite(lr) & lr > 0),
    has_neg_inf = any(is.infinite(lr) & lr < 0)
  )
}

## Discrete inverse-CDF, R `quantile` type 1: Q(p) = inf{x : F(x) >= p}
## over the active (weight > 0) atoms, ascending, with a 1e-12
## cumulative slack so a probability landing exactly on a breakpoint
## resolves to that atom.
#' @noRd
lr_dist_quantile_R <- function(d, probs, under_h1 = TRUE) {
  lr <- d$log10_lr
  w <- if (under_h1) d$p_h1 else d$p_h2
  keep <- w > 0
  if (!any(keep)) {
    stop("no active support under the requested hypothesis.", call. = FALSE)
  }
  lr <- lr[keep]
  w <- w[keep]
  ord <- order(lr)
  lr <- lr[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  slack <- 1e-12
  vapply(probs, function(p) {
    if (p <= 0) return(lr[1])
    hit <- which(cw + slack >= p)
    if (length(hit) == 0L) lr[length(lr)] else lr[hit[1]]
  }, numeric(1))
}

## Error rates of the `log10 LR > threshold` classifier. Mass at
## exactly == threshold is indeterminate (excluded from both tallies),
## the analytic analogue of `threshold_rates()` over the exact
## distribution.
#' @noRd
decision_rates_R <- function(d, threshold) {
  lr <- d$log10_lr
  p1 <- d$p_h1
  p2 <- d$p_h2
  below <- lr < threshold
  above <- lr > threshold
  fnr <- sum(p1[below])
  tnr <- sum(p2[below])
  tpr <- sum(p1[above])
  fpr <- sum(p2[above])
  TP <- tpr; TN <- tnr; FP <- fpr; FN <- fnr
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- if (denom == 0) 0 else (TP * TN - FP * FN) / denom
  list(
    threshold = threshold,
    fpr = fpr, fnr = fnr, tpr = tpr, tnr = tnr, mcc = mcc
  )
}

#' @noRd
roc_curve_R <- function(d) {
  lr <- d$log10_lr
  p1 <- d$p_h1
  p2 <- d$p_h2
  thr <- sort(unique(lr))
  fpr <- numeric(length(thr))
  tpr <- numeric(length(thr))
  fnr <- numeric(length(thr))
  tnr <- numeric(length(thr))
  for (t in seq_along(thr)) {
    below <- lr < thr[t]
    above <- lr > thr[t]
    fnr[t] <- sum(p1[below])
    tnr[t] <- sum(p2[below])
    tpr[t] <- sum(p1[above])
    fpr[t] <- sum(p2[above])
  }
  ## Concordance AUC with the same i (outer) / j (inner) nesting as
  ## the C++ kernel.
  auc <- 0
  n <- length(lr)
  for (i in seq_len(n)) {
    if (p1[i] <= 0) next
    for (j in seq_len(n)) {
      if (p2[j] <= 0) next
      ind <- if (lr[i] > lr[j]) 1 else if (lr[i] == lr[j]) 0.5 else 0
      auc <- auc + p1[i] * p2[j] * ind
    }
  }
  list(threshold = thr, fpr = fpr, tpr = tpr, fnr = fnr, tnr = tnr,
       auc = auc)
}

## Weighted-Euclidean optimal threshold over the distinct atom values:
## minimise D = sqrt(fnr^2 + (weight*fpr)^2), first minimiser on a tie
## (ascending scan). Analytic analogue of `decision_threshold()`.
#' @noRd
choose_threshold_weighted_R <- function(d, weight) {
  if (!(weight > 0)) stop("weight must be positive.", call. = FALSE)
  lr <- d$log10_lr
  p1 <- d$p_h1
  p2 <- d$p_h2
  thr <- sort(unique(lr))
  best <- NULL
  for (t in thr) {
    fnr <- sum(p1[lr < t])
    fpr <- sum(p2[lr > t])
    dist <- sqrt(fnr^2 + (weight * fpr)^2)
    if (is.null(best) || dist < best$distance) {
      best <- list(threshold = t, fpr = fpr, fnr = fnr, distance = dist)
    }
  }
  best
}

#' @noRd
aggregate_lr_dist <- function(d) {
  if (nrow(d) == 0L) return(d)
  ord <- order(d$log10_lr)
  d <- d[ord, , drop = FALSE]
  rownames(d) <- NULL
  k <- d$log10_lr
  if (length(k) == 1L) {
    grp <- 1L
  } else {
    grp <- cumsum(c(TRUE, k[-1] != k[-length(k)]))
  }
  ph1 <- tapply(d$p_h1, grp, sum)
  ph2 <- tapply(d$p_h2, grp, sum)
  key <- tapply(d$log10_lr, grp, function(x) x[1])
  data.frame(
    log10_lr = unname(as.numeric(key)),
    p_h1 = unname(as.numeric(ph1)),
    p_h2 = unname(as.numeric(ph2))
  )
}
