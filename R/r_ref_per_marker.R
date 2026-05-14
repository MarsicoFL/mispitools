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
