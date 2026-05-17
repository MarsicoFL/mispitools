## F4.3 — S3 layer over a sparse LR distribution.
##
## `lr_dist` wraps the (log10_lr, p_h1, p_h2) atoms produced by the
## C++ kernel (per_marker_lr_dist / lr_dist_compose). The decision-
## theoretic primitives (moments, quantiles, ROC/FPR/FNR) live in
## src/core/decision.cpp and are surfaced here through `summary()` and
## `plot()`. The model-aware constructor `lr_distribution(model_list)`
## and a `quantile()` method arrive in F4.4 and build on this class.

#' @noRd
new_lr_dist <- function(log10_lr, p_h1, p_h2, marker_id = NA_character_,
                        poi = NULL) {
  out <- data.frame(
    log10_lr = as.numeric(log10_lr),
    p_h1 = as.numeric(p_h1),
    p_h2 = as.numeric(p_h2)
  )
  attr(out, "marker_id") <- marker_id
  attr(out, "poi") <- poi
  class(out) <- c("lr_dist", "data.frame")
  out
}

#' Coerce to a Likelihood-Ratio Distribution Object
#'
#' @description
#' Wraps a sparse \eqn{(\log_{10}\mathrm{LR}, P(\cdot\mid H_1),
#' P(\cdot\mid H_2))} table in an `lr_dist` object so the decision-
#' theoretic `summary()` and `plot()` methods can be applied. This is
#' the low-level constructor; the model-aware `lr_distribution()`
#' wrapper is added in a later release and produces the same class.
#'
#' @param x A data.frame or list with numeric components `log10_lr`,
#'   `p_h1` and `p_h2`, all of the same length. Typically the output of
#'   the per-marker LR distribution engine.
#'
#' @return An object of class `lr_dist` (a data.frame with the three
#'   columns), suitable for [summary.lr_dist()] and [plot.lr_dist()].
#'
#' @details
#' Each row is an atom of the discrete LR distribution: `log10_lr` is
#' the log10 likelihood ratio, `p_h1` its probability under \eqn{H_1}
#' (the person of interest is the missing person) and `p_h2` its
#' probability under \eqn{H_2} (unrelated). `+Inf` / `-Inf` atoms
#' (which arise under `mutation = "none"`) are allowed and handled by
#' the summaries with the limit \eqn{0\log 0 = 0}.
#'
#' @seealso [summary.lr_dist()], [plot.lr_dist()]
#'
#' @export
#' @examples
#' d <- data.frame(
#'   log10_lr = c(-1, 0, 2),
#'   p_h1 = c(0.1, 0.3, 0.6),
#'   p_h2 = c(0.6, 0.3, 0.1)
#' )
#' x <- as_lr_dist(d)
#' summary(x)
as_lr_dist <- function(x) {
  if (inherits(x, "lr_dist")) {
    return(x)
  }
  if (!is.list(x)) {
    stop("`x` must be a data.frame or list.", call. = FALSE)
  }
  nm <- names(x)
  if (is.null(nm) || !all(c("log10_lr", "p_h1", "p_h2") %in% nm)) {
    stop("`x` must have components 'log10_lr', 'p_h1' and 'p_h2'.",
         call. = FALSE)
  }
  lr <- x[["log10_lr"]]
  p1 <- x[["p_h1"]]
  p2 <- x[["p_h2"]]
  if (length(lr) != length(p1) || length(lr) != length(p2)) {
    stop("'log10_lr', 'p_h1' and 'p_h2' must have the same length.",
         call. = FALSE)
  }
  if (any(p1 < 0) || any(p2 < 0)) {
    stop("probabilities must be non-negative.", call. = FALSE)
  }
  new_lr_dist(lr, p1, p2,
              marker_id = attr(x, "marker_id") %||% NA_character_,
              poi = attr(x, "poi"))
}

#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @export
print.lr_dist <- function(x, ...) {
  mid <- attr(x, "marker_id")
  hdr <- "<lr_dist>"
  if (!is.null(mid) && !is.na(mid) && nzchar(mid)) {
    hdr <- paste0(hdr, " marker: ", mid)
  }
  cat(hdr, "\n", sep = "")
  cat(nrow(x), " atoms; ",
      "E[log10 LR | H1] = ",
      format(stats::weighted.mean(
        x$log10_lr[x$p_h1 > 0], x$p_h1[x$p_h1 > 0]), digits = 4),
      "\n", sep = "")
  print(utils::head(as.data.frame(x), 10L), row.names = FALSE)
  if (nrow(x) > 10L) cat("... (", nrow(x) - 10L, " more)\n", sep = "")
  invisible(x)
}

#' Summarise a Likelihood-Ratio Distribution
#'
#' @description
#' Computes the decision-theoretic summary of an `lr_dist` object: the
#' mean, variance and standard deviation of \eqn{\log_{10}\mathrm{LR}}
#' under each hypothesis, the median and quartiles, the area under the
#' ROC curve, and the total probability mass (a sanity check on the
#' input joint).
#'
#' @param object An `lr_dist` object (see [as_lr_dist()]).
#' @param ... Unused.
#'
#' @return A list of class `summary.lr_dist` with components `mean_h1`,
#'   `mean_h2`, `var_h1`, `var_h2`, `sd_h1`, `sd_h2`, `mass_h1`,
#'   `mass_h2`, `auc`, `quantiles_h1`, `quantiles_h2`, `has_pos_inf`
#'   and `has_neg_inf`. Printed in a compact table.
#'
#' @details
#' `mean_h1` is \eqn{E[\log_{10}\mathrm{LR}\mid H_1]} and equals the
#' per-marker KL-derived expectation. A `+Inf` atom drives the H1
#' moments to `+Inf` (and a `-Inf` atom the H2 moments to `-Inf`),
#' matching the engine's `0\log 0 = 0` convention. `auc` is the exact
#' concordance statistic \eqn{P(\mathrm{LR}_{H_1} > \mathrm{LR}_{H_2})
#' + \tfrac12 P(\mathrm{LR}_{H_1} = \mathrm{LR}_{H_2})}.
#'
#' @seealso [plot.lr_dist()], [as_lr_dist()]
#'
#' @export
#' @examples
#' d <- as_lr_dist(data.frame(
#'   log10_lr = c(-1, 0, 2),
#'   p_h1 = c(0.1, 0.3, 0.6),
#'   p_h2 = c(0.6, 0.3, 0.1)
#' ))
#' s <- summary(d)
#' s$mean_h1
#' s$auc
summary.lr_dist <- function(object, ...) {
  s <- cpp_lr_dist_summary(object$log10_lr, object$p_h1, object$p_h2)
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  q1 <- cpp_lr_dist_quantile(object$log10_lr, object$p_h1, object$p_h2,
                             probs, under_h1 = TRUE)
  q2 <- cpp_lr_dist_quantile(object$log10_lr, object$p_h1, object$p_h2,
                             probs, under_h1 = FALSE)
  names(q1) <- names(q2) <- paste0(probs * 100, "%")
  roc <- cpp_lr_dist_roc(object$log10_lr, object$p_h1, object$p_h2)
  out <- c(s, list(auc = roc$auc, quantiles_h1 = q1, quantiles_h2 = q2))
  class(out) <- "summary.lr_dist"
  out
}

#' @export
print.summary.lr_dist <- function(x, ...) {
  cat("Likelihood-ratio distribution summary\n")
  tab <- data.frame(
    H1 = c(x$mean_h1, x$var_h1, x$sd_h1, x$mass_h1),
    H2 = c(x$mean_h2, x$var_h2, x$sd_h2, x$mass_h2),
    row.names = c("E[log10 LR]", "Var", "SD", "mass")
  )
  print(round(tab, 6))
  cat("\nAUC: ", format(x$auc, digits = 6), "\n", sep = "")
  cat("Quantiles of log10 LR | H1:\n")
  print(round(x$quantiles_h1, 4))
  cat("Quantiles of log10 LR | H2:\n")
  print(round(x$quantiles_h2, 4))
  if (isTRUE(x$has_pos_inf)) {
    cat("Note: a +Inf atom is present (H2 lacks an H1 state).\n")
  }
  if (isTRUE(x$has_neg_inf)) {
    cat("Note: a -Inf atom is present (H1 lacks an H2 state).\n")
  }
  invisible(x)
}

#' Plot a Likelihood-Ratio Distribution
#'
#' @description
#' Draws the discrete \eqn{\log_{10}\mathrm{LR}} distribution under both
#' hypotheses as an overlaid lollipop plot: \eqn{H_1} (related, blue)
#' and \eqn{H_2} (unrelated, red). Less overlap means stronger
#' discrimination.
#'
#' @param x An `lr_dist` object (see [as_lr_dist()]).
#' @param ... Unused.
#'
#' @return A `ggplot2` object.
#'
#' @details
#' `+Inf` / `-Inf` atoms cannot be placed on a finite axis and are
#' dropped from the plot with a message; their mass is still reported by
#' [summary.lr_dist()]. The x-axis is \eqn{\log_{10}\mathrm{LR}}
#' (`0` = neutral, `> 0` favours \eqn{H_1}).
#'
#' @seealso [summary.lr_dist()], [as_lr_dist()]
#'
#' @import ggplot2
#' @export
#' @examples
#' d <- as_lr_dist(data.frame(
#'   log10_lr = c(-1, 0, 2),
#'   p_h1 = c(0.1, 0.3, 0.6),
#'   p_h2 = c(0.6, 0.3, 0.1)
#' ))
#' plot(d)
plot.lr_dist <- function(x, ...) {
  finite <- is.finite(x$log10_lr)
  if (!all(finite)) {
    message(sum(!finite),
            " infinite atom(s) dropped from the plot (see summary()).")
  }
  d <- x[finite, , drop = FALSE]
  if (nrow(d) == 0L) {
    stop("no finite atoms to plot.", call. = FALSE)
  }
  long <- data.frame(
    log10_lr = rep(d$log10_lr, 2L),
    prob = c(d$p_h1, d$p_h2),
    Hypothesis = factor(rep(c("H1 (related)", "H2 (unrelated)"),
                            each = nrow(d)),
                        levels = c("H1 (related)", "H2 (unrelated)"))
  )
  ggplot2::ggplot(
    long,
    ggplot2::aes(x = log10_lr, y = prob, colour = Hypothesis)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = log10_lr, yend = 0),
      linewidth = 0.6
    ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(
      values = c("H1 (related)" = "#2C7FB8",
                 "H2 (unrelated)" = "#D7301F")
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::labs(x = expression(log[10] * "(LR)"),
                  y = "Probability", colour = NULL) +
    ggplot2::theme_minimal()
}
