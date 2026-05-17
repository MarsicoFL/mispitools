## F2.2 — R-side wrapper for the C++ kernel cpt_marker_joint_cpp().
##
## Mirrors cpt_marker_joint_R(): consumes a `marker_model`, builds the
## POD inputs that the C++ engine expects, runs the kernel, and decorates
## the output with genotype-label strings + attributes. The returned
## data.frame is row-aligned with the R-reference output so the F2.6
## cross-check can compare element-by-element.

#' @noRd
cpt_marker_joint_cpp_wrap <- function(model, poi = NULL) {
  if (!inherits(model, "marker_model")) {
    stop("`model` must be a 'marker_model' object.", call. = FALSE)
  }
  if (!model$mutation$model %in% c("none", "equal", "stepwise")) {
    stop("`cpt_marker_joint_cpp_wrap()` supports mutation models ",
         "\"none\", \"equal\", \"stepwise\". The asymmetric (Dawid) ",
         "matrix builder exists in the core (F5.1) but is not yet ",
         "routed through this high-level wrapper.", call. = FALSE)
  }
  if (!is.null(model$linkage)) {
    stop("`cpt_marker_joint_cpp_wrap()` does not handle linked markers; ",
         "linkage support arrives in F5.", call. = FALSE)
  }
  if (!requireNamespace("pedtools", quietly = TRUE)) {
    stop("Package 'pedtools' is required.", call. = FALSE)
  }

  ped <- model$ped
  members <- as.character(labels(ped))
  n <- length(members)
  name_to_idx <- stats::setNames(seq_len(n) - 1L, members)

  father <- integer(n)
  mother <- integer(n)
  for (i in seq_len(n)) {
    fa <- as.character(pedtools::father(ped, id = members[i]))
    mo <- as.character(pedtools::mother(ped, id = members[i]))
    father[i] <- if (length(fa) == 0L || !nzchar(fa) || is.na(fa)) -1L else
      unname(name_to_idx[[fa]])
    mother[i] <- if (length(mo) == 0L || !nzchar(mo) || is.na(mo)) -1L else
      unname(name_to_idx[[mo]])
  }

  poi_id <- resolve_poi(ped, poi)
  poi_idx <- unname(name_to_idx[[poi_id]])

  freqs <- model$freqs
  alleles <- model$alleles
  K <- length(alleles)

  mut <- model$mutation
  mut_kind   <- switch(mut$model, none = 0L, equal = 1L, stepwise = 2L)
  mut_rate   <- if (mut$model == "none") 0.0 else as.numeric(mut$rate)
  mut_range  <- if (mut$model == "stepwise") as.numeric(mut$ratio) else 0.0
  mut_labels <- if (mut$model == "stepwise") {
    s <- suppressWarnings(as.numeric(alleles))
    if (anyNA(s)) {
      stop("Stepwise mutation requires numeric allele labels; got non-",
           "numeric: ", paste(alleles[is.na(s)], collapse = ", "),
           call. = FALSE)
    }
    s
  } else {
    rep(NA_real_, K)
  }

  res <- cpt_marker_joint_cpp(
    father = as.integer(father),
    mother = as.integer(mother),
    poi = as.integer(poi_idx),
    freqs = as.numeric(unname(freqs)),
    mutation_kind = mut_kind,
    mutation_rate = mut_rate,
    mutation_range = mut_range,
    mutation_rate2 = 0.0,
    mutation_bias = 0.5,
    numeric_labels = mut_labels
  )

  G <- res$n_genotypes
  geno_labels <- character(G)
  g <- 0L
  for (j in seq_len(K)) {
    for (i in seq_len(j)) {
      g <- g + 1L
      geno_labels[g] <- paste0(alleles[i], "/", alleles[j])
    }
  }

  states <- res$states
  n_rows <- nrow(states)
  out <- data.frame(matrix(NA_character_, nrow = n_rows, ncol = n),
                    stringsAsFactors = FALSE)
  names(out) <- members
  for (mi in seq_len(n)) {
    out[[mi]] <- geno_labels[states[, mi]]
  }
  out$P_H1 <- res$P_H1
  out$P_H2 <- res$P_H2

  attr(out, "marker_id") <- model$marker_id
  attr(out, "alleles") <- alleles
  attr(out, "poi") <- poi_id
  out
}
