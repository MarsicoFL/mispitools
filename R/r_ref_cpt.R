## R reference engine — joint CPT for a single marker (F1.2 / F1.4 / F1.5).
##
## Pure-R, brute-force-but-pruned implementation used as the oracle
## against which the C++ engine (F2.x) is later validated. Internal:
## not exported, no public API surface. Slow on large pedigrees by
## design — F1 scope is correctness, not performance.
##
## F1.4 adds equal-rate mutation. Founders are still drawn from HWE on
## the (pre-mutation) population allele frequencies; transmissions from
## parent to child go through a K×K mutation matrix M with
## M[i,i] = 1-R, M[i,j] = R/(K-1) for i != j. The "transmission-aware"
## genotype table T_trans is right-multiplied by M to obtain the
## mutated transmission probabilities, after which the rest of the
## propagation logic is unchanged.
##
## F1.5 adds stepwise mutation, parameterized by `rate` (overall
## per-meiosis mutation probability R) and `ratio` (geometric step
## ratio r in (0, 1)). Allele labels must be coercible to numeric to
## define step distances. M[i, i] = 1 - R; for j != i,
## M[i, j] = (R / sum_{k != i} r^|s_i - s_k|) * r^|s_i - s_j|, where
## s_i are the numeric allele labels. Matches the canonical Familias /
## pedmut / fbnet stepwise mutation matrix when there are no
## microvariant groups (rate2 = 0).
##
## P_H1: joint over all pedigree members under the pedigree topology
## (HWE founders + Mendelian + mutation propagation).
##
## P_H2: P_HWE(g_POI) × marginal over typed members under the
## sub-pedigree obtained by removing POI; if POI was a parent, the
## missing-parent allele on its children is integrated against HWE.
## This matches `P_2(.) = P_1(typed under sub-pedigree) * P_HWE(POI)`
## from the F1 spec and gives a sample space wider than P_H1's support.

#' @noRd
cpt_marker_joint_R <- function(model, poi = NULL) {
  if (!inherits(model, "marker_model")) {
    stop("`model` must be a 'marker_model' object.", call. = FALSE)
  }
  if (!model$mutation$model %in% c("none", "equal", "stepwise")) {
    stop("`cpt_marker_joint_R()` (F1.5 reference engine) supports ",
         "mutation models \"none\", \"equal\", \"stepwise\". The ",
         "asymmetric model arrives in F5.1.",
         call. = FALSE)
  }
  if (!is.null(model$linkage)) {
    stop("`cpt_marker_joint_R()` does not handle linked markers; linkage ",
         "support is added in F5.", call. = FALSE)
  }
  if (!requireNamespace("pedtools", quietly = TRUE)) {
    stop("Package 'pedtools' is required.", call. = FALSE)
  }

  ped <- model$ped
  freqs <- model$freqs
  alleles <- model$alleles
  K <- length(alleles)

  members <- as.character(labels(ped))
  founder_ids <- as.character(pedtools::founders(ped))
  nonfounder_ids <- as.character(pedtools::nonfounders(ped))

  poi_id <- resolve_poi(ped, poi)
  others <- setdiff(members, poi_id)

  mut_matrix <- mutation_matrix_R(model$mutation, K, alleles = alleles)
  prim <- precompute_genotype_tables(freqs, mut_matrix = mut_matrix)
  G <- prim$G

  h1 <- build_joint(ped, prim, founder_ids, nonfounder_ids,
                    excluded = character(0))
  h2 <- build_joint(ped, prim, founder_ids, nonfounder_ids,
                    excluded = poi_id)

  if (length(others) > 0L) {
    h2_keys <- do.call(paste, c(h2$states[, others, drop = FALSE],
                                list(sep = "_")))
    h2_lookup <- stats::setNames(h2$prob, h2_keys)
    h1_others_keys <- do.call(paste, c(h1$states[, others, drop = FALSE],
                                       list(sep = "_")))
    P_H2_for_h1 <- prim$hwe[h1$states[[poi_id]]] *
      ifelse(h1_others_keys %in% names(h2_lookup),
             unname(h2_lookup[h1_others_keys]), 0)

    n_h2 <- nrow(h2$states)
    h2_expanded <- h2$states[rep(seq_len(n_h2), each = G), , drop = FALSE]
    h2_expanded[[poi_id]] <- rep(seq_len(G), times = n_h2)
    h2_expanded <- h2_expanded[, members, drop = FALSE]
    rownames(h2_expanded) <- NULL
    P_H2_h2 <- prim$hwe[h2_expanded[[poi_id]]] * rep(h2$prob, each = G)

    h1_keys <- do.call(paste, c(h1$states[, members, drop = FALSE],
                                list(sep = "_")))
    h2_full_keys <- do.call(paste, c(h2_expanded, list(sep = "_")))
    new_in_h2 <- !(h2_full_keys %in% h1_keys)

    combined_states <- rbind(h1$states, h2_expanded[new_in_h2, , drop = FALSE])
    combined_P_H1 <- c(h1$prob, rep(0.0, sum(new_in_h2)))
    combined_P_H2 <- c(P_H2_for_h1, P_H2_h2[new_in_h2])
  } else {
    combined_states <- h1$states
    combined_P_H1 <- h1$prob
    combined_P_H2 <- prim$hwe[h1$states[[poi_id]]]
  }

  keep <- (combined_P_H1 > 0) | (combined_P_H2 > 0)
  combined_states <- combined_states[keep, , drop = FALSE]
  combined_P_H1 <- combined_P_H1[keep]
  combined_P_H2 <- combined_P_H2[keep]

  out <- data.frame(matrix(NA_character_, nrow = nrow(combined_states),
                           ncol = length(members)),
                    stringsAsFactors = FALSE)
  names(out) <- members
  for (m in members) {
    out[[m]] <- prim$geno_labels[combined_states[[m]]]
  }
  out$P_H1 <- combined_P_H1
  out$P_H2 <- combined_P_H2

  ord_keys <- lapply(members, function(m) combined_states[[m]])
  ord <- do.call(order, ord_keys)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  attr(out, "marker_id") <- model$marker_id
  attr(out, "alleles") <- alleles
  attr(out, "poi") <- poi_id
  out
}

#' @noRd
mutation_matrix_R <- function(mutation, K, alleles = NULL) {
  if (mutation$model == "none") {
    return(diag(K))
  }
  if (mutation$model == "equal") {
    R <- mutation$rate
    if (K < 2L) {
      stop("Equal-rate mutation requires at least 2 alleles.", call. = FALSE)
    }
    M <- matrix(R / (K - 1L), nrow = K, ncol = K)
    diag(M) <- 1 - R
    return(M)
  }
  if (mutation$model == "stepwise") {
    R <- mutation$rate
    r <- mutation$ratio
    if (is.null(r)) {
      stop("Stepwise mutation requires `mutation$ratio` (geometric step ",
           "ratio in (0, 1)).", call. = FALSE)
    }
    if (K < 2L) {
      stop("Stepwise mutation requires at least 2 alleles.", call. = FALSE)
    }
    if (is.null(alleles) || length(alleles) != K) {
      stop("Stepwise mutation requires `alleles` (the K allele labels) ",
           "to compute step distances.", call. = FALSE)
    }
    s <- suppressWarnings(as.numeric(alleles))
    if (anyNA(s)) {
      stop("Stepwise mutation requires numeric allele labels; got non-",
           "numeric: ", paste(alleles[is.na(s)], collapse = ", "),
           call. = FALSE)
    }
    M <- matrix(0.0, nrow = K, ncol = K)
    for (i in seq_len(K)) {
      steps <- abs(s - s[i])
      w <- r ^ steps
      w[i] <- 0
      sw <- sum(w)
      if (!is.finite(sw) || sw <= 0) {
        stop("Stepwise mutation: row weights for allele \"", alleles[i],
             "\" sum to 0; cannot normalize.", call. = FALSE)
      }
      M[i, ] <- (R / sw) * w
      M[i, i] <- 1 - R
    }
    return(M)
  }
  stop("Internal error: mutation model \"", mutation$model,
       "\" not supported by mutation_matrix_R().", call. = FALSE)
}

#' @noRd
precompute_genotype_tables <- function(freqs, mut_matrix = NULL) {
  K <- length(freqs)
  alleles <- names(freqs)
  geno_idx <- expand.grid(a1 = seq_len(K), a2 = seq_len(K),
                          KEEP.OUT.ATTRS = FALSE)
  geno_idx <- geno_idx[geno_idx$a1 <= geno_idx$a2, , drop = FALSE]
  rownames(geno_idx) <- NULL
  G <- nrow(geno_idx)
  geno_labels <- paste0(alleles[geno_idx$a1], "/", alleles[geno_idx$a2])

  hwe <- ifelse(geno_idx$a1 == geno_idx$a2,
                freqs[geno_idx$a1]^2,
                2 * freqs[geno_idx$a1] * freqs[geno_idx$a2])
  hwe <- as.numeric(hwe)

  T_trans <- matrix(0.0, nrow = G, ncol = K)
  for (g in seq_len(G)) {
    i <- geno_idx$a1[g]
    j <- geno_idx$a2[g]
    if (i == j) {
      T_trans[g, i] <- 1.0
    } else {
      T_trans[g, i] <- 0.5
      T_trans[g, j] <- 0.5
    }
  }
  if (!is.null(mut_matrix)) {
    if (!is.matrix(mut_matrix) || nrow(mut_matrix) != K ||
        ncol(mut_matrix) != K) {
      stop("`mut_matrix` must be a K x K numeric matrix.", call. = FALSE)
    }
    T_trans <- T_trans %*% mut_matrix
  }

  child_dist <- array(0.0, dim = c(G, G, G))
  for (gp in seq_len(G)) {
    for (gm in seq_len(G)) {
      for (gc in seq_len(G)) {
        cA <- geno_idx$a1[gc]
        cB <- geno_idx$a2[gc]
        if (cA == cB) {
          child_dist[gp, gm, gc] <- T_trans[gp, cA] * T_trans[gm, cA]
        } else {
          child_dist[gp, gm, gc] <- T_trans[gp, cA] * T_trans[gm, cB] +
                                    T_trans[gp, cB] * T_trans[gm, cA]
        }
      }
    }
  }

  child_dist_one_missing <- matrix(0.0, nrow = G, ncol = G)
  for (g_K in seq_len(G)) {
    for (g_C in seq_len(G)) {
      child_dist_one_missing[g_K, g_C] <- sum(hwe * child_dist[g_K, , g_C])
    }
  }

  list(K = K, G = G, geno_idx = geno_idx, geno_labels = geno_labels,
       hwe = hwe, T_trans = T_trans, child_dist = child_dist,
       child_dist_one_missing = child_dist_one_missing)
}

#' @noRd
build_joint <- function(ped, prim, founder_ids, nonfounder_ids, excluded) {
  G <- prim$G
  hwe <- prim$hwe
  child_dist <- prim$child_dist
  child_dist_one_missing <- prim$child_dist_one_missing

  founder_ids_use <- setdiff(founder_ids, excluded)
  nonfounder_ids_use <- setdiff(nonfounder_ids, excluded)

  if (length(founder_ids_use) == 0L) {
    states <- data.frame(row.names = seq_len(0L))
    prob <- numeric(0)
  } else {
    grids <- lapply(founder_ids_use, function(.) seq_len(G))
    names(grids) <- founder_ids_use
    states <- do.call(expand.grid,
                      c(grids,
                        list(KEEP.OUT.ATTRS = FALSE,
                             stringsAsFactors = FALSE)))
    prob <- rep(1.0, nrow(states))
    for (fid in founder_ids_use) {
      prob <- prob * hwe[states[[fid]]]
    }
  }

  nf_order <- ancestral_order(ped, nonfounder_ids_use, excluded = excluded)
  for (nf in nf_order) {
    fa <- as.character(pedtools::father(ped, id = nf))
    mo <- as.character(pedtools::mother(ped, id = nf))
    n_states <- nrow(states)
    if (n_states == 0L) break

    expanded <- states[rep(seq_len(n_states), each = G), , drop = FALSE]
    expanded[[nf]] <- rep(seq_len(G), times = n_states)
    rownames(expanded) <- NULL
    gc <- expanded[[nf]]

    fa_known <- !(fa %in% excluded)
    mo_known <- !(mo %in% excluded)

    if (fa_known && mo_known) {
      gp <- expanded[[fa]]
      gm <- expanded[[mo]]
      flat_idx <- gp + (gm - 1L) * G + (gc - 1L) * G * G
      cond <- child_dist[flat_idx]
    } else if (fa_known && !mo_known) {
      gp <- expanded[[fa]]
      flat_idx <- gp + (gc - 1L) * G
      cond <- child_dist_one_missing[flat_idx]
    } else if (!fa_known && mo_known) {
      gm <- expanded[[mo]]
      flat_idx <- gm + (gc - 1L) * G
      cond <- child_dist_one_missing[flat_idx]
    } else {
      stop("Internal error: nonfounder \"", nf, "\" has both parents ",
           "excluded; cannot impute.", call. = FALSE)
    }

    new_prob <- rep(prob, each = G) * cond
    keep <- new_prob > 0
    states <- expanded[keep, , drop = FALSE]
    prob <- new_prob[keep]
    rownames(states) <- NULL
  }

  list(states = states, prob = prob)
}

#' @noRd
ancestral_order <- function(ped, nonfounder_ids, excluded = character(0)) {
  if (length(nonfounder_ids) == 0L) return(character(0))
  ## Treat excluded members as already-resolved (their alleles are HWE-imputed
  ## for any descendants).
  assigned <- c(setdiff(as.character(pedtools::founders(ped)), excluded),
                as.character(excluded))
  remaining <- as.character(nonfounder_ids)
  result <- character(0)
  max_iter <- length(remaining) + 1L
  for (iter in seq_len(max_iter)) {
    if (length(remaining) == 0L) break
    progress <- FALSE
    for (nf in remaining) {
      fa <- as.character(pedtools::father(ped, id = nf))
      mo <- as.character(pedtools::mother(ped, id = nf))
      if (fa %in% assigned && mo %in% assigned) {
        result <- c(result, nf)
        assigned <- c(assigned, nf)
        remaining <- setdiff(remaining, nf)
        progress <- TRUE
      }
    }
    if (!progress) {
      stop("Cannot determine ancestral order; pedigree may have a cycle.",
           call. = FALSE)
    }
  }
  result
}

#' @noRd
resolve_poi <- function(ped, poi) {
  members <- as.character(labels(ped))
  if (!is.null(poi)) {
    if (!is.character(poi) || length(poi) != 1L || is.na(poi) ||
        !nzchar(poi)) {
      stop("`poi` must be a single non-empty character string.",
           call. = FALSE)
    }
    if (!poi %in% members) {
      stop(sprintf("`poi` = \"%s\" is not a member of the pedigree.",
                   poi), call. = FALSE)
    }
    return(poi)
  }
  nonfounders <- as.character(pedtools::nonfounders(ped))
  typed <- tryCatch(as.character(pedtools::typedMembers(ped)),
                    error = function(e) character(0))
  candidates <- setdiff(nonfounders, typed)
  if (length(candidates) >= 1L) {
    return(candidates[length(candidates)])
  }
  if (length(nonfounders) >= 1L) {
    return(nonfounders[length(nonfounders)])
  }
  members[length(members)]
}
