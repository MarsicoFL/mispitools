library(shiny)
library(shinythemes)
library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)

Hypothesis <- logLR <- w <- Var1 <- Var2 <- value <- label <- TPR <- FPR <- FNR <- NULL

DEFAULT_HAIR_FREQ <- c(
  Black = 0.08,
  Brown = 0.45,
  Blonde = 0.30,
  Red = 0.07,
  Gray = 0.10
)

hair_labels <- c("1" = "Black", "2" = "Brown", "3" = "Blonde",
                 "4" = "Red", "5" = "Gray/White")

hair_choices <- c("Black" = 1, "Brown" = 2, "Blonde" = 3, "Red" = 4, "Gray/White" = 5)

sex_choices <- c("Female" = "F", "Male" = "M")

hair_error_help <- "Total probability of observing a wrong color, split equally among the other 4 colors."

check_num <- function(x, label, lo, hi, whole = FALSE) {
  if (is.null(x) || length(x) != 1 || !is.numeric(x) || !is.finite(x)) {
    return(paste0(label, ": please enter a number."))
  }
  if (x < lo || x > hi) {
    if (is.infinite(hi)) return(sprintf("%s: must be at least %s.", label, lo))
    return(sprintf("%s: must be between %s and %s.", label, lo, hi))
  }
  if (whole && x != round(x)) {
    return(paste0(label, ": must be a whole number."))
  }
  NULL
}

age_window <- function(MPa, MPr) {
  MPmin <- max(1, MPa - MPr)
  MPmax <- min(80, MPa + MPr)
  list(min = MPmin, max = MPmax, T1p = (MPmax - MPmin + 1) / 80)
}

format_lr <- function(LR, digits = 4) {
  if (LR == 0) "LR = 0" else paste("LR =", sprintf(paste0("%.", digits, "g"), LR))
}

format_log_lr <- function(LR) {
  if (LR == 0) "LR = 0 (exclusion)" else sprintf("log10(LR) = %.3f", log10(LR))
}

build_cpt_pop <- function(propF, MPa, MPr, propC) {
  T1p <- age_window(MPa, MPr)$T1p
  T0p <- 1 - T1p
  propS <- c(propF, 1 - propF)
  jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")
  jointprob <- c(propS[1]*T1p, propS[1]*T0p, propS[2]*T1p, propS[2]*T0p)
  names(jointprob) <- jointname
  CPT_POP <- outer(jointprob, propC)
  rownames(CPT_POP) <- jointname
  colnames(CPT_POP) <- names(propC)
  CPT_POP
}

build_cpt_mp <- function(MPs, MPc, eps, epa, epc) {
  jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")
  if (MPs == "F") {
    jointprob <- c((1-eps)*(1-epa), (1-eps)*epa, eps*(1-epa), eps*epa)
  } else {
    jointprob <- c(eps*(1-epa), eps*epa, (1-eps)*(1-epa), (1-eps)*epa)
  }
  names(jointprob) <- jointname
  n_colors <- 5
  errorMat <- matrix(epc / (n_colors - 1), nrow = n_colors, ncol = n_colors)
  diag(errorMat) <- 1 - epc
  probC <- errorMat[MPc, ]
  CPT_MP <- outer(jointprob, probC)
  rownames(CPT_MP) <- jointname
  colnames(CPT_MP) <- names(DEFAULT_HAIR_FREQ)
  CPT_MP
}

create_cpt_plots <- function(CPT_POP, CPT_MP, narrow = FALSE) {
  dfPOP <- reshape2::melt(CPT_POP)
  dfMP <- reshape2::melt(CPT_MP)
  dfLR <- reshape2::melt(log10(CPT_MP / CPT_POP))
  dfLR$label <- ifelse(is.finite(dfLR$value), sprintf("%.2f", dfLR$value), "LR = 0")
  dfLR$value[!is.finite(dfLR$value)] <- NA

  base_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, color = "gray40"),
      axis.text = ggplot2::element_text(size = 9),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.9, "cm"),
      plot.title.position = "plot"
    )

  p1 <- ggplot2::ggplot(dfPOP, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient(low = "#ffffff", high = "#7f93ad",
                                 limits = c(0, max(dfPOP$value) * 1.1)) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)),
                       size = 2.8, color = "black") +
    ggplot2::labs(title = "P(D | H2)", subtitle = "Reference population",
                  x = "Hair Color", y = "", fill = "Probability") +
    base_theme

  p2 <- ggplot2::ggplot(dfMP, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient(low = "#ffffff", high = "#7f93ad",
                                 limits = c(0, max(dfMP$value) * 1.1)) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)),
                       size = 2.8, color = "black") +
    ggplot2::labs(title = "P(D | H1)", subtitle = "If person IS the MP",
                  x = "Hair Color", y = "", fill = "Probability") +
    base_theme

  p3 <- ggplot2::ggplot(dfLR, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(low = "#c9827f", mid = "#ffffff", high = "#7f93ad",
                                  midpoint = 0, na.value = "gray80",
                                  labels = function(x) sprintf("%.1f", x)) +
    ggplot2::geom_text(ggplot2::aes(label = label),
                       size = 2.8, color = "black") +
    ggplot2::labs(title = expression(log[10](LR)),
                  subtitle = "Evidence weight (grey: LR = 0)",
                  x = "Hair Color", y = "", fill = expression(log[10](LR))) +
    base_theme

  p1 + p2 + p3 + patchwork::plot_layout(ncol = if (narrow) 1 else 3)
}

build_lr_distributions <- function(CPT_POP, CPT_MP) {
  pop_vec <- as.vector(CPT_POP)
  mp_vec <- as.vector(CPT_MP)
  log_lr <- round(log10(mp_vec / pop_vec), 10)
  df_H1 <- data.frame(x = log_lr, w = mp_vec)
  df_H2 <- data.frame(x = log_lr, w = pop_vec)
  list(H1 = df_H1[df_H1$w > 0, ], H2 = df_H2[df_H2$w > 0, ])
}

create_dist_plot <- function(distH1, distH2, binwidth = 0.25, base_size = 12) {
  if (nrow(distH1) == 0 || nrow(distH2) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") +
             ggplot2::theme_void())
  }

  dH1 <- data.frame(logLR = distH1$x, w = distH1$w, Hypothesis = "H1: Person IS the MP")
  dH2 <- data.frame(logLR = distH2$x, w = distH2$w, Hypothesis = "H2: Person is NOT the MP")
  dAll <- rbind(dH1, dH2)
  dFinite <- dAll[is.finite(dAll$logLR), ]
  dZero <- dAll[dAll$logLR == -Inf, ]

  p <- ggplot2::ggplot(dFinite, ggplot2::aes(x = logLR, weight = w, fill = Hypothesis)) +
    ggplot2::geom_histogram(binwidth = binwidth, position = "identity",
                            alpha = 0.65, color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c("H1: Person IS the MP" = "#1f3a5f",
                                          "H2: Person is NOT the MP" = "#999999")) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
    ggplot2::labs(title = expression("Distribution of " * log[10](LR)),
                  subtitle = "Vertical line: LR = 1 (uninformative)",
                  x = expression(log[10](LR)), y = "Probability") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title.position = "plot",
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (nrow(dZero) > 0) {
    dZero <- stats::aggregate(w ~ Hypothesis, data = dZero, FUN = sum)
    dZero$logLR <- min(dFinite$logLR) - 3 * binwidth
    p <- p +
      ggplot2::geom_col(data = dZero, ggplot2::aes(x = logLR, y = w, fill = Hypothesis),
                        inherit.aes = FALSE, width = binwidth, alpha = 0.65,
                        color = "white", linewidth = 0.3) +
      ggplot2::geom_text(data = dZero,
                         ggplot2::aes(x = logLR, y = w, label = sprintf("LR = 0\nP = %.3f", w)),
                         inherit.aes = FALSE, vjust = -0.3, size = 3.2) +
      ggplot2::geom_vline(xintercept = dZero$logLR[1] + 1.5 * binwidth,
                          linetype = "dotted", color = "gray50") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
      ggplot2::labs(subtitle = "Vertical line: LR = 1 (uninformative). Left bar: LR = 0 (exclusion)")
  }
  p
}

calc_metrics <- function(distH1, distH2, threshold) {
  if (nrow(distH1) == 0 || nrow(distH2) == 0) {
    return(list(TPR = NA, FNR = NA, FPR = NA, TNR = NA, MCC = NA))
  }

  x1 <- distH1$x
  w1 <- distH1$w
  x2 <- distH2$x
  w2 <- distH2$w

  TP <- sum(w1[x1 >= threshold])
  FN <- sum(w1[x1 < threshold])
  FP <- sum(w2[x2 >= threshold])
  TN <- sum(w2[x2 < threshold])

  denom <- sqrt(as.numeric(TP+FP) * as.numeric(TP+FN) *
                as.numeric(TN+FP) * as.numeric(TN+FN))
  MCC <- if (denom == 0) 0 else (TP*TN - FP*FN) / denom

  list(TPR = TP, FNR = FN, FPR = FP, TNR = TN, MCC = MCC)
}

build_roc <- function(distH1, distH2) {
  if (nrow(distH1) == 0 || nrow(distH2) == 0) {
    return(list(data = data.frame(th = numeric(0), TPR = numeric(0), FPR = numeric(0)),
                AUC = NA))
  }

  thresholds <- c(Inf, sort(unique(c(distH1$x, distH2$x)), decreasing = TRUE))
  TPR_vec <- vapply(thresholds, function(t) sum(distH1$w[distH1$x >= t]), numeric(1))
  FPR_vec <- vapply(thresholds, function(t) sum(distH2$w[distH2$x >= t]), numeric(1))
  n <- length(thresholds)
  AUC_value <- sum(diff(FPR_vec) * (TPR_vec[-1] + TPR_vec[-n]) / 2)

  list(data = data.frame(th = thresholds, TPR = TPR_vec, FPR = FPR_vec),
       AUC = AUC_value)
}

calc_lr_sex <- function(MPs, obsSex, eps, propF) {
  propS <- c(F = propF, M = 1 - propF)
  if (MPs == obsSex) {
    LR <- (1 - eps) / propS[obsSex]
  } else {
    LR <- eps / propS[obsSex]
  }
  unname(LR)
}

calc_lr_age <- function(MPa, MPr, obsAge, epa, T1p = age_window(MPa, MPr)$T1p) {
  aw <- age_window(MPa, MPr)
  if (obsAge >= aw$min && obsAge <= aw$max) {
    (1 - epa) / T1p
  } else {
    epa / (1 - T1p)
  }
}

calc_lr_hair <- function(MPc, obsC, epc, propC) {
  n_colors <- length(propC)
  if (MPc == obsC) {
    p_obs_h1 <- 1 - epc
  } else {
    p_obs_h1 <- epc / (n_colors - 1)
  }
  unname(p_obs_h1 / propC[obsC])
}

parse_sim_upload <- function(filepath) {
  if (is.null(filepath) || !file.exists(filepath)) return(NULL)

  head_bytes <- readBin(filepath, "raw", n = 4096)
  if (length(head_bytes) < 2 || any(head_bytes == as.raw(0)) ||
      identical(head_bytes[1:2], charToRaw("PK"))) return(NULL)

  lines <- tryCatch(readLines(filepath, warn = FALSE, encoding = "UTF-8"),
                    error = function(e) NULL)
  if (is.null(lines) || length(lines) < 3) return(NULL)
  if (!all(validUTF8(lines))) lines <- iconv(lines, "latin1", "UTF-8")
  if (anyNA(lines)) return(NULL)

  lines <- gsub("\r", "", lines, fixed = TRUE)
  while (length(lines) > 0 && !nzchar(trimws(lines[length(lines)]))) {
    lines <- lines[-length(lines)]
  }

  result <- .try_famlink(lines)
  if (!is.null(result)) return(result)

  result <- .try_familias(lines)
  if (!is.null(result)) return(result)

  result <- .try_csv_generic(filepath, lines)
  if (!is.null(result)) return(result)

  NULL
}

.try_famlink <- function(lines) {
  head5 <- lines[1:min(5, length(lines))]
  if (!any(grepl("Raw data LR|Pedigree \\(Numerator\\)", head5, ignore.case = TRUE)))
    return(NULL)

  num_idx <- grep("^Pedigree \\(Numerator\\)", lines, ignore.case = TRUE)
  den_idx <- grep("^Versus \\(Den", lines, ignore.case = TRUE)
  true_idx <- grep("^True:", lines, ignore.case = TRUE)

  if (length(num_idx) == 0 || length(den_idx) == 0 || length(true_idx) == 0) return(NULL)
  num_idx <- num_idx[1]; den_idx <- den_idx[1]; true_idx <- true_idx[1]

  num_parts <- strsplit(lines[num_idx], "\t")[[1]]
  den_parts <- strsplit(lines[den_idx], "\t")[[1]]
  true_parts <- strsplit(lines[true_idx], "\t")[[1]]

  max_len <- max(length(num_parts), length(den_parts), length(true_parts))
  num_parts  <- c(num_parts,  rep("", max_len - length(num_parts)))
  den_parts  <- c(den_parts,  rep("", max_len - length(den_parts)))
  true_parts <- c(true_parts, rep("", max_len - length(true_parts)))

  data_col_indices <- c()
  for (j in 2:max_len) {
    nv <- trimws(num_parts[j])
    if (nzchar(nv) && nv != "|") data_col_indices <- c(data_col_indices, j)
  }
  if (length(data_col_indices) < 2) return(NULL)

  groups <- list()
  current <- c(data_col_indices[1])
  for (k in 2:length(data_col_indices)) {
    if (data_col_indices[k] > data_col_indices[k - 1] + 1) {
      groups <- c(groups, list(current))
      current <- c()
    }
    current <- c(current, data_col_indices[k])
  }
  groups <- c(groups, list(current))

  data_start <- true_idx + 1
  if (data_start > length(lines)) return(NULL)
  data_lines <- lines[data_start:length(lines)]
  data_lines <- data_lines[nzchar(trimws(data_lines))]
  if (length(data_lines) == 0) return(NULL)

  data_split <- strsplit(data_lines, "\t")
  n_rows <- length(data_split)
  data_mat <- matrix(NA_real_, nrow = n_rows, ncol = max_len)
  for (i in seq_along(data_split)) {
    vals <- data_split[[i]]
    for (j in seq_along(vals)) {
      if (j <= max_len) {
        v <- trimws(vals[j])
        if (nzchar(v) && v != "|") suppressWarnings(data_mat[i, j] <- as.numeric(v))
      }
    }
  }

  results <- list()
  for (g_idx in seq_along(groups)) {
    cols <- groups[[g_idx]]
    if (length(cols) < 2) next

    related_cols <- c()
    unrelated_cols <- c()
    for (col in cols) {
      true_val <- trimws(true_parts[col])
      num_val  <- trimws(num_parts[col])
      if (nzchar(true_val) && true_val == num_val) {
        related_cols <- c(related_cols, col)
      } else {
        unrelated_cols <- c(unrelated_cols, col)
      }
    }

    if (length(related_cols) > 0 && length(unrelated_cols) > 0) {
      rel_data   <- as.vector(data_mat[, related_cols, drop = FALSE])
      unrel_data <- as.vector(data_mat[, unrelated_cols, drop = FALSE])
      df <- data.frame(Unrelated = unrel_data, Related = rel_data)
      df <- df[complete.cases(df), ]
      if (nrow(df) > 0) {
        num_name <- trimws(num_parts[cols[1]])
        den_name <- trimws(den_parts[cols[1]])
        comp_name <- paste0(num_name, " vs ", den_name)
        results[[comp_name]] <- df
      }
    }
  }

  if (length(results) == 0) return(NULL)
  results
}

.try_familias <- function(lines) {
  if (length(lines) < 3) return(NULL)

  sep <- "\t"
  row2_parts_tab <- strsplit(lines[2], "\t")[[1]]
  row2_parts_comma <- strsplit(lines[2], ",")[[1]]
  if (length(row2_parts_comma) > length(row2_parts_tab) && length(row2_parts_tab) <= 2) {
    sep <- ","
  }

  row2_parts <- strsplit(lines[2], sep)[[1]]
  row2_clean <- trimws(row2_parts)
  row2_non_empty <- row2_clean[nzchar(row2_clean)]

  if (length(row2_non_empty) < 2) return(NULL)
  suppressWarnings(nums <- as.numeric(row2_non_empty))
  if (!all(is.na(nums))) return(NULL)

  row3_parts <- strsplit(lines[3], sep)[[1]]
  row3_ne <- trimws(row3_parts)
  row3_ne <- row3_ne[nzchar(row3_ne)]
  suppressWarnings(row3_nums <- as.numeric(row3_ne))
  if (length(row3_nums) < 2 || all(is.na(row3_nums))) return(NULL)

  row1_parts <- trimws(strsplit(lines[1], sep)[[1]])

  n_cols <- max(length(row1_parts), length(row2_clean))
  row1_parts <- c(row1_parts, rep("", n_cols - length(row1_parts)))
  row2_clean <- c(row2_clean, rep("", n_cols - length(row2_clean)))

  data_lines <- lines[3:length(lines)]
  data_lines <- data_lines[nzchar(trimws(data_lines))]
  if (length(data_lines) == 0) return(NULL)

  data_split <- strsplit(data_lines, sep)
  n_rows <- length(data_split)
  max_data_cols <- max(vapply(data_split, length, integer(1)), na.rm = TRUE)
  n_cols <- max(n_cols, max_data_cols)
  row1_parts <- c(row1_parts, rep("", n_cols - length(row1_parts)))
  row2_clean <- c(row2_clean, rep("", n_cols - length(row2_clean)))

  data_mat <- matrix(NA_real_, nrow = n_rows, ncol = n_cols)
  for (i in seq_along(data_split)) {
    vals <- trimws(data_split[[i]])
    for (j in seq_along(vals)) {
      if (j <= n_cols) suppressWarnings(data_mat[i, j] <- as.numeric(vals[j]))
    }
  }

  valid_cols <- which(nzchar(row1_parts) | nzchar(row2_clean))
  if (length(valid_cols) < 2) return(NULL)

  has_ref <- any(grepl("^Reference pedigree$", row2_clean, ignore.case = TRUE))

  if (has_ref) {
    families <- list()
    for (col in valid_cols) {
      fam_name <- row1_parts[col]
      col_type <- row2_clean[col]
      if (!nzchar(fam_name)) next
      if (is.null(families[[fam_name]]))
        families[[fam_name]] <- list(ref_col = NULL, rel_col = NULL, rel_label = NULL)

      if (grepl("^Reference pedigree$", col_type, ignore.case = TRUE)) {
        families[[fam_name]]$ref_col <- col
      } else if (nzchar(col_type)) {
        families[[fam_name]]$rel_col <- col
        families[[fam_name]]$rel_label <- col_type
      }
    }

    results <- list()
    for (fam_name in names(families)) {
      fam <- families[[fam_name]]
      if (is.null(fam$ref_col) || is.null(fam$rel_col)) next
      df <- data.frame(
        Unrelated = data_mat[, fam$ref_col],
        Related   = data_mat[, fam$rel_col]
      )
      df <- df[complete.cases(df), ]
      if (nrow(df) > 0) {
        label <- paste0(fam_name, " (", fam$rel_label, ")")
        results[[label]] <- df
      }
    }
    if (length(results) > 0) return(results)
  }

  df <- data.frame(
    Unrelated = data_mat[, valid_cols[1]],
    Related   = data_mat[, valid_cols[2]]
  )
  df <- df[complete.cases(df), ]
  if (nrow(df) == 0) return(NULL)

  fam_name  <- row1_parts[valid_cols[1]]
  rel_label <- row2_clean[valid_cols[2]]
  if (nzchar(fam_name) && nzchar(rel_label)) {
    label <- paste0(fam_name, " (", rel_label, ")")
  } else {
    label <- "Default"
  }
  setNames(list(df), label)
}

.try_csv_generic <- function(filepath, lines) {
  first <- lines[1]
  sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else if (grepl(";", first, fixed = TRUE)) ";" else ","
  dec <- if (sep == ";") "," else "."
  header_fields <- trimws(gsub("\"", "", strsplit(first, sep, fixed = TRUE)[[1]]))
  has_header <- any(is.na(suppressWarnings(as.numeric(sub(",", ".", header_fields, fixed = TRUE)))))

  tryCatch({
    data <- utils::read.csv(filepath, sep = sep, dec = dec, header = has_header,
                            stringsAsFactors = FALSE)
    if (ncol(data) < 2) return(NULL)

    names_lower <- tolower(trimws(names(data)))
    unrel_col <- NULL; rel_col <- NULL

    unrel_patterns <- c("unrelated", "h2", "rped", "random", "reference",
                        "x1", "lr_unrelated", "lr.unrelated", "false")
    for (pat in unrel_patterns) {
      idx <- which(names_lower == pat)
      if (length(idx) > 0) { unrel_col <- idx[1]; break }
    }

    rel_patterns <- c("related", "h1", "tped", "true_ped", "match",
                      "x2", "lr_related", "lr.related")
    for (pat in rel_patterns) {
      idx <- which(names_lower == pat)
      if (length(idx) > 0) { rel_col <- idx[1]; break }
    }

    if (is.null(unrel_col) || is.null(rel_col)) {
      num_cols <- which(vapply(data, function(x) {
        is.numeric(x) || !all(is.na(suppressWarnings(as.numeric(x))))
      }, logical(1)))
      if (length(num_cols) >= 2) {
        if (is.null(unrel_col)) unrel_col <- num_cols[1]
        if (is.null(rel_col))   rel_col   <- num_cols[2]
      }
    }

    if (is.null(unrel_col)) unrel_col <- 1
    if (is.null(rel_col))   rel_col   <- 2

    df <- data.frame(
      Unrelated = suppressWarnings(as.numeric(data[[unrel_col]])),
      Related   = suppressWarnings(as.numeric(data[[rel_col]]))
    )
    df <- df[complete.cases(df), ]
    if (nrow(df) > 0) list("Default" = df) else NULL
  }, error = function(e) NULL)
}

genetic_examples <- c(
  "Full sibling" = "full_sibling.txt",
  "Half sibling" = "half_sibling.txt",
  "Grandparent" = "grandparent.txt",
  "Uncle/Aunt" = "uncle.txt",
  "Two uncles/aunts" = "two_uncles.txt",
  "First cousin" = "first_cousin.txt",
  "Two first cousins" = "two_first_cousins.txt",
  "Uncle/aunt and first cousin" = "uncle_first_cousin.txt",
  "Family tree (Familias example)" = "family_tree.txt"
)

genetic_example_help <- function(file) {
  if (identical(file, "family_tree.txt")) {
    "Familias example dataset, 10,000 simulations."
  } else {
    paste("Simulated with forrel 1.9.1: 10,000 replicates under H1 (the person of interest",
          "is the missing person) and H2 (unrelated), 24 autosomal STRs, Argentina allele",
          "frequencies, no mutation model. References are genotyped relatives of the missing person.")
  }
}

candidates_example_path <- file.path("examples", "nongenetic", "candidates_example.csv")

is_binary_file <- function(filepath) {
  head_bytes <- readBin(filepath, "raw", n = 4096)
  length(head_bytes) < 2 || any(head_bytes == as.raw(0)) ||
    identical(head_bytes[1:2], charToRaw("PK"))
}

row_list <- function(rows) {
  paste0(paste(head(rows, 5) + 1, collapse = ", "), if (length(rows) > 5) ", ..." else "")
}

read_candidates <- function(filepath) {
  if (is.null(filepath) || !file.exists(filepath)) return("The file could not be read.")
  if (is_binary_file(filepath)) {
    return("The file is not a text CSV file. Excel and other binary files are not supported.")
  }
  first <- tryCatch(readLines(filepath, n = 1, warn = FALSE), error = function(e) character(0))
  if (length(first) == 0) return("The file is empty.")
  if (!validUTF8(first)) first <- iconv(first, "latin1", "UTF-8")
  sep <- if (grepl("\t", first, fixed = TRUE)) "\t" else if (grepl(";", first, fixed = TRUE)) ";" else ","
  data <- tryCatch(
    utils::read.table(filepath, sep = sep, header = TRUE, colClasses = "character",
                      quote = "\"", comment.char = "", strip.white = TRUE,
                      check.names = FALSE, fileEncoding = if (validUTF8(first)) "" else "latin1"),
    error = function(e) conditionMessage(e))
  if (is.character(data)) return(paste("Could not read the file:", data))
  names(data) <- trimws(names(data))
  required <- c("ID", "Sex", "Age", "HairColor")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    return(sprintf("Missing required column(s): %s. Required columns: ID, Sex, Age, HairColor.",
                   paste(missing_cols, collapse = ", ")))
  }
  if (nrow(data) == 0) return("The file has no candidate rows.")
  data <- data[, required]
  data$Sex <- toupper(data$Sex)
  bad <- which(!nzchar(data$ID))
  if (length(bad) > 0) return(paste("Column ID is empty in line(s)", row_list(bad)))
  bad <- which(duplicated(data$ID))
  if (length(bad) > 0) return(paste("Column ID has duplicated values in line(s)", row_list(bad)))
  bad <- which(!data$Sex %in% c("F", "M"))
  if (length(bad) > 0) return(paste("Column Sex must be F or M; invalid value in line(s)", row_list(bad)))
  age <- suppressWarnings(as.numeric(data$Age))
  bad <- which(is.na(age) | age < 1 | age > 80 | age != round(age))
  if (length(bad) > 0) {
    return(paste("Column Age must be a whole number from 1 to 80; invalid value in line(s)", row_list(bad)))
  }
  hair <- match(tolower(data$HairColor), tolower(names(DEFAULT_HAIR_FREQ)))
  bad <- which(is.na(hair))
  if (length(bad) > 0) {
    return(paste("Column HairColor must be one of Black, Brown, Blonde, Red, Gray; invalid value in line(s)",
                 row_list(bad)))
  }
  data$Age <- age
  data$HairColor <- names(DEFAULT_HAIR_FREQ)[hair]
  data
}

candidate_lrs <- function(cands, MPs, MPa, MPr, MPc, eps, epa, epc, empirical) {
  aw <- age_window(MPa, MPr)
  propF <- 0.51
  hair_freq <- DEFAULT_HAIR_FREQ
  T1p <- aw$T1p
  notes <- character(0)
  if (empirical) {
    propF <- mean(cands$Sex == "F")
    hair_freq <- vapply(names(DEFAULT_HAIR_FREQ), function(h) mean(cands$HairColor == h), numeric(1))
    T1p <- mean(cands$Age >= aw$min & cands$Age <= aw$max)
    if (propF %in% c(0, 1)) {
      propF <- 0.51
      notes <- c(notes, "Only one sex is present in the database: the default female proportion (0.51) is used.")
    }
    for (h in names(hair_freq)[hair_freq == 0]) {
      hair_freq[h] <- DEFAULT_HAIR_FREQ[h]
      notes <- c(notes, sprintf("No candidate has hair color %s: the default population frequency (%.2f) is used for it.",
                                h, DEFAULT_HAIR_FREQ[h]))
    }
    if (T1p %in% c(0, 1)) {
      T1p <- aw$T1p
      notes <- c(notes, sprintf("The share of candidates inside the age range is %s: the uniform 1-80 age model is used.",
                                if (mean(cands$Age >= aw$min & cands$Age <= aw$max) == 0) "0" else "1"))
    }
  }
  hair_idx <- match(cands$HairColor, names(DEFAULT_HAIR_FREQ))
  LR <- vapply(seq_len(nrow(cands)), function(i) {
    calc_lr_sex(MPs, cands$Sex[i], eps, propF) *
      calc_lr_age(MPa, MPr, cands$Age[i], epa, T1p) *
      calc_lr_hair(MPc, hair_idx[i], epc, hair_freq)
  }, numeric(1))
  res <- data.frame(cands, LR = LR, log10LR = log10(LR), stringsAsFactors = FALSE)
  res <- res[order(-res$LR, res$ID), ]
  rownames(res) <- NULL
  list(results = res, propF = propF, hair_freq = hair_freq, T1p = T1p, age_min = aw$min,
       age_max = aw$max, notes = notes)
}

frac_below <- function(x, t) {
  findInterval(t, sort(x), left.open = TRUE) / length(x)
}

compute_decision_curve <- function(sim_df, n_points = 1000) {
  max_lr <- max(1, min(max(c(sim_df$Unrelated, sim_df$Related)), 10000))
  lr_values <- 10^seq(0, log10(max_lr), length.out = n_points)
  data.frame(LR = lr_values,
             FPR = 1 - frac_below(sim_df$Unrelated, lr_values),
             FNR = frac_below(sim_df$Related, lr_values))
}

compute_optimal_dt <- function(sim_df, weight = 10) {
  max_lr <- max(1, min(max(c(sim_df$Unrelated, sim_df$Related)), 10000))
  observed <- c(sim_df$Unrelated, sim_df$Related)
  lr_values <- sort(unique(c(1, observed[observed > 1 & observed <= max_lr], max_lr)))
  FPR <- 1 - frac_below(sim_df$Unrelated, lr_values)
  FNR <- frac_below(sim_df$Related, lr_values)
  idx <- which.min(sqrt(FNR^2 + (weight * FPR)^2))
  list(DT = lr_values[idx], FPR = FPR[idx], FNR = FNR[idx])
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
        font-size: 14px;
        color: #222;
        background-color: #ffffff;
      }

      h1, h2, h3, h4, h5, h6 {
        font-weight: 600;
        color: #222;
      }

      h2 { font-size: 22px; }
      h3 { font-size: 19px; }
      h4 { font-size: 16px; }
      h5 { font-size: 14px; }

      a { color: #1f3a5f; }

      .main-header {
        background: #1f3a5f;
        color: #ffffff;
        padding: 10px 20px;
        margin: 0 -15px 15px -15px;
        display: flex;
        align-items: center;
      }

      .main-header img {
        height: 40px;
        margin-right: 14px;
        background: #ffffff;
        border-radius: 2px;
      }

      .main-header h1 {
        margin: 0;
        font-size: 20px;
        font-weight: 600;
        color: #ffffff;
      }

      .main-header p {
        margin: 2px 0 0 0;
        font-size: 13px;
        font-weight: normal;
      }

      .info-box, .warning-box, .success-box, .tutorial-section {
        background: #ffffff;
        border: 1px solid #ddd;
        padding: 12px 15px;
        margin: 12px 0;
        border-radius: 2px;
      }

      .warning-box {
        background: #fafafa;
      }

      .result-box {
        background: #ffffff;
        border: 1px solid #ddd;
        padding: 15px;
        margin: 12px 0;
        border-radius: 2px;
        text-align: center;
      }

      .result-box h2 {
        margin: 6px 0;
        font-size: 24px;
        font-weight: 600;
        color: #222;
      }

      .result-box h4 {
        margin: 0;
        font-weight: normal;
        color: #555;
      }

      .formula-box {
        background: #f7f7f7;
        color: #222;
        border: 1px solid #ddd;
        padding: 12px 15px;
        border-radius: 2px;
        margin: 12px 0;
      }

      .formula-box pre, .shiny-text-output, pre {
        border-radius: 2px;
      }

      .well {
        background: #ffffff;
        border: 1px solid #ddd;
        border-radius: 2px;
        box-shadow: none;
      }

      .interpretation-text {
        font-size: 15px;
        line-height: 1.5;
        color: #222;
      }

      .metric-card {
        background: #ffffff;
        border: 1px solid #ddd;
        border-radius: 2px;
        padding: 10px;
        text-align: center;
        margin-bottom: 12px;
      }

      .metric-card .metric-value {
        font-size: 22px;
        font-weight: 600;
        color: #222;
        margin: 6px 0;
      }

      .metric-card .metric-label {
        font-size: 13px;
        color: #555;
      }

      .upload-zone {
        padding: 0;
      }

      .nav-tabs {
        border-bottom: 1px solid #ddd;
      }

      .nav-tabs > li > a, .nav-tabs > li > a:hover, .nav-tabs > li > a:focus {
        color: #444;
        border-radius: 2px 2px 0 0;
        background: transparent;
      }

      .nav-tabs > li.active > a, .nav-tabs > li.active > a:hover, .nav-tabs > li.active > a:focus {
        color: #1f3a5f;
        font-weight: 600;
        border: 1px solid #ddd;
        border-bottom-color: #ffffff;
      }

      .nav-pills {
        border-bottom: 1px solid #ddd;
        margin-bottom: 8px;
      }

      .nav-pills > li > a, .nav-pills > li > a:hover, .nav-pills > li > a:focus {
        color: #444;
        background: transparent;
        border-radius: 0;
        border-bottom: 2px solid transparent;
        padding: 6px 12px;
      }

      .nav-pills > li.active > a, .nav-pills > li.active > a:hover, .nav-pills > li.active > a:focus {
        color: #1f3a5f;
        background: transparent;
        font-weight: 600;
        border-bottom: 2px solid #1f3a5f;
      }

      .btn {
        border-radius: 2px;
        box-shadow: none;
        transition: none;
      }

      .btn-primary, .btn-primary:hover, .btn-primary:focus, .btn-primary:active {
        background: #1f3a5f;
        border-color: #1f3a5f;
        color: #ffffff;
      }

      .btn-default {
        background: #ffffff;
        border-color: #bbb;
        color: #222;
      }

      .form-control, .selectize-input, .selectize-dropdown {
        border-radius: 2px;
        box-shadow: none;
        border-color: #ccc;
      }

      .irs-bar, .irs-bar-edge, .irs-single, .irs-from, .irs-to,
      .irs--shiny .irs-bar, .irs--shiny .irs-bar-edge, .irs--shiny .irs-single,
      .irs--shiny .irs-from, .irs--shiny .irs-to {
        background: #1f3a5f;
        border-color: #1f3a5f;
      }

      .table {
        color: #222;
      }

      .progress {
        border-radius: 2px;
        box-shadow: none;
      }
    "))
  ),

  div(class = "main-header",
    tags$img(src = "MispiIcon.png", alt = "mispitools logo"),
    div(
      h1("mispitools"),
      p("Likelihood Ratio Analysis for Kinship in Forensic Science")
    )
  ),

  tabsetPanel(
    id = "main_tabs",
    type = "tabs",

    tabPanel(
      "Overview",
      fluidRow(
        column(10, offset = 1,
          br(),
          h2("Evidence Evaluation for Kinship in Forensic Science"),
          p(class = "interpretation-text",
            "This application implements a Bayesian framework for evaluating
            evidence in kinship cases for forensic science. It includes two main sections:
            non-genetic evidence (sex, age, hair color) with simulation capabilities,
            and genetic evidence analysis (upload pre-computed simulation results)."
          ),
          hr(),

          fluidRow(
            column(6,
              div(class = "info-box",
                h4("The Likelihood Ratio"),
                p("The LR compares two hypotheses:"),
                tags$ul(
                  tags$li(strong("H1:"), " The individuals are related as stated in the pedigree"),
                  tags$li(strong("H2:"), " The individuals are unrelated")
                ),
                div(class = "formula-box",
                  "LR = P(Evidence | H1) / P(Evidence | H2)"
                ),
                p("The LR does ", strong("not"), " provide a probability of identity.
                         It quantifies how much the evidence changes the odds.")
              )
            ),
            column(6,
              div(class = "info-box",
                h4("Using the LR"),
                p("The LR quantifies how many times more likely the evidence is under H1 vs H2."),
                tags$ul(
                  tags$li("LR > 1: Evidence more probable under H1"),
                  tags$li("LR < 1: Evidence more probable under H2"),
                  tags$li("LR = 1: Evidence equally probable under both hypotheses")
                ),
                p("The numeric LR value should be reported without verbal qualifiers.")
              )
            )
          ),

          h4("Application Sections"),
          fluidRow(
            column(6,
              div(class = "well", style = "text-align: center; padding: 20px;",
                h5("Non-Genetic Evidence"),
                p(class = "text-muted", "Calculate and simulate LRs from non-genetic evidence:
                  biological sex, age, and hair color. Includes LR simulation,
                  evidence combination, and decision analysis.")
              )
            ),
            column(6,
              div(class = "well", style = "text-align: center; padding: 20px;",
                h5("Genetic Evidence"),
                p(class = "text-muted", "Analyze pre-computed genetic simulation results
                  from Familias, forrel, mispitools (R), or FamLink.
                  Evaluate statistical power and optimize LR decision thresholds.")
              )
            )
          ),

          hr(),
          div(style = "color: #444; padding: 10px 0;",
            p(strong("References")),
            p("Egeland T., Marsico F.L. (2026). Using all available evidence to solve kinship cases.
              International Journal of Legal Medicine 140(4):1981-1996. ",
              tags$a(href = "https://doi.org/10.1007/s00414-026-03773-6", target = "_blank",
                     "https://doi.org/10.1007/s00414-026-03773-6")),
            p("Marsico F.L., Caridi I. (2023). Incorporating non-genetic evidence in large scale
              missing person searches: A general approach beyond filtering.
              Forensic Science International: Genetics 66:102891"),
            p("Marsico F.L., Vigeland M.D., et al. (2021). Forensic Science International:
              Genetics 54:102519")
          )
        )
      )
    ),

    tabPanel(
      "Non-Genetic Evidence",
      br(),
      fluidRow(column(10, offset = 1,
        p(class = "interpretation-text",
          "Calculate and simulate likelihood ratios from non-genetic evidence types
          (sex, age, hair color). All analyses in this section can be run directly
          in the web application.")
      )),
      tabsetPanel(
        id = "nongen_tabs",
        type = "pills",

    tabPanel(
      "Individual Evidence",
      fluidRow(
        column(12,
          br(),
          h3("Single Evidence Type Analysis"),
          p("Calculate the likelihood ratio for one type of evidence.")
        )
      ),
      fluidRow(
        column(3,
          wellPanel(
            selectInput("indiv_type", "Evidence Type:",
              choices = c("Biological Sex" = "sex",
                          "Age" = "age",
                          "Hair Color" = "hair")),
            hr(),

            conditionalPanel(
              condition = "input.indiv_type == 'sex'",
              h5("Case Parameters"),
              selectInput("indiv_MPs", "MP reported sex:", sex_choices),
              selectInput("indiv_obsSex", "Observed sex:", sex_choices),
              hr(),
              h5("Model Parameters"),
              sliderInput("indiv_eps_sex", "Observation error rate:",
                0, 0.3, 0.05, 0.01),
              sliderInput("indiv_propF", "Female proportion in population:",
                0.3, 0.7, 0.51, 0.01)
            ),

            conditionalPanel(
              condition = "input.indiv_type == 'age'",
              h5("Case Parameters"),
              numericInput("indiv_MPa", "MP estimated age (years):", 35, 1, 80),
              numericInput("indiv_MPr", "Estimation uncertainty (years):", 5, 1, 20),
              numericInput("indiv_obsAge", "Observed age:", 37, 1, 80),
              hr(),
              h5("Model Parameters"),
              sliderInput("indiv_epa", "Age category error rate:",
                0, 0.3, 0.05, 0.01)
            ),

            conditionalPanel(
              condition = "input.indiv_type == 'hair'",
              h5("Case Parameters"),
              selectInput("indiv_MPc", "MP hair color:", hair_choices, selected = 2),
              selectInput("indiv_obsC", "Observed hair color:", hair_choices, selected = 2),
              hr(),
              h5("Model Parameters"),
              sliderInput("indiv_epc", "Color misclassification rate (total):",
                0, 0.3, 0.05, 0.01),
              p(class = "text-muted small", hair_error_help),
              p(class = "text-muted small",
                "Black: 8%, Brown: 45%, Blonde: 30%, Red: 7%, Gray: 10%")
            )
          )
        ),
        column(9,
          div(class = "result-box",
            h4("Likelihood Ratio"),
            h2(textOutput("indiv_lr_value")),
            p(textOutput("indiv_lr_log"))
          ),
          div(class = "formula-box",
            h5(style = "margin-top: 0;", "Calculation Details"),
            verbatimTextOutput("indiv_formula")
          )
        )
      )
    ),

    tabPanel(
      "Probability Tables",
      fluidRow(
        column(12,
          br(),
          h3("Conditional Probability Tables"),
          p("Visualization of evidence probabilities under each hypothesis.")
        )
      ),
      fluidRow(
        column(2,
          wellPanel(
            h5("Missing Person"),
            selectInput("cpt_MPs", "Sex:", sex_choices),
            numericInput("cpt_MPa", "Age:", 40, 1, 80),
            numericInput("cpt_MPr", "Age range:", 6, 1, 20),
            selectInput("cpt_MPc", "Hair color:", hair_choices, selected = 2)
          )
        ),
        column(2,
          wellPanel(
            h5("Error Rates"),
            numericInput("cpt_eps", "Sex:", 0.05, 0, 0.5, 0.01),
            numericInput("cpt_epa", "Age:", 0.05, 0, 0.5, 0.01),
            numericInput("cpt_epc", "Hair (total):", 0.05, 0, 0.5, 0.01),
            p(class = "text-muted small", hair_error_help)
          ),
          wellPanel(
            h5("Population"),
            sliderInput("cpt_propF", "Female %:", 0.3, 0.7, 0.51, 0.01)
          )
        ),
        column(8,
          plotOutput("cpt_plots", height = "auto"),
          div(class = "info-box",
            fluidRow(
              column(4,
                strong("P(D|H2):"), " Phenotype probability in reference population"
              ),
              column(4,
                strong("P(D|H1):"), " Probability if person IS the MP"
              ),
              column(4,
                strong("log10(LR):"),
                span(" Blue = supports H1,"),
                span(" red = supports H2, grey = LR = 0")
              )
            )
          )
        )
      )
    ),

    tabPanel(
      "LR Distribution",
      fluidRow(
        column(12,
          br(),
          h3("Likelihood Ratio Distributions"),
          p("Comparison of LR distributions under H1 and H2.")
        )
      ),
      fluidRow(
        column(6,
          plotOutput("dist_plot", height = "350px")
        ),
        column(6,
          plotOutput("dist_roc", height = "350px")
        )
      ),
      fluidRow(
        column(6,
          div(class = "success-box",
            h5("Under H1 (Person IS the MP)"),
            verbatimTextOutput("dist_stats_h1")
          )
        ),
        column(6,
          div(class = "warning-box",
            h5("Under H2 (Person is NOT the MP)"),
            verbatimTextOutput("dist_stats_h2")
          )
        )
      )
    ),

    tabPanel(
      "Combine Evidence",
      fluidRow(
        column(12,
          br(),
          h3("Evidence Combination"),
          p("Under conditional independence, individual LRs multiply to give a combined LR.")
        )
      ),
      fluidRow(
        column(12,
          div(class = "warning-box",
            strong("Independence Assumption: "),
            "LR multiplication assumes conditional independence of evidence sources."
          )
        )
      ),
      fluidRow(
        column(3,
          wellPanel(
            h5("Sex Evidence"),
            checkboxInput("comb_use_sex", "Include", TRUE),
            conditionalPanel(
              condition = "input.comb_use_sex",
              selectInput("comb_MPs", "MP:", sex_choices),
              selectInput("comb_obsSex", "Observed:", sex_choices),
              numericInput("comb_eps", "Error:", 0.05, 0, 0.3, 0.01),
              sliderInput("comb_propF", "Female proportion in population:",
                0.3, 0.7, 0.51, 0.01)
            ),
            hr(),
            strong(textOutput("comb_lr_sex"))
          )
        ),
        column(3,
          wellPanel(
            h5("Age Evidence"),
            checkboxInput("comb_use_age", "Include", TRUE),
            conditionalPanel(
              condition = "input.comb_use_age",
              numericInput("comb_MPa", "MP age:", 35, 1, 80),
              numericInput("comb_MPr", "Range:", 5, 1, 20),
              numericInput("comb_obsAge", "Observed:", 37, 1, 80),
              numericInput("comb_epa", "Error:", 0.05, 0, 0.3, 0.01)
            ),
            hr(),
            strong(textOutput("comb_lr_age"))
          )
        ),
        column(3,
          wellPanel(
            h5("Hair Color Evidence"),
            checkboxInput("comb_use_hair", "Include", TRUE),
            conditionalPanel(
              condition = "input.comb_use_hair",
              selectInput("comb_MPc", "MP:", hair_choices, selected = 2),
              selectInput("comb_obsC", "Observed:", hair_choices, selected = 2),
              numericInput("comb_epc", "Error (total):", 0.05, 0, 0.3, 0.01),
              p(class = "text-muted small", hair_error_help)
            ),
            hr(),
            strong(textOutput("comb_lr_hair"))
          )
        ),
        column(3,
          div(class = "result-box",
            h4("Combined Result"),
            h2(textOutput("comb_lr_total")),
            p(textOutput("comb_lr_log")),
            hr(),
            verbatimTextOutput("comb_formula")
          )
        )
      )
    ),

    tabPanel(
      "Database Search",
      fluidRow(
        column(12,
          br(),
          h3("Database Search"),
          p("Compute the non-genetic LR of every candidate in a database against the missing
            person profile, and count candidates at or above a threshold (expected false
            positives in a database search, Marsico & Caridi 2023).")
        )
      ),
      fluidRow(
        column(3,
          wellPanel(
            h5("Candidate Database"),
            fileInput("dbs_file", NULL, accept = c(".csv", ".tsv", ".txt", "text/csv", "text/plain")),
            p(class = "text-muted small",
              "CSV with columns ID, Sex (F/M), Age (whole number 1-80) and HairColor
              (Black, Brown, Blonde, Red, Gray). Comma, semicolon or tab separated."),
            actionButton("dbs_load_example", "Load example", class = "btn-primary btn-block"),
            downloadButton("dbs_download_example", "Download example",
              class = "btn-default btn-block", style = "margin-top: 6px;"),
            p(class = "text-muted small", style = "margin-top: 6px;",
              "Synthetic example: 2,001 candidates. C0001 matches the default MP profile, but every candidate in the same sex, age-range and hair-colour class gets the same LR: non-genetic evidence ranks classes, it does not identify individuals."),
            uiOutput("dbs_msg"),
            hr(),
            h5("Missing Person"),
            selectInput("dbs_MPs", "Sex:", sex_choices),
            numericInput("dbs_MPa", "Age:", 40, 1, 80),
            numericInput("dbs_MPr", "Age range:", 6, 1, 20),
            selectInput("dbs_MPc", "Hair color:", hair_choices, selected = 2),
            hr(),
            h5("Error Rates"),
            numericInput("dbs_eps", "Sex:", 0.05, 0, 0.5, 0.01),
            numericInput("dbs_epa", "Age:", 0.05, 0, 0.5, 0.01),
            numericInput("dbs_epc", "Hair (total):", 0.05, 0, 0.5, 0.01),
            p(class = "text-muted small", hair_error_help),
            hr(),
            checkboxInput("dbs_empirical", "Estimate population frequencies from the database", TRUE),
            p(class = "text-muted small",
              "Unchecked: female proportion 0.51, default hair color frequencies and the
              uniform 1-80 age model.")
          )
        ),
        column(9,
          wellPanel(
            sliderInput("dbs_threshold", HTML("log<sub>10</sub>(LR) threshold:"),
              min = -3, max = 3, value = 0, step = 0.1),
            verbatimTextOutput("dbs_summary")
          ),
          plotOutput("dbs_hist", height = "350px"),
          fluidRow(
            column(6, numericInput("dbs_top", "Top candidates shown:", 20, 1, 500)),
            column(6, br(), downloadButton("dbs_download_results", "Download results"))
          ),
          div(style = "overflow-x: auto;", tableOutput("dbs_table"))
        )
      )
    ),

    tabPanel(
      "Decision Analysis",
      fluidRow(
        column(12,
          br(),
          h3("Decision Threshold Selection"),
          p("Evaluate the trade-off between false positive and false negative rates.")
        )
      ),
      fluidRow(
        column(4,
          wellPanel(
            h5("Threshold Selection"),
            sliderInput("dec_threshold", HTML("log<sub>10</sub>(LR) Threshold:"),
                               min = -3, max = 3, value = 0, step = 0.1),
            hr(),
            h5("Performance Metrics"),
            verbatimTextOutput("dec_metrics")
          ),
          div(class = "info-box",
            h6("Metric Definitions"),
            tags$ul(class = "small",
              tags$li(strong("TPR:"), " True Positive Rate (sensitivity)"),
              tags$li(strong("FPR:"), " False Positive Rate (1 - specificity)"),
              tags$li(strong("Rule:"), " a person is declared a match when log10(LR) >= threshold"),
              tags$li(strong("MCC:"), " Matthews Correlation Coefficient (assumes balanced classes)")
            )
          )
        ),
        column(8,
          plotOutput("dec_roc_plot", height = "400px")
        )
      )
    )
    )
    ),

    tabPanel(
      "Genetic Evidence",
      br(),
      fluidRow(column(10, offset = 1,
        div(class = "info-box",
          h4("About Genetic Simulations"),
          p("The web version of mispitools cannot perform genetic simulations
            due to computational requirements. Simulations must be run using
            dedicated software and then uploaded here for analysis:"),
          tags$ul(
            tags$li(strong("Familias"), " - standalone software for forensic genetics"),
            tags$li(strong("forrel"), " - R package for forensic pedigree analysis"),
            tags$li(strong("mispitools"), " - R package (command-line version)"),
            tags$li(strong("FamLink"), " - software for kinship analysis")
          ),
          p("Once you have simulation results, upload them below to analyze
            statistical power, LR distributions, and optimize decision thresholds.")
        )
      )),
      tabsetPanel(
        id = "gen_tabs",
        type = "pills",

    tabPanel(
      "Simulation Analysis",
      fluidRow(
        column(12,
          br(),
          h3("Simulation Results Analysis"),
          p(class = "interpretation-text",
            "Upload simulation data from forensic genetics software (Familias, Forrel, FamLink, mispitools)
            to evaluate statistical power and optimize LR decision thresholds."
          )
        )
      ),
      fluidRow(
        column(4,
          wellPanel(
            div(class = "upload-zone",
              h5("Upload Simulation Data"),
              fileInput("sim_file", NULL,
                accept = c(".csv", ".tsv", ".txt",
                           "text/csv", "text/tab-separated-values", "text/plain")),
              p(class = "text-muted small",
                "Accepted: Familias, FamLink, forrel, mispitools, or any CSV/TSV
                 with Unrelated & Related columns."),
              hr(),
              selectInput("sim_example", "Example dataset", genetic_examples),
              actionButton("sim_load_example", "Load example",
                class = "btn-primary btn-block"),
              downloadButton("sim_download_example", "Download example",
                class = "btn-default btn-block", style = "margin-top: 6px;"),
              p(class = "text-muted small", style = "margin-top: 6px;",
                textOutput("sim_example_help", inline = TRUE)),
              uiOutput("sim_example_pedigree")
            ),
            uiOutput("sim_parse_msg"),
            uiOutput("sim_family_ui"),
            hr(),
            conditionalPanel(
              condition = "output.sim_has_data",
              div(class = "success-box", style = "padding: 10px 15px;",
                textOutput("sim_data_info")
              ),
              hr(),
              h5("LR Threshold"),
              numericInput("sim_threshold", "LR threshold value (match if LR >= threshold):",
                value = 100, min = 0, max = 1000000, step = 1),
              hr(),
              h5("Optimal Decision Threshold"),
              sliderInput("sim_weight", "FP/FN cost weight:",
                min = 1, max = 50, value = 10, step = 1),
              p(class = "text-muted small",
                "Weight > 1 penalizes false positives more. Default: 10 (FP is 10x worse than FN)."),
              actionButton("sim_calc_dt", "Compute Optimal DT",
                class = "btn-primary btn-block"),
              br(), br(),
              verbatimTextOutput("sim_dt_result")
            )
          )
        ),
        column(8,
          conditionalPanel(
            condition = "output.sim_has_data",
            div(class = "sim-subtab",
              tabsetPanel(
                id = "sim_subtabs",
                type = "pills",

                tabPanel(
                  "Summary",
                  br(),
                  fluidRow(
                    column(6,
                      div(class = "metric-card",
                        span(class = "metric-label", "Simulations"),
                        div(class = "metric-value", textOutput("sim_n_sims"))
                      )
                    ),
                    column(3,
                      div(class = "metric-card",
                        span(class = "metric-label", "Inclusion Power (LR >= 10,000)"),
                        div(class = "metric-value", textOutput("sim_ip"))
                      )
                    ),
                    column(3,
                      div(class = "metric-card",
                        span(class = "metric-label", "Exclusion Power"),
                        div(class = "metric-value", textOutput("sim_ep"))
                      )
                    )
                  ),
                  fluidRow(
                    column(6,
                      div(class = "success-box",
                        h5("Related (H1) - Summary"),
                        verbatimTextOutput("sim_summary_rel")
                      )
                    ),
                    column(6,
                      div(class = "warning-box",
                        h5("Unrelated (H2) - Summary"),
                        verbatimTextOutput("sim_summary_unrel")
                      )
                    )
                  ),
                  plotOutput("sim_violin", height = "350px")
                ),

                tabPanel(
                  "Distributions",
                  br(),
                  plotOutput("sim_density_plot", height = "400px"),
                  fluidRow(
                    column(6,
                      div(class = "success-box",
                        h5("Related (H1) - log10(LR)"),
                        verbatimTextOutput("sim_dist_stats_rel")
                      )
                    ),
                    column(6,
                      div(class = "warning-box",
                        h5("Unrelated (H2) - log10(LR)"),
                        verbatimTextOutput("sim_dist_stats_unrel")
                      )
                    )
                  )
                ),

                tabPanel(
                  "Decision Plot",
                  br(),
                  div(class = "info-box",
                    p("Each point represents a different LR threshold. The curve shows
                      the trade-off between false negative rate (missing a true match) and
                      false positive rate (incorrectly flagging a non-match).
                      Points closer to the origin (0,0) indicate better thresholds.
                      Grey labels show the LR value at selected points.")
                  ),
                  plotOutput("sim_decision_plot", height = "450px")
                ),

                tabPanel(
                  "Rates at Threshold",
                  br(),
                  div(class = "info-box",
                    p("Classification metrics at the selected LR threshold. A comparison is
                      declared a match when LR >= threshold. Expected false positives and PPV
                      assume one true match and N unrelated candidates in the search.
                      Adjust the threshold in the left panel to see how error rates change.")
                  ),
                  br(),
                  fluidRow(
                    column(3,
                      div(class = "metric-card",
                        span(class = "metric-label", "True Positive Rate"),
                        div(class = "metric-value", textOutput("sim_rate_tpr"))
                      )
                    ),
                    column(3,
                      div(class = "metric-card",
                        span(class = "metric-label", "False Positive Rate"),
                        div(class = "metric-value", textOutput("sim_rate_fpr"))
                      )
                    ),
                    column(3,
                      div(class = "metric-card",
                        span(class = "metric-label", "False Negative Rate"),
                        div(class = "metric-value", textOutput("sim_rate_fnr"))
                      )
                    ),
                    column(3,
                      div(class = "metric-card",
                        span(class = "metric-label", "True Negative Rate"),
                        div(class = "metric-value", textOutput("sim_rate_tnr"))
                      )
                    )
                  ),
                  fluidRow(
                    column(4,
                      div(class = "metric-card",
                        span(class = "metric-label", "PPV (1 true match among N candidates)"),
                        div(class = "metric-value", textOutput("sim_rate_ppv"))
                      )
                    ),
                    column(4,
                      div(class = "metric-card",
                        span(class = "metric-label", "Accuracy (balanced classes)"),
                        div(class = "metric-value", textOutput("sim_rate_acc"))
                      )
                    ),
                    column(4,
                      div(class = "metric-card",
                        span(class = "metric-label", "Expected false positives (FPR x N)"),
                        div(style = "margin-top: 8px;",
                          numericInput("sim_db_size", "N = number of unrelated candidates compared:",
                            value = 10000, min = 1, step = 100),
                          div(class = "metric-value", style = "font-size: 22px;",
                            textOutput("sim_expected_fp"))
                        )
                      )
                    )
                  ),
                  br(),
                  tableOutput("sim_rates_table")
                )
              )
            )
          ),
          conditionalPanel(
            condition = "!output.sim_has_data",
            div(style = "padding: 10px 0; color: #555;",
              h4("Upload simulation data to begin analysis")
            ),
            div(class = "info-box",
              h5("Accepted file formats:"),
              tags$ol(
                tags$li(tags$b("Familias (single family):"),
                  " Tab-separated .txt file with 2 header rows (case name + pedigree label), then numeric LR data. ",
                  "First column = Unrelated (Reference Pedigree), second = Related."),
                tags$li(tags$b("Familias (multi-family):"),
                  " Same structure but with column pairs per family. Row 1 repeats the family name, ",
                  "row 2 has 'Reference pedigree' for H2 and a relationship label (Father, Brother, etc.) for H1. ",
                  "Can contain one or more families."),
                tags$li(tags$b("FamLink:"),
                  " Starts with 'Raw data LR:s...', then 'Pedigree (Numerator)', 'Versus (Denominator)' and 'True:' header rows. ",
                  "Groups are separated by '|'. The 'True' row indicates which column is Related (H1). ",
                  "Can contain one or more pedigree comparisons."),
                tags$li(tags$b("Generic CSV/TSV:"),
                  " Any .csv, .tsv, or .txt file with columns named 'Unrelated' and 'Related', ",
                  "or two numeric columns (first = Unrelated, second = Related). ",
                  "Supports comma, tab, or semicolon separators.")
              ),
              hr(),
              h5("Example format (Familias, tab-separated):"),
              tags$pre(style = "background: #f7f7f7; padding: 10px; font-size: 12px; overflow-x: auto;",
                "FamilyName\\tFamilyName\nReference pedigree\\tFather\n3.33333\\t5\n2.5\\t3.33333\n1.25\\t2.08333"
              ),
              h5("Example format (Generic CSV):"),
              tags$pre(style = "background: #f7f7f7; padding: 10px; font-size: 12px; overflow-x: auto;",
                "Unrelated,Related\n0.5,3.2\n0.8,5.1\n0.3,2.7"
              ),
              p(class = "text-muted small",
                "Note: Files must be text-based (.txt, .csv, .tsv). Excel files (.xlsx) are not supported. ",
                "Simulations can be generated using Familias, forrel, mispitools (R version), or FamLink.")
            )
          )
        )
      )
    )
    )
    )
  )
)

server <- function(input, output, session) {

  is_narrow <- function(id) {
    w <- session$clientData[[paste0("output_", id, "_width")]]
    !is.null(w) && w < 600
  }

  plot_size <- function(id, wide) {
    if (is_narrow(id)) 10 else wide
  }

  indiv_lr <- reactive({
    if (input$indiv_type == "sex") {
      validate(
        check_num(input$indiv_eps_sex, "Observation error rate", 0, 0.3),
        check_num(input$indiv_propF, "Female proportion in population", 0.3, 0.7)
      )
      calc_lr_sex(input$indiv_MPs, input$indiv_obsSex,
                  input$indiv_eps_sex, input$indiv_propF)
    } else if (input$indiv_type == "age") {
      validate(
        check_num(input$indiv_MPa, "MP estimated age", 1, 80, whole = TRUE),
        check_num(input$indiv_MPr, "Estimation uncertainty", 1, 20, whole = TRUE),
        check_num(input$indiv_obsAge, "Observed age", 1, 80, whole = TRUE),
        check_num(input$indiv_epa, "Age category error rate", 0, 0.3)
      )
      calc_lr_age(input$indiv_MPa, input$indiv_MPr,
                  input$indiv_obsAge, input$indiv_epa)
    } else {
      validate(check_num(input$indiv_epc, "Color misclassification rate", 0, 0.3))
      calc_lr_hair(as.numeric(input$indiv_MPc), as.numeric(input$indiv_obsC),
                   input$indiv_epc, DEFAULT_HAIR_FREQ)
    }
  })

  output$indiv_lr_value <- renderText({
    format_lr(indiv_lr())
  })

  output$indiv_lr_log <- renderText({
    format_log_lr(indiv_lr())
  })

  output$indiv_formula <- renderPrint({
    LR <- indiv_lr()
    if (input$indiv_type == "sex") {
      match_status <- if (input$indiv_MPs == input$indiv_obsSex) "concordant" else "discordant"
      cat("Sex Evidence (", match_status, ")\n", sep = "")
      cat("----------------------------------------\n")
      if (match_status == "concordant") {
        cat("P(observed | H1) = 1 - eps =", 1 - input$indiv_eps_sex, "\n")
      } else {
        cat("P(observed | H1) = eps =", input$indiv_eps_sex, "\n")
      }
      obs_freq <- if (input$indiv_obsSex == "F") input$indiv_propF else (1 - input$indiv_propF)
      cat("P(observed | H2) =", round(obs_freq, 3), "\n")

    } else if (input$indiv_type == "age") {
      aw <- age_window(input$indiv_MPa, input$indiv_MPr)
      in_range <- input$indiv_obsAge >= aw$min && input$indiv_obsAge <= aw$max
      cat("Age Evidence (", if (in_range) "in range" else "out of range", ")\n", sep = "")
      cat("----------------------------------------\n")
      cat("MP age:", input$indiv_MPa, "+/-", input$indiv_MPr, "years\n")
      cat("Range used (ages 1-80):", aw$min, "to", aw$max,
          "=", aw$max - aw$min + 1, "integer ages\n")
      cat("Observed:", input$indiv_obsAge, "years\n")
      if (in_range) {
        cat("P(in range | H1) = 1 - eps =", 1 - input$indiv_epa, "\n")
        cat("P(in range | H2) =", aw$max - aw$min + 1, "/ 80 =", round(aw$T1p, 4), "\n")
      } else {
        cat("P(out of range | H1) = eps =", input$indiv_epa, "\n")
        cat("P(out of range | H2) = 1 -", aw$max - aw$min + 1, "/ 80 =",
            round(1 - aw$T1p, 4), "\n")
      }

    } else {
      match_status <- if (input$indiv_MPc == input$indiv_obsC) "concordant" else "discordant"
      cat("Hair Color Evidence (", match_status, ")\n", sep = "")
      cat("----------------------------------------\n")
      cat("MP:", hair_labels[input$indiv_MPc], "\n")
      cat("Observed:", hair_labels[input$indiv_obsC], "\n")
      if (match_status == "concordant") {
        cat("P(observed | H1) = 1 - eps =", 1 - input$indiv_epc, "\n")
      } else {
        cat("P(observed | H1) = eps / 4 =", input$indiv_epc / 4, "\n")
      }
      cat("P(observed | H2) = population frequency =",
          DEFAULT_HAIR_FREQ[as.numeric(input$indiv_obsC)], "\n")
    }
    cat(if (LR == 0) "LR = 0 (exclusion)" else paste("LR =", round(LR, 4)), "\n")
  })

  cpt_pop <- reactive({
    validate(
      check_num(input$cpt_MPa, "Age", 1, 80, whole = TRUE),
      check_num(input$cpt_MPr, "Age range", 1, 20, whole = TRUE),
      check_num(input$cpt_propF, "Female proportion", 0.3, 0.7)
    )
    build_cpt_pop(input$cpt_propF, input$cpt_MPa, input$cpt_MPr, DEFAULT_HAIR_FREQ)
  })

  cpt_mp <- reactive({
    validate(
      check_num(input$cpt_eps, "Sex error rate", 0, 0.5),
      check_num(input$cpt_epa, "Age error rate", 0, 0.5),
      check_num(input$cpt_epc, "Hair error rate", 0, 0.5)
    )
    build_cpt_mp(input$cpt_MPs, as.numeric(input$cpt_MPc),
                 input$cpt_eps, input$cpt_epa, input$cpt_epc)
  })

  output$cpt_plots <- renderPlot({
    create_cpt_plots(cpt_pop(), cpt_mp(), narrow = is_narrow("cpt_plots"))
  }, res = 100, height = function() if (is_narrow("cpt_plots")) 1000 else 380)

  lr_dists <- reactive({
    build_lr_distributions(cpt_pop(), cpt_mp())
  })

  output$dist_plot <- renderPlot({
    dists <- lr_dists()
    create_dist_plot(dists$H1, dists$H2, base_size = plot_size("dist_plot", 12))
  }, res = 100)

  dist_stats <- function(d) {
    if (nrow(d) == 0) {
      cat("No data available\n")
    } else {
      p0 <- sum(d$w[d$x == -Inf])
      xf <- d$x[is.finite(d$x)]
      if (p0 > 0) {
        cat("Weighted mean log10(LR): undefined, P(LR = 0) =", sprintf("%.3f", p0), "\n")
        cat("Range: LR = 0 to log10(LR) =", sprintf("%.3f", max(xf)), "\n")
      } else {
        cat("Weighted mean log10(LR):", sprintf("%.3f", sum(d$x * d$w)), "\n")
        cat("Range: [", sprintf("%.3f", min(xf)), ",", sprintf("%.3f", max(xf)), "]\n")
      }
    }
  }

  output$dist_stats_h1 <- renderPrint({
    dist_stats(lr_dists()$H1)
  })

  output$dist_stats_h2 <- renderPrint({
    dist_stats(lr_dists()$H2)
  })

  output$dist_roc <- renderPlot({
    dists <- lr_dists()
    roc_obj <- build_roc(dists$H1, dists$H2)
    df_roc <- roc_obj$data
    auc_val <- roc_obj$AUC

    if (nrow(df_roc) == 0) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") +
               ggplot2::theme_void())
    }

    ggplot2::ggplot(df_roc, ggplot2::aes(x = FPR, y = TPR)) +
      ggplot2::geom_path(color = "#1f3a5f", linewidth = 1) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                           color = "gray50", linewidth = 0.8) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(title = sprintf("ROC Curve (AUC = %.3f)", auc_val),
                    subtitle = "Discrimination between H1 and H2",
                    x = "False Positive Rate (1 - Specificity)",
                    y = "True Positive Rate (Sensitivity)") +
      ggplot2::theme_minimal(base_size = plot_size("dist_roc", 12)) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }, res = 100)

  comb_lr_sex_val <- reactive({
    if (input$comb_use_sex) {
      validate(
        check_num(input$comb_eps, "Sex error", 0, 0.3),
        check_num(input$comb_propF, "Female proportion in population", 0.3, 0.7)
      )
      calc_lr_sex(input$comb_MPs, input$comb_obsSex, input$comb_eps, input$comb_propF)
    } else {
      1
    }
  })

  comb_lr_age_val <- reactive({
    if (input$comb_use_age) {
      validate(
        check_num(input$comb_MPa, "MP age", 1, 80, whole = TRUE),
        check_num(input$comb_MPr, "Age range", 1, 20, whole = TRUE),
        check_num(input$comb_obsAge, "Observed age", 1, 80, whole = TRUE),
        check_num(input$comb_epa, "Age error", 0, 0.3)
      )
      calc_lr_age(input$comb_MPa, input$comb_MPr, input$comb_obsAge, input$comb_epa)
    } else {
      1
    }
  })

  comb_lr_hair_val <- reactive({
    if (input$comb_use_hair) {
      validate(check_num(input$comb_epc, "Hair error", 0, 0.3))
      calc_lr_hair(as.numeric(input$comb_MPc), as.numeric(input$comb_obsC),
                   input$comb_epc, DEFAULT_HAIR_FREQ)
    } else {
      1
    }
  })

  output$comb_lr_sex <- renderText({
    if (input$comb_use_sex) format_lr(comb_lr_sex_val()) else "Not included"
  })

  output$comb_lr_age <- renderText({
    if (input$comb_use_age) format_lr(comb_lr_age_val()) else "Not included"
  })

  output$comb_lr_hair <- renderText({
    if (input$comb_use_hair) format_lr(comb_lr_hair_val()) else "Not included"
  })

  comb_total <- reactive({
    comb_lr_sex_val() * comb_lr_age_val() * comb_lr_hair_val()
  })

  output$comb_lr_total <- renderText({
    format_lr(comb_total())
  })

  output$comb_lr_log <- renderText({
    format_log_lr(comb_total())
  })

  output$comb_formula <- renderPrint({
    vals <- c()
    if (input$comb_use_sex) vals <- c(vals, comb_lr_sex_val())
    if (input$comb_use_age) vals <- c(vals, comb_lr_age_val())
    if (input$comb_use_hair) vals <- c(vals, comb_lr_hair_val())

    if (length(vals) == 0) {
      cat("No evidence selected\n")
    } else {
      cat(paste(sprintf("%.4g", vals), collapse = " x "), "\n")
      cat("=", sprintf("%.4g", prod(vals)), "\n")
    }
  })

  output$dec_metrics <- renderPrint({
    dists <- lr_dists()
    m <- calc_metrics(dists$H1, dists$H2, input$dec_threshold)

    cat("At threshold log10(LR) =", input$dec_threshold, "\n")
    cat("Match declared when log10(LR) >= threshold\n")
    cat("========================================\n\n")
    cat("True Positive Rate:  ", sprintf("%.3f", m$TPR), "\n")
    cat("False Positive Rate: ", sprintf("%.3f", m$FPR), "\n")
    cat("False Negative Rate: ", sprintf("%.3f", m$FNR), "\n")
    cat("True Negative Rate:  ", sprintf("%.3f", m$TNR), "\n")
    cat("\nMCC (balanced classes): ", sprintf("%.3f", m$MCC), "\n")
  })

  output$dec_roc_plot <- renderPlot({
    dists <- lr_dists()
    roc_obj <- build_roc(dists$H1, dists$H2)
    df_roc <- roc_obj$data
    auc_val <- roc_obj$AUC

    if (nrow(df_roc) == 0) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") +
               ggplot2::theme_void())
    }

    m <- calc_metrics(dists$H1, dists$H2, input$dec_threshold)
    current_point <- data.frame(FPR = m$FPR, TPR = m$TPR)

    ggplot2::ggplot(df_roc, ggplot2::aes(x = FPR, y = TPR)) +
      ggplot2::geom_path(color = "#1f3a5f", linewidth = 1) +
      ggplot2::geom_point(data = current_point, ggplot2::aes(x = FPR, y = TPR),
                         color = "#a33a36", size = 5, shape = 16) +
      ggplot2::geom_point(data = current_point, ggplot2::aes(x = FPR, y = TPR),
                         color = "white", size = 3, shape = 16) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                           color = "gray50", linewidth = 0.8) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(title = sprintf("ROC Curve (AUC = %.3f)", auc_val),
                    subtitle = sprintf("Red point: current threshold (log10(LR) >= %.1f)",
                                      input$dec_threshold),
                    x = "False Positive Rate", y = "True Positive Rate") +
      ggplot2::theme_minimal(base_size = plot_size("dec_roc_plot", 14)) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }, res = 100)

  sim_families_val <- reactiveVal(NULL)
  sim_parse_error  <- reactiveVal(NULL)

  observeEvent(input$sim_file, {
    result <- tryCatch(parse_sim_upload(input$sim_file$datapath), error = function(e) NULL)
    if (is.null(result)) {
      sim_families_val(NULL)
      sim_parse_error(paste0(
        "Could not parse the uploaded file ('", input$sim_file$name, "'). ",
        "Please check the accepted formats below.\n\n",
        "Common issues:\n",
        "- The file must be a text file (.txt, .csv, .tsv), not Excel (.xlsx)\n",
        "- For Familias format: Row 1 = family names, Row 2 = pedigree labels ",
        "(e.g., 'Reference pedigree', 'Father'), then numeric LR values\n",
        "- For FamLink format: must start with 'Raw data LR' and have ",
        "'Pedigree (Numerator)', 'Versus (Denominator)' and 'True:' header rows\n",
        "- For generic CSV/TSV: need at least 2 numeric columns\n",
        "- Check that the separator is tab, comma, or semicolon"
      ))
    } else {
      sim_families_val(result)
      sim_parse_error(NULL)
    }
  })

  output$sim_example_help <- renderText({
    genetic_example_help(input$sim_example)
  })

  output$sim_example_pedigree <- renderUI({
    req(input$sim_example %in% genetic_examples)
    if (identical(input$sim_example, "family_tree.txt")) return(NULL)
    img_file <- sub("\\.txt$", ".png", input$sim_example)
    div(style = "margin-top: 8px;",
      tags$img(src = file.path("pedigrees", img_file),
               alt = paste("Pedigree:", names(genetic_examples)[genetic_examples == input$sim_example]),
               style = "width: 100%; height: auto; border: 1px solid #ddd; background: white;"))
  })

  sim_example_path <- reactive({
    req(input$sim_example %in% genetic_examples)
    file.path("examples", "genetic", input$sim_example)
  })

  observeEvent(input$sim_load_example, {
    path <- sim_example_path()
    if (file.exists(path)) {
      sim_families_val(parse_sim_upload(path))
      sim_parse_error(NULL)
    } else {
      sim_parse_error("The example dataset could not be found on the server.")
    }
  })

  output$sim_download_example <- downloadHandler(
    filename = function() input$sim_example,
    content = function(file) file.copy(sim_example_path(), file)
  )

  output$sim_family_ui <- renderUI({
    fams <- sim_families_val()
    if (is.null(fams) || length(fams) <= 1) return(NULL)
    div(style = "margin-top: 10px;",
      hr(),
      h5("Select Comparison"),
      selectInput("sim_family", NULL, choices = names(fams), selected = names(fams)[1]),
      p(class = "text-muted small",
        sprintf("%d comparisons detected in this file.", length(fams)))
    )
  })

  output$sim_parse_msg <- renderUI({
    err <- sim_parse_error()
    if (!is.null(err)) {
      err_lines <- strsplit(err, "\n")[[1]]
      return(div(class = "warning-box", style = "margin-top: 10px; padding: 10px;",
        tags$b("File format error"),
        tags$pre(style = "white-space: pre-wrap; margin-top: 8px; font-size: 12px;
                         background: transparent; border: none; padding: 0;",
          paste(err_lines, collapse = "\n")
        )
      ))
    }
    NULL
  })

  sim_data <- reactive({
    fams <- sim_families_val()
    if (is.null(fams) || length(fams) == 0) return(NULL)
    if (length(fams) == 1) return(fams[[1]])
    sel <- input$sim_family
    if (is.null(sel) || !(sel %in% names(fams))) return(fams[[1]])
    fams[[sel]]
  })

  output$sim_has_data <- reactive({
    !is.null(sim_data())
  })
  outputOptions(output, "sim_has_data", suspendWhenHidden = FALSE)

  output$sim_data_info <- renderText({
    df <- sim_data()
    req(df)
    fams <- sim_families_val()
    n_fams <- if (!is.null(fams)) length(fams) else 1
    if (n_fams > 1) {
      sprintf("Data loaded: %s simulations (%d comparisons available)",
              format(nrow(df), big.mark = ","), n_fams)
    } else {
      sprintf("Data loaded: %s simulations", format(nrow(df), big.mark = ","))
    }
  })

  output$sim_n_sims <- renderText({
    df <- sim_data()
    req(df)
    format(nrow(df), big.mark = ",")
  })

  output$sim_ip <- renderText({
    df <- sim_data()
    req(df)
    sprintf("%.3f", mean(df$Related >= 10000))
  })

  output$sim_ep <- renderText({
    df <- sim_data()
    req(df)
    ep <- sum(df$Unrelated == 0 | df$Unrelated < 1e-10) / nrow(df)
    sprintf("%.3f", ep)
  })

  sim_summary <- function(x) {
    cat("Raw LR values:\n")
    print(summary(x))
    cat("\nlog10(LR) values (LR > 0 only):\n")
    valid <- x[x > 0]
    if (length(valid) > 0) print(summary(log10(valid)))
    if (any(x == 0)) {
      cat("\nLR = 0 in", sum(x == 0), "of", length(x), "simulations\n")
    }
  }

  output$sim_summary_rel <- renderPrint({
    df <- sim_data()
    req(df)
    sim_summary(df$Related)
  })

  output$sim_summary_unrel <- renderPrint({
    df <- sim_data()
    req(df)
    sim_summary(df$Unrelated)
  })

  zero_caption <- function(df) {
    if (any(df$Related == 0) || any(df$Unrelated == 0)) {
      sprintf("Not shown (LR = 0): %d related and %d unrelated simulations",
              sum(df$Related == 0), sum(df$Unrelated == 0))
    } else {
      NULL
    }
  }

  output$sim_violin <- renderPlot({
    df <- sim_data()
    req(df)

    rel_valid <- df$Related[df$Related > 0]
    unrel_valid <- df$Unrelated[df$Unrelated > 0]

    df_long <- data.frame(
      Hypothesis = c(
        rep("Related (H1)", length(rel_valid)),
        rep("Unrelated (H2)", length(unrel_valid))
      ),
      log10LR = c(log10(rel_valid), log10(unrel_valid))
    )

    ggplot2::ggplot(df_long, ggplot2::aes(x = Hypothesis, y = log10LR, fill = Hypothesis)) +
      ggplot2::geom_violin(alpha = 0.7, color = "white", linewidth = 0.3) +
      ggplot2::geom_boxplot(width = 0.15, fill = "white", alpha = 0.8,
                            outlier.size = 0.5, outlier.alpha = 0.3) +
      ggplot2::scale_fill_manual(values = c("Related (H1)" = "#1f3a5f",
                                            "Unrelated (H2)" = "#999999")) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      ggplot2::labs(title = "LR Distribution Comparison",
                    subtitle = "Dashed line: LR = 1 (neutral evidence)",
                    caption = zero_caption(df),
                    x = "", y = expression(log[10](LR))) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = plot_size("sim_violin", 13)) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank()
      )
  }, res = 100)

  output$sim_density_plot <- renderPlot({
    df <- sim_data()
    req(df)

    rel_valid <- df$Related[df$Related > 0]
    unrel_valid <- df$Unrelated[df$Unrelated > 0]

    if (length(rel_valid) < 2 || length(unrel_valid) < 2) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Insufficient positive LR values") +
               ggplot2::theme_void())
    }

    dens_rel <- stats::density(log10(rel_valid))
    dens_unrel <- stats::density(log10(unrel_valid))

    df_dens <- rbind(
      data.frame(x = dens_rel$x, y = dens_rel$y, Hypothesis = "Related (H1)"),
      data.frame(x = dens_unrel$x, y = dens_unrel$y, Hypothesis = "Unrelated (H2)")
    )

    ggplot2::ggplot(df_dens, ggplot2::aes(x = x, y = y, fill = Hypothesis)) +
      ggplot2::geom_area(alpha = 0.55, position = "identity", color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = c("Related (H1)" = "#1f3a5f",
                                            "Unrelated (H2)" = "#999999")) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
      ggplot2::labs(title = "LR Distribution Density",
                    subtitle = "Vertical dashed line: LR = 1 (neutral). Less overlap = better discrimination.",
                    caption = zero_caption(df),
                    x = expression(log[10](LR)), y = "Density") +
      ggplot2::theme_minimal(base_size = plot_size("sim_density_plot", 13)) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  }, res = 100)

  sim_log_stats <- function(x) {
    valid <- x[x > 0]
    if (length(valid) == 0) {
      cat("No positive values\n")
    } else {
      lv <- log10(valid)
      cat("Mean:   ", sprintf("%.3f", mean(lv)), "\n")
      cat("Median: ", sprintf("%.3f", stats::median(lv)), "\n")
      cat("SD:     ", sprintf("%.3f", stats::sd(lv)), "\n")
      cat("Range:  [", sprintf("%.3f", min(lv)), ",", sprintf("%.3f", max(lv)), "]\n")
    }
    if (any(x == 0)) {
      cat("Excluded (LR = 0):", sum(x == 0), "of", length(x), "\n")
    }
  }

  output$sim_dist_stats_rel <- renderPrint({
    df <- sim_data()
    req(df)
    sim_log_stats(df$Related)
  })

  output$sim_dist_stats_unrel <- renderPrint({
    df <- sim_data()
    req(df)
    sim_log_stats(df$Unrelated)
  })

  output$sim_decision_plot <- renderPlot({
    df <- sim_data()
    req(df)
    validate(check_num(input$sim_threshold, "LR threshold", 0, 1e6))

    curve_data <- compute_decision_curve(df, n_points = 500)

    opt <- sim_opt()
    opt_point <- data.frame(LR = opt$DT, FPR = opt$FPR, FNR = opt$FNR)

    sel_point <- data.frame(FPR = mean(df$Unrelated >= input$sim_threshold),
                            FNR = mean(df$Related < input$sim_threshold))

    label_indices <- unique(c(1,
                              round(nrow(curve_data) * 0.25),
                              round(nrow(curve_data) * 0.5),
                              round(nrow(curve_data) * 0.75),
                              nrow(curve_data)))
    label_data <- curve_data[label_indices, ]

    ggplot2::ggplot(curve_data, ggplot2::aes(x = FNR, y = FPR, color = log10(LR))) +
      ggplot2::geom_point(size = 2, alpha = 0.8) +
      ggplot2::scale_color_gradient2(
        low = "#cccccc", mid = "#7f93ad", high = "#1f3a5f",
        midpoint = stats::median(log10(curve_data$LR)),
        name = expression(log[10](LR~threshold))
      ) +
      ggplot2::geom_point(data = opt_point, ggplot2::aes(x = FNR, y = FPR),
                         color = "#222222", size = 5, shape = 18, inherit.aes = FALSE) +
      ggplot2::annotate("text",
        x = opt_point$FNR, y = opt_point$FPR,
        label = sprintf("Optimal DT\nLR = %.2f", opt_point$LR),
        hjust = 0.5, vjust = -0.6, size = 3.5, fontface = "bold", color = "#222222") +
      ggplot2::geom_point(data = sel_point, ggplot2::aes(x = FNR, y = FPR),
                         color = "#a33a36", size = 5, shape = 16, inherit.aes = FALSE) +
      ggplot2::geom_point(data = sel_point, ggplot2::aes(x = FNR, y = FPR),
                         color = "white", size = 2.5, shape = 16, inherit.aes = FALSE) +
      ggplot2::annotate("text",
        x = sel_point$FNR, y = sel_point$FPR,
        label = sprintf("Selected\nLR = %g", input$sim_threshold),
        hjust = 0.5, vjust = 1.5, size = 3.2, color = "#a33a36") +
      ggplot2::geom_text(data = label_data,
        ggplot2::aes(x = FNR, y = FPR, label = round(LR)),
        size = 2.8, vjust = -1, color = "gray40", inherit.aes = FALSE,
        check_overlap = TRUE) +
      ggplot2::scale_x_continuous(n.breaks = if (is_narrow("sim_decision_plot")) 4 else 6) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.15))) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = "LR Decision Threshold Plot",
        subtitle = if (is_narrow("sim_decision_plot")) {
          "Black diamond: optimal DT\nRed dot: selected threshold\n(match if LR >= threshold)"
        } else {
          "Black diamond: optimal DT | Red dot: selected threshold\n(match if LR >= threshold)"
        },
        x = "False Negative Rate (FNR)",
        y = "False Positive Rate (FPR)"
      ) +
      ggplot2::theme_minimal(base_size = plot_size("sim_decision_plot", 13)) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = if (is_narrow("sim_decision_plot")) "bottom" else "right",
        panel.grid.minor = ggplot2::element_blank()
      )
  }, res = 100)

  sim_rates <- reactive({
    df <- sim_data()
    req(df)
    threshold <- input$sim_threshold
    validate(check_num(threshold, "LR threshold", 0, 1e6))

    n <- nrow(df)
    TP <- sum(df$Related >= threshold)
    FN <- sum(df$Related < threshold)
    FP <- sum(df$Unrelated >= threshold)
    TN <- sum(df$Unrelated < threshold)

    TPR <- TP / n
    TNR <- TN / n
    FPR_val <- FP / n
    FNR_val <- FN / n
    Accuracy <- (TP + TN) / (TP + TN + FP + FN)

    denom <- sqrt(as.numeric(TP + FP) * as.numeric(TP + FN) *
                  as.numeric(TN + FP) * as.numeric(TN + FN))
    MCC <- if (denom == 0) 0 else (TP * TN - FP * FN) / denom

    list(TPR = TPR, TNR = TNR, FPR = FPR_val, FNR = FNR_val,
         MCC = MCC, Accuracy = Accuracy, FP_count = FP, FN_count = FN,
         TP_count = TP, TN_count = TN)
  })

  sim_n_message <- reactive({
    check_num(input$sim_db_size, "Number of unrelated candidates (N)", 1, Inf, whole = TRUE)
  })

  sim_n_candidates <- reactive({
    req(is.null(sim_n_message()))
    input$sim_db_size
  })

  sim_ppv <- reactive({
    r <- sim_rates()
    denom <- r$TPR + r$FPR * sim_n_candidates()
    if (denom == 0) NA else r$TPR / denom
  })

  output$sim_rate_tpr <- renderText({
    sprintf("%.4f", sim_rates()$TPR)
  })

  output$sim_rate_fpr <- renderText({
    sprintf("%.4f", sim_rates()$FPR)
  })

  output$sim_rate_fnr <- renderText({
    sprintf("%.4f", sim_rates()$FNR)
  })

  output$sim_rate_tnr <- renderText({
    sprintf("%.4f", sim_rates()$TNR)
  })

  output$sim_rate_ppv <- renderText({
    if (is.na(sim_ppv())) "undefined" else sprintf("%.4f", sim_ppv())
  })

  output$sim_rate_acc <- renderText({
    sprintf("%.4f", sim_rates()$Accuracy)
  })

  output$sim_expected_fp <- renderText({
    validate(sim_n_message())
    sprintf("%.1f", sim_rates()$FPR * sim_n_candidates())
  })

  output$sim_rates_table <- renderTable({
    r <- sim_rates()
    data.frame(
      Metric = c("True Positive Rate (TPR)",
                 "False Positive Rate (FPR)",
                 "True Negative Rate (TNR)",
                 "False Negative Rate (FNR)",
                 sprintf("Expected false positives (FPR x N, N = %s)",
                         format(sim_n_candidates(), big.mark = ",")),
                 "PPV (1 true match among N candidates)",
                 "Accuracy (balanced classes)",
                 "Matthews Correlation Coefficient (MCC, balanced classes)"),
      Value = c(sprintf("%.6f", r$TPR),
                sprintf("%.6f", r$FPR),
                sprintf("%.6f", r$TNR),
                sprintf("%.6f", r$FNR),
                sprintf("%.1f", r$FPR * sim_n_candidates()),
                if (is.na(sim_ppv())) "undefined" else sprintf("%.6f", sim_ppv()),
                sprintf("%.6f", r$Accuracy),
                sprintf("%.6f", r$MCC))
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  sim_opt <- reactive({
    df <- sim_data()
    req(df)
    compute_optimal_dt(df, weight = input$sim_weight)
  })

  sim_dt_shown <- reactiveVal(FALSE)

  observeEvent(sim_data(), {
    sim_dt_shown(FALSE)
  })

  observeEvent(input$sim_calc_dt, {
    req(sim_data())
    sim_dt_shown(TRUE)
  })

  output$sim_dt_result <- renderPrint({
    if (!sim_dt_shown()) {
      cat("Click 'Compute Optimal DT' to calculate.\n")
    } else {
      result <- sim_opt()
      cat("Optimal Decision Threshold\n")
      cat("========================================\n\n")
      cat("  LR threshold: ", sprintf("%.2f", result$DT), "\n")
      cat("  FPR at DT:    ", sprintf("%.6f", result$FPR), "\n")
      cat("  FNR at DT:    ", sprintf("%.6f", result$FNR), "\n")
      cat("  Weight used:  ", input$sim_weight, "\n")
      cat("\n  Formula: D = sqrt(FNR^2 + (w*FPR)^2)\n")
      cat("  Match declared when LR >= threshold\n")
    }
  })

  dbs_data <- reactiveVal(NULL)
  dbs_error <- reactiveVal(NULL)

  dbs_set <- function(result) {
    if (is.character(result)) {
      dbs_data(NULL)
      dbs_error(result)
    } else {
      dbs_data(result)
      dbs_error(NULL)
    }
  }

  observeEvent(input$dbs_file, {
    dbs_set(tryCatch(read_candidates(input$dbs_file$datapath),
                     error = function(e) "The file could not be read."))
  })

  observeEvent(input$dbs_load_example, {
    dbs_set(read_candidates(candidates_example_path))
  })

  output$dbs_download_example <- downloadHandler(
    filename = function() "candidates_example.csv",
    content = function(file) file.copy(candidates_example_path, file)
  )

  output$dbs_msg <- renderUI({
    err <- dbs_error()
    if (!is.null(err)) {
      div(class = "warning-box", style = "padding: 10px;", tags$b("File error: "), err)
    } else if (!is.null(dbs_data())) {
      div(class = "success-box", style = "padding: 10px;",
          sprintf("Loaded %s candidates.", format(nrow(dbs_data()), big.mark = ",")))
    }
  })

  dbs_lr <- reactive({
    cands <- dbs_data()
    validate(need(!is.null(cands), "Upload a candidate database or load the example."))
    validate(
      check_num(input$dbs_MPa, "Age", 1, 80, whole = TRUE),
      check_num(input$dbs_MPr, "Age range", 1, 20, whole = TRUE),
      check_num(input$dbs_eps, "Sex error rate", 0, 0.5),
      check_num(input$dbs_epa, "Age error rate", 0, 0.5),
      check_num(input$dbs_epc, "Hair error rate", 0, 0.5)
    )
    candidate_lrs(cands, input$dbs_MPs, input$dbs_MPa, input$dbs_MPr, as.numeric(input$dbs_MPc),
                  input$dbs_eps, input$dbs_epa, input$dbs_epc, isTRUE(input$dbs_empirical))
  })

  output$dbs_summary <- renderPrint({
    out <- dbs_lr()
    res <- out$results
    n <- nrow(res)
    above <- sum(res$log10LR >= input$dbs_threshold)
    cat("Candidates:", n, "\n")
    cat(sprintf("LR >= %s (log10(LR) >= %.1f): %d (%.4f of candidates)\n",
                format(signif(10^input$dbs_threshold, 4)), input$dbs_threshold, above, above / n))
    cat("LR = 0 (exclusions):", sum(res$LR == 0), "\n\n")
    cat(if (isTRUE(input$dbs_empirical)) "Frequencies estimated from the database:\n" else
          "Default population frequencies:\n")
    cat(sprintf("  P(F) = %.4f\n", out$propF))
    cat(sprintf("  P(age %g-%g | H2) = %.4f\n", out$age_min, out$age_max, out$T1p))
    cat("  Hair:", paste(sprintf("%s %.4f", names(out$hair_freq), out$hair_freq), collapse = ", "), "\n")
    if (length(out$notes) > 0) cat("\nNote:", paste(out$notes, collapse = "\nNote: "), "\n")
  })

  output$dbs_hist <- renderPlot({
    res <- dbs_lr()$results
    res$w <- 1 / nrow(res)
    fin <- res[is.finite(res$log10LR), ]
    p0 <- sum(res$LR == 0) / nrow(res)
    binwidth <- 0.25
    validate(need(nrow(fin) > 0, "All candidates have LR = 0."))
    p <- ggplot2::ggplot(fin, ggplot2::aes(x = log10LR, weight = w)) +
      ggplot2::geom_histogram(binwidth = binwidth, fill = "#999999", alpha = 0.8,
                              color = "white", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = input$dbs_threshold, linetype = "dashed",
                          color = "gray30", linewidth = 0.8) +
      ggplot2::labs(title = "Candidate log10(LR)",
                    subtitle = "Dashed line: threshold",
                    x = expression(log[10](LR)), y = "Probability") +
      ggplot2::theme_minimal(base_size = plot_size("dbs_hist", 12)) +
      ggplot2::theme(plot.title.position = "plot",
                     plot.title = ggplot2::element_text(face = "bold"),
                     panel.grid.minor = ggplot2::element_blank())
    if (p0 > 0) {
      zx <- min(fin$log10LR) - 3 * binwidth
      zero <- data.frame(log10LR = zx, w = p0)
      p <- p +
        ggplot2::geom_col(data = zero, ggplot2::aes(x = log10LR, y = w), inherit.aes = FALSE,
                          width = binwidth, fill = "#999999", alpha = 0.8) +
        ggplot2::geom_text(data = zero, ggplot2::aes(x = log10LR, y = w,
                           label = sprintf("LR = 0\nP = %.3f", w)),
                           inherit.aes = FALSE, vjust = -0.3, size = 3.2) +
        ggplot2::geom_vline(xintercept = zx + 1.5 * binwidth, linetype = "dotted", color = "gray50") +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
        ggplot2::labs(subtitle = "Dashed line: threshold. Left bar: LR = 0 (exclusion)")
    }
    p
  }, res = 100)

  output$dbs_table <- renderTable({
    res <- dbs_lr()$results
    validate(check_num(input$dbs_top, "Top candidates shown", 1, 500, whole = TRUE))
    top <- utils::head(res, input$dbs_top)
    data.frame(
      Rank = seq_len(nrow(top)),
      ID = top$ID, Sex = top$Sex, Age = as.integer(top$Age), HairColor = top$HairColor,
      LR = sub("LR = ", "", vapply(top$LR, format_lr, character(1)), fixed = TRUE),
      log10LR = ifelse(top$LR == 0, "LR = 0 (exclusion)", sprintf("%.3f", top$log10LR)),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  output$dbs_download_results <- downloadHandler(
    filename = function() "database_search_results.csv",
    content = function(file) {
      res <- dbs_lr()$results
      utils::write.csv(data.frame(Rank = seq_len(nrow(res)), res), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
