library(shiny)
library(shinythemes)
library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)
library(shiny.i18n)

i18n <- Translator$new(translation_json_path = file.path("translations", "translation.json"))
i18n$set_translation_language("en")
i18n$use_js()

tr <- function(key) i18n$t(key)

trf <- function(key, ...) sprintf(i18n$t(key), ...)

tr_span <- function(key) span(class = "i18n", `data-key` = key, isolate(i18n$t(key)))

big_mark <- function() {
  s <- shiny::getDefaultReactiveDomain()
  lang <- if (is.null(s) || is.null(s$userData$shiny.i18n)) "en" else s$userData$shiny.i18n$lang()
  if (identical(lang, "es")) "\u00a0" else ","
}

tr_names <- function(choices) setNames(choices, vapply(names(choices), tr, character(1)))

fmt_msg <- function(m) do.call(sprintf, c(list(tr(m[1])), as.list(m[-1])))

COL_H1 <- "#0072B2"
COL_H2 <- "#E69F00"
COL_THRESHOLD <- "#000000"
COL_OPTIMAL <- "#009E73"
COL_MP <- "#CC79A7"
COL_NEUTRAL <- "#999999"
COL_CASE <- "#6A3D9A"
COL_REFERENCE_LINE <- "grey60"

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

hair_error_help <- "Total error, split equally among the other four colours."

hair_freq_help <- "Population frequencies: Black 0.08, Brown 0.45, Blonde 0.30, Red 0.07, Gray 0.10."

defaults_help <- "Defaults: Marsico & Caridi 2023."

interpretation_templates <- list(
  single = c("The observation is %s times more probable under H1 than under H2.",
             "The observation is %s times more probable under H2 than under H1."),
  combined = c("The observations are %s times more probable under H1 than under H2.",
               "The observations are %s times more probable under H2 than under H1.")
)

mcc_help <- "Matthews correlation between decision and truth, from -1 to 1; balanced classes."

check_num <- function(x, label, lo, hi, whole = FALSE) {
  if (is.null(x) || length(x) != 1 || !is.numeric(x) || !is.finite(x)) {
    return(trf("%s: enter a number.", tr(label)))
  }
  if (x < lo || x > hi) {
    if (is.infinite(hi)) return(trf("%s: must be at least %s.", tr(label), lo))
    return(trf("%s: must be between %s and %s.", tr(label), lo, hi))
  }
  if (whole && x != round(x)) {
    return(trf("%s: must be a whole number.", tr(label)))
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
  if (LR == 0) tr("LR = 0 (exclusion)") else sprintf("log10(LR) = %.3f", log10(LR))
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

create_cpt_plots <- function(CPT_POP, CPT_MP, narrow = FALSE, case_cell = NULL) {
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

  hair_axis <- list(ggplot2::scale_x_discrete(labels = function(x) vapply(x, tr, character(1))),
                    ggplot2::guides(x = ggplot2::guide_axis(n.dodge = 2)))

  p1 <- ggplot2::ggplot(dfPOP, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient(low = "#ffffff", high = "#969696",
                                 limits = c(0, max(dfPOP$value) * 1.1)) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)),
                       size = 2.8, color = "black") +
    hair_axis +
    ggplot2::labs(title = "P(D | H2)", subtitle = tr("H2: unrelated person (population)"),
                  x = tr("Hair Color"), y = "", fill = tr("Probability")) +
    base_theme

  p2 <- ggplot2::ggplot(dfMP, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient(low = "#ffffff", high = "#969696",
                                 limits = c(0, max(dfMP$value) * 1.1)) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)),
                       size = 2.8, color = "black") +
    hair_axis +
    ggplot2::labs(title = "P(D | H1)", subtitle = tr("H1: the person is the MP (related)"),
                  x = tr("Hair Color"), y = "", fill = tr("Probability")) +
    base_theme

  p3 <- ggplot2::ggplot(dfLR, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_tile(data = dfLR[is.na(dfLR$value), ], fill = "white",
                       color = "grey50", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(low = COL_H2, mid = "#ffffff", high = COL_H1,
                                  midpoint = 0, na.value = "white",
                                  labels = function(x) sprintf("%.1f", x)) +
    ggplot2::geom_text(ggplot2::aes(label = label),
                       size = 2.8, color = "black") +
    hair_axis +
    ggplot2::labs(title = expression(log[10](LR)),
                  subtitle = tr("Evidence weight (white cell: LR = 0)"),
                  x = tr("Hair Color"), y = "", fill = expression(log[10](LR))) +
    base_theme

  if (!is.null(case_cell)) {
    outline <- ggplot2::geom_tile(data = data.frame(Var1 = case_cell[1], Var2 = case_cell[2]),
                                  ggplot2::aes(x = Var2, y = Var1), inherit.aes = FALSE,
                                  fill = NA, color = COL_CASE, linewidth = 1.4)
    p1 <- p1 + outline
    p2 <- p2 + outline
    p3 <- p3 + outline
  }

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

weighted_median <- function(x, w) {
  o <- order(x)
  x[o][which(cumsum(w[o]) / sum(w) >= 0.5)[1]]
}

format_lr_value <- function(x) {
  if (x >= 1000) format(round(x), big.mark = big_mark(), scientific = FALSE) else format(signif(x, 3), scientific = FALSE)
}

format_log_value <- function(x) {
  if (x == -Inf) "LR = 0" else sprintf("%.2f", x)
}

log_axis_label <- function() {
  bquote(log[10](LR) ~ .(tr("(1 = LR 10, 2 = LR 100)")))
}

one_line <- function(x) gsub("\n", " ", x)

roc_title <- function(auc) {
  trf("Discrimination: AUC = %.3f\n(1 = perfect, 0.5 = chance)", auc)
}

dist_title <- function(distH1, distH2) {
  if (nrow(distH1) == 0 || nrow(distH2) == 0) return(tr("Insufficient data"))
  trf("Typical log10(LR): %s if the person is the MP,\n%s if unrelated",
      format_log_value(weighted_median(distH1$x, distH1$w)),
      format_log_value(weighted_median(distH2$x, distH2$w)))
}

create_dist_plot <- function(distH1, distH2, binwidth = 0.25, base_size = 12,
                             threshold = NULL, case_log = NULL) {
  if (nrow(distH1) == 0 || nrow(distH2) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::annotate("text", x = 0.5, y = 0.5, label = tr("Insufficient data")) +
             ggplot2::theme_void())
  }

  lab_h1 <- tr("H1: the person is the MP (related)")
  lab_h2 <- tr("H2: unrelated person")
  dH1 <- data.frame(logLR = distH1$x, w = distH1$w, Hypothesis = lab_h1)
  dH2 <- data.frame(logLR = distH2$x, w = distH2$w, Hypothesis = lab_h2)
  dAll <- rbind(dH1, dH2)
  dFinite <- dAll[is.finite(dAll$logLR), ]
  dZero <- dAll[dAll$logLR == -Inf, ]

  p <- ggplot2::ggplot(dFinite, ggplot2::aes(x = logLR, weight = w, fill = Hypothesis)) +
    ggplot2::geom_histogram(binwidth = binwidth, position = "identity",
                            alpha = 0.65, color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = setNames(c(COL_H1, COL_H2), c(lab_h1, lab_h2)),
                               breaks = c(lab_h1, lab_h2)) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::geom_vline(xintercept = 0, color = COL_REFERENCE_LINE, linewidth = 0.8) +
    ggplot2::labs(title = dist_title(distH1, distH2),
                  subtitle = tr("Grey line: LR = 1 (uninformative)"),
                  x = log_axis_label(), y = tr("Probability")) +
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
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.45))) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(subtitle = tr("Grey line: LR = 1 (uninformative).\nLeft bar: LR = 0 (exclusion)"))
  }
  marks <- character(0)
  if (!is.null(threshold)) {
    p <- p + ggplot2::geom_vline(xintercept = threshold, color = COL_THRESHOLD, linewidth = 0.8,
                                 linetype = "dashed")
    marks <- c(marks, tr("Black dashed line: threshold."))
  }
  if (!is.null(case_log)) {
    p <- p + ggplot2::geom_vline(xintercept = case_log, color = COL_CASE, linewidth = 1)
    marks <- c(marks, tr("Purple line: Case LR."))
  }
  if (length(marks) > 0) {
    p <- p + ggplot2::labs(caption = paste(marks, collapse = " "))
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
    tr("Familias example, 10,000 simulations.")
  } else {
    tr("forrel 1.9.1: 10,000 simulations per hypothesis, 24 autosomal STRs, Argentina allele frequencies, no mutation. Ref = genotyped relative.")
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
  if (is.character(data)) return(c("Could not read the file: %s", data))
  names(data) <- trimws(names(data))
  required <- c("ID", "Sex", "Age", "HairColor")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    return(c("Missing required column(s): %s. Required columns: ID, Sex, Age, HairColor.",
             paste(missing_cols, collapse = ", ")))
  }
  if (nrow(data) == 0) return("The file has no candidate rows.")
  data <- data[, required]
  data$Sex <- toupper(data$Sex)
  bad <- which(!nzchar(data$ID))
  if (length(bad) > 0) return(c("Column ID is empty in line(s) %s", row_list(bad)))
  bad <- which(duplicated(data$ID))
  if (length(bad) > 0) return(c("Column ID has duplicated values in line(s) %s", row_list(bad)))
  bad <- which(!data$Sex %in% c("F", "M"))
  if (length(bad) > 0) return(c("Column Sex must be F or M; invalid value in line(s) %s", row_list(bad)))
  age <- suppressWarnings(as.numeric(data$Age))
  bad <- which(is.na(age) | age < 1 | age > 80 | age != round(age))
  if (length(bad) > 0) {
    return(c("Column Age must be a whole number from 1 to 80; invalid value in line(s) %s", row_list(bad)))
  }
  hair <- match(tolower(data$HairColor), tolower(names(DEFAULT_HAIR_FREQ)))
  bad <- which(is.na(hair))
  if (length(bad) > 0) {
    return(c("Column HairColor must be one of Black, Brown, Blonde, Red, Gray; invalid value in line(s) %s",
             row_list(bad)))
  }
  data$Age <- age
  data$HairColor <- names(DEFAULT_HAIR_FREQ)[hair]
  data
}

candidate_lrs <- function(cands, MPs, MPa, MPr, MPc, eps, epa, epc, propF_default, empirical) {
  aw <- age_window(MPa, MPr)
  propF <- propF_default
  hair_freq <- DEFAULT_HAIR_FREQ
  T1p <- aw$T1p
  notes <- character(0)
  if (empirical) {
    propF <- mean(cands$Sex == "F")
    hair_freq <- vapply(names(DEFAULT_HAIR_FREQ), function(h) mean(cands$HairColor == h), numeric(1))
    T1p <- mean(cands$Age >= aw$min & cands$Age <= aw$max)
    if (propF %in% c(0, 1)) {
      propF <- propF_default
      notes <- c(notes, tr("Only one sex is present in the database: the female proportion of the MP panel is used."))
    }
    for (h in names(hair_freq)[hair_freq == 0]) {
      hair_freq[h] <- DEFAULT_HAIR_FREQ[h]
      notes <- c(notes, trf("No candidate has hair color %s: the default population frequency (%.2f) is used for it.",
                            tr(h), DEFAULT_HAIR_FREQ[h]))
    }
    if (T1p %in% c(0, 1)) {
      T1p <- aw$T1p
      notes <- c(notes, trf("The share of candidates inside the age range is %s: the uniform 1-80 age model is used.",
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
  usei18n(i18n),
  lang = "en",
  tags$head(
    tags$script(src = "lang.js"),
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

      a { color: #2b2b2b; }

      .main-header {
        background: #2b2b2b;
        color: #ffffff;
        padding: 10px 20px;
        margin: 0 -15px 15px -15px;
        display: flex;
        flex-wrap: wrap;
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

      .lang-switch {
        margin-left: auto;
        font-size: 13px;
        white-space: nowrap;
      }

      .lang-switch a, .lang-switch a:hover, .lang-switch a:focus {
        color: #ffffff;
        opacity: 0.75;
        text-decoration: none;
      }

      .lang-switch a.active {
        opacity: 1;
        font-weight: 600;
        text-decoration: underline;
      }

      .lang-switch .lang-sep {
        margin: 0 6px;
        opacity: 0.6;
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

      .h1-box {
        border-left: 4px solid #0072B2;
      }

      .h2-box {
        border-left: 4px solid #E69F00;
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

      pre.shiny-text-output {
        white-space: pre-wrap;
        word-break: normal;
        overflow-wrap: anywhere;
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

      .ped-legend {
        font-size: 12px;
        color: #222;
        margin-top: 6px;
        line-height: 1.8;
      }

      .ped-legend .ped-item {
        display: inline-block;
        margin-right: 12px;
        white-space: nowrap;
      }

      .ped-legend .swatch {
        display: inline-block;
        width: 12px;
        height: 12px;
        border: 1px solid #000000;
        vertical-align: -1px;
        margin-right: 5px;
        background: #ffffff;
      }

      .ped-legend .swatch-ref {
        background: repeating-linear-gradient(45deg, #000000 0, #000000 1.5px, #ffffff 1.5px, #ffffff 4px);
      }

      .ped-legend .swatch-mp {
        background: #CC79A7;
        border-color: #CC79A7;
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
        color: #2b2b2b;
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
        color: #2b2b2b;
        background: transparent;
        font-weight: 600;
        border-bottom: 2px solid #2b2b2b;
      }

      .btn {
        border-radius: 2px;
        box-shadow: none;
        transition: none;
      }

      .btn-primary, .btn-primary:hover, .btn-primary:focus, .btn-primary:active {
        background: #2b2b2b;
        border-color: #2b2b2b;
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
        background: #2b2b2b;
        border-color: #2b2b2b;
      }

      .table {
        color: #222;
      }

      .progress {
        border-radius: 2px;
        box-shadow: none;
      }

      .text-muted.small {
        color: #555;
        font-size: 13px;
      }

      .profile-banner {
        background: #f7f7f7;
        border: 1px solid #ddd;
        border-radius: 2px;
        padding: 8px 12px;
        margin: 10px 0;
      }

      .plain-sentence {
        font-size: 15px;
        line-height: 1.5;
        margin: 8px 0 0 0;
      }

      .entry-card h5 {
        line-height: 1.4;
      }

      .entry-card a {
        display: block;
        margin: 6px 0;
        text-decoration: underline;
      }

      details > summary {
        cursor: pointer;
        font-weight: 600;
      }

      details[open] > summary {
        margin-bottom: 8px;
      }
    "))
  ),

  div(class = "main-header",
    tags$img(src = "MispiIcon.png", alt = "mispitools logo"),
    div(
      h1("mispitools"),
      p(tr("Likelihood Ratio Analysis for Kinship in Forensic Science"))
    ),
    div(class = "lang-switch",
      tags$a(href = "#", class = "lang-link active", `data-lang` = "en", lang = "en", "English"),
      span(class = "lang-sep", "|"),
      tags$a(href = "#", class = "lang-link", `data-lang` = "es", lang = "es", "Español")
    )
  ),

  tabsetPanel(
    id = "main_tabs",
    type = "tabs",

    tabPanel(
      tr("Overview"), value = "Overview",
      fluidRow(
        column(10, offset = 1,
          br(),
          h2(tr("Evidence evaluation in kinship cases")),
          p(class = "interpretation-text",
            tr("Likelihood ratios for missing-person identification from non-genetic evidence (sex, age, hair colour) and from genetic simulation results. Non-genetic analyses run in the browser; genetic simulations are run elsewhere and uploaded.")
          ),
          div(class = "info-box",
            tags$ul(
              tags$li(strong("H1:"), " ", tr("the person is the missing person (MP)")),
              tags$li(strong("H2:"), " ", tr("the person is unrelated to the MP"))
            ),
            div(class = "formula-box",
              tr("LR = P(Evidence | H1) / P(Evidence | H2)")
            ),
            p(style = "margin: 0;", tr("The LR is not the probability that the person is the MP."))
          ),
          fluidRow(
            column(4,
              div(class = "well entry-card",
                h5(tr("LR of one observed person")),
                actionLink("go_case", tr("Non-genetic: case LR"))
              )
            ),
            column(4,
              div(class = "well entry-card",
                h5(tr("Expected false positives")),
                actionLink("go_gen_rates", tr("Genetic: rates at threshold")),
                actionLink("go_dbs", tr("Non-genetic: database search"))
              )
            ),
            column(4,
              div(class = "well entry-card",
                h5(tr("Choosing a threshold")),
                actionLink("go_choose", tr("Non-genetic: profile performance")),
                actionLink("go_gen_optimal", tr("Genetic: optimal threshold"))
              )
            )
          ),
          hr(),
          div(style = "color: #444; padding: 10px 0;",
            p(strong(tr("References"))),
            p("Egeland T., Marsico F.L. (2026). Using all available evidence to solve kinship cases.
              International Journal of Legal Medicine 140(4):1981-1996. ",
              tags$a(href = "https://doi.org/10.1007/s00414-026-03773-6", target = "_blank",
                     "https://doi.org/10.1007/s00414-026-03773-6")),
            p("Marsico F.L., Caridi I. (2023). Incorporating non-genetic evidence in large scale
              missing person searches: A general approach beyond filtering.
              Forensic Science International: Genetics 66:102891"),
            p("Marsico F.L., Vigeland M.D., Egeland T., Piñero M.H. (2021). Making decisions in
              missing person identification cases with low statistical power.
              Forensic Science International: Genetics 54:102519")
          )
        )
      )
    ),

    tabPanel(
      tr("Non-Genetic Evidence"), value = "Non-Genetic Evidence",
      br(),
      fluidRow(
        column(3,
          wellPanel(
            h5(tr("Missing person (MP)")),
            selectInput("ng_MPs", tr("Sex:"), sex_choices, selectize = FALSE),
            numericInput("ng_MPa", tr("Age (years):"), 40, 1, 80),
            numericInput("ng_MPr", tr("Age uncertainty (± years):"), 6, 1, 20),
            selectInput("ng_MPc", tr("Hair color:"), hair_choices, selected = 2, selectize = FALSE),
            hr(),
            h5(tr("Observation error rates (0–1)")),
            numericInput("ng_eps", tr("Sex:"), 0.05, 0, 0.5, 0.01),
            numericInput("ng_epa", tr("Age:"), 0.05, 0, 0.5, 0.01),
            numericInput("ng_epc", tr("Hair (total):"), 0.05, 0, 0.5, 0.01),
            p(class = "text-muted small", tr(hair_error_help)),
            hr(),
            h5(tr("Population")),
            sliderInput("ng_propF", tr("Female proportion (0–1):"), 0.3, 0.7, 0.51, 0.01),
            p(class = "text-muted small", tr(hair_freq_help)),
            p(class = "text-muted small", tr(defaults_help))
          )
        ),
        column(9,
          tabsetPanel(
            id = "nongen_tabs",
            type = "pills",

            tabPanel(
              tr("Case LR"), value = "Case LR",
              br(),
              p(class = "interpretation-text",
                tr("One observed person compared with the MP profile.")),
              fluidRow(
                column(4,
                  wellPanel(
                    h5(tr("Observed person")),
                    checkboxInput("case_use_sex", tr("Use sex"), TRUE),
                    conditionalPanel(
                      condition = "input.case_use_sex",
                      selectInput("case_obsSex", tr("Observed sex:"), sex_choices, selectize = FALSE)
                    ),
                    checkboxInput("case_use_age", tr("Use age"), TRUE),
                    conditionalPanel(
                      condition = "input.case_use_age",
                      numericInput("case_obsAge", tr("Observed age (years):"), 42, 1, 80)
                    ),
                    checkboxInput("case_use_hair", tr("Use hair colour"), TRUE),
                    conditionalPanel(
                      condition = "input.case_use_hair",
                      selectInput("case_obsC", tr("Observed hair color:"), hair_choices,
                                  selected = 2, selectize = FALSE)
                    )
                  )
                ),
                column(8,
                  div(class = "result-box",
                    h4(tr("Likelihood Ratio")),
                    h2(textOutput("case_lr_value")),
                    p(textOutput("case_lr_log")),
                    p(class = "plain-sentence", textOutput("case_interp", inline = TRUE))
                  ),
                  div(style = "overflow-x: auto;", tableOutput("case_table")),
                  p(class = "text-muted small",
                    tr("Sex, age and hair colour are assumed conditionally independent, so their LRs multiply.")),
                  uiOutput("case_cell_note")
                )
              )
            ),

            tabPanel(
              tr("Profile performance"), value = "Profile performance",
              br(),
              div(class = "info-box",
                p(style = "margin: 0;",
                  tr("Not a case result. Every possible observation (sex, age inside or outside the MP range, hair colour) is scored against the MP profile. Weighting those LRs by P(D | H1) and by P(D | H2) gives the LRs to expect if the person is the MP or unrelated, and the error rates of a threshold."))
              ),
              plotOutput("cpt_plots", height = "auto"),
              p(class = "text-muted small",
                tr("log10(LR): blue supports H1, orange supports H2. Rows: F/M = sex; T1/T0 = age inside/outside the MP range."),
                " ", textOutput("cpt_case_note", inline = TRUE)),
              hr(),
              fluidRow(
                column(4,
                  wellPanel(
                    sliderInput("dec_threshold", tr("log10(LR) decision threshold:"),
                                min = -3, max = 3, value = 0, step = 0.1),
                    p(class = "text-muted small", textOutput("dec_threshold_lr", inline = TRUE)),
                    hr(),
                    verbatimTextOutput("dec_metrics"),
                    tags$details(
                      tags$summary(tr("Advanced metrics")),
                      p(class = "small", strong("MCC: "), textOutput("dec_mcc", inline = TRUE)),
                      p(class = "text-muted small", tr(mcc_help))
                    )
                  )
                ),
                column(8,
                  plotOutput("dist_plot", height = "350px")
                )
              ),
              fluidRow(
                column(6,
                  plotOutput("dec_roc_plot", height = "350px")
                ),
                column(6,
                  div(class = "success-box h1-box",
                    h5(tr("H1: the person is the MP (related)")),
                    verbatimTextOutput("dist_stats_h1")
                  ),
                  div(class = "warning-box h2-box",
                    h5(tr("H2: unrelated person")),
                    verbatimTextOutput("dist_stats_h2")
                  )
                )
              )
            ),

            tabPanel(
              tr("Database Search"), value = "Database Search",
              br(),
              p(tr("Non-genetic LR of each candidate against the MP profile (Marsico & Caridi 2023)."), " ",
                tr("Genetic false-positive rates:"), " ",
                actionLink("go_gen_rates_dbs", tr("Genetic Evidence > Rates at Threshold"))),
              fluidRow(
                column(4,
                  wellPanel(
                    h5(tr("Candidate Database")),
                    fileInput("dbs_file", NULL, accept = c(".csv", ".tsv", ".txt", "text/csv", "text/plain"),
                      buttonLabel = tr("Browse...")),
                    p(class = "text-muted small",
                      tr("CSV columns: ID, Sex (F/M), Age (1-80), HairColor (Black, Brown, Blonde, Red, Gray). Comma, semicolon or tab.")),
                    actionButton("dbs_load_example", tr("Load example"), class = "btn-primary btn-block"),
                    downloadButton("dbs_download_example", tr("Download example"),
                      class = "btn-default btn-block", style = "margin-top: 6px;"),
                    p(class = "text-muted small", style = "margin-top: 6px;",
                      tr("Synthetic example: 2,001 candidates; C0001 matches the default MP profile.")),
                    uiOutput("dbs_msg"),
                    hr(),
                    checkboxInput("dbs_empirical", tr("Estimate population frequencies from the database"), TRUE),
                    p(class = "text-muted small",
                      tr("Checked, candidate LRs can differ from Case LR, which uses the population values of the MP panel. Unchecked, they are the same."))
                  )
                ),
                column(8,
                  wellPanel(
                    sliderInput("dbs_threshold", tr("log10(LR) decision threshold:"),
                      min = -3, max = 3, value = 0, step = 0.1),
                    p(class = "text-muted small", textOutput("dbs_threshold_lr", inline = TRUE)),
                    p(class = "plain-sentence", textOutput("dbs_sentence", inline = TRUE)),
                    br(),
                    verbatimTextOutput("dbs_summary")
                  )
                )
              ),
              plotOutput("dbs_hist", height = "350px"),
              fluidRow(
                column(6, numericInput("dbs_top", tr("Top candidates shown:"), 20, 1, 500)),
                column(6, br(), downloadButton("dbs_download_results", tr("Download results")))
              ),
              div(style = "overflow-x: auto;", tableOutput("dbs_table"))
            )
          )
        )
      )
    ),

    tabPanel(
      tr("Genetic Evidence"), value = "Genetic Evidence",
      br(),
      tags$details(class = "info-box",
        tags$summary(tr("About")),
        p(style = "margin: 0;", tr("Genetic simulations are not run in the web version. Run them in Familias, forrel, mispitools (R) or FamLink and upload the results here."))
      ),
      fluidRow(
        column(4,
          wellPanel(
            div(class = "upload-zone",
              h5(tr("Simulation results")),
              fileInput("sim_file", NULL,
                accept = c(".csv", ".tsv", ".txt",
                           "text/csv", "text/tab-separated-values", "text/plain"),
                buttonLabel = tr("Browse...")),
              hr(),
              selectInput("sim_example", tr("Example dataset"), genetic_examples, selectize = FALSE),
              actionButton("sim_load_example", tr("Load example"),
                class = "btn-primary btn-block"),
              downloadButton("sim_download_example", tr("Download example"),
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
              h5(tr("Decision threshold")),
              numericInput("sim_threshold", tr("LR threshold (match if LR ≥ threshold):"),
                value = 100, min = 0, max = 1000000, step = 1),
              numericInput("sim_db_size", tr("Unrelated candidates compared (N):"),
                value = 10000, min = 1, step = 100),
              hr(),
              h5(tr("Optimal threshold")),
              sliderInput("sim_weight", tr("FP/FN cost ratio (×):"),
                min = 1, max = 50, value = 10, step = 1),
              actionButton("sim_calc_dt", tr("Compute optimal threshold"),
                class = "btn-primary btn-block"),
              uiOutput("sim_dt_sentence")
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
                  tr("Summary"), value = "Summary",
                  br(),
                  fluidRow(
                    column(6,
                      div(class = "metric-card",
                        span(class = "metric-label", tr("Simulations")),
                        div(class = "metric-value", textOutput("sim_n_sims"))
                      )
                    ),
                    column(6,
                      div(class = "metric-card",
                        span(class = "metric-label", tr("Relatives with LR ≥ 10,000")),
                        div(class = "metric-value", textOutput("sim_ip"))
                      )
                    )
                  ),
                  div(style = "overflow-x: auto;", tableOutput("sim_summary_table")),
                  tags$details(
                    tags$summary(tr("Advanced metrics")),
                    p(class = "small", strong(tr("Exclusion power: ")), textOutput("sim_ep", inline = TRUE)),
                    p(class = "text-muted small", tr("Share of unrelated simulations with LR = 0."))
                  ),
                  plotOutput("sim_violin", height = "350px")
                ),

                tabPanel(
                  tr("Distributions"), value = "Distributions",
                  br(),
                  plotOutput("sim_density_plot", height = "400px"),
                  fluidRow(
                    column(6,
                      div(class = "success-box h1-box",
                        h5(tr("H1: the person is the MP (related)"), "- log10(LR)"),
                        verbatimTextOutput("sim_dist_stats_rel")
                      )
                    ),
                    column(6,
                      div(class = "warning-box h2-box",
                        h5(tr("H2: unrelated person"), "- log10(LR)"),
                        verbatimTextOutput("sim_dist_stats_unrel")
                      )
                    )
                  )
                ),

                tabPanel(
                  tr("Decision Plot"), value = "Decision Plot",
                  br(),
                  p(class = "text-muted small",
                    tr("Each point is an LR threshold; grey labels give its value. Points nearer (0, 0) are better.")),
                  plotOutput("sim_decision_plot", height = "450px")
                ),

                tabPanel(
                  tr("Rates at Threshold"), value = "Rates at Threshold",
                  br(),
                  div(class = "info-box",
                    p(class = "plain-sentence", style = "margin: 0;", textOutput("sim_rates_sentence", inline = TRUE))
                  ),
                  fluidRow(
                    column(3,
                      div(class = "metric-card",
                        span(class = "metric-label", tr("True Positive Rate")),
                        div(class = "metric-value", style = "font-size: 18px;", textOutput("sim_rate_tpr"))
                      )
                    ),
                    column(3,
                      div(class = "metric-card",
                        span(class = "metric-label", tr("False Positive Rate")),
                        div(class = "metric-value", style = "font-size: 18px;", textOutput("sim_rate_fpr"))
                      )
                    ),
                    column(3,
                      div(class = "metric-card",
                        span(class = "metric-label", tr("False Negative Rate")),
                        div(class = "metric-value", style = "font-size: 18px;", textOutput("sim_rate_fnr"))
                      )
                    ),
                    column(3,
                      div(class = "metric-card",
                        span(class = "metric-label", tr("True Negative Rate")),
                        div(class = "metric-value", style = "font-size: 18px;", textOutput("sim_rate_tnr"))
                      )
                    )
                  ),
                  tags$details(
                    tags$summary(tr("Advanced metrics")),
                    uiOutput("sim_advanced")
                  )
                )
              )
            )
          ),
          conditionalPanel(
            condition = "!output.sim_has_data",
            div(class = "info-box",
              h5(tr("Accepted formats")),
              tags$ul(
                tags$li(tags$b(tr("Familias, one family:")), " ",
                  tr("tab-separated; two header rows (case, pedigree label); column 1 Unrelated, column 2 Related.")),
                tags$li(tags$b(tr("Familias, several families:")), " ",
                  tr("one column pair per family; row 2 holds 'Reference pedigree' (H2) and a relationship label (H1).")),
                tags$li(tags$b("FamLink:"), " ",
                  tr("'Raw data LR' export with 'Pedigree (Numerator)', 'Versus (Denominator)' and 'True:' rows; groups separated by '|'.")),
                tags$li(tags$b(tr("Generic CSV/TSV:")), " ",
                  tr("columns 'Unrelated' and 'Related', or the first two numeric columns; comma, tab or semicolon."))
              ),
              p(class = "text-muted small", style = "margin: 0;",
                tr("Text files only (.txt, .csv, .tsv); Excel is not supported."))
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {

  session$userData$shiny.i18n <- list(lang = reactiveVal("en"))

  select_choices <- list(
    ng_MPs = sex_choices, case_obsSex = sex_choices,
    ng_MPc = hair_choices, case_obsC = hair_choices,
    sim_example = genetic_examples
  )

  observeEvent(input$lang, {
    req(input$lang %in% c("en", "es"))
    session$userData$shiny.i18n$lang(input$lang)
    for (id in names(select_choices)) {
      updateSelectInput(session, id, choices = tr_names(select_choices[[id]]),
                        selected = input[[id]])
    }
    session$sendCustomMessage("i18n-placeholder", tr("No file selected"))
  }, priority = 10)

  go_to <- function(main, sub = NULL, sim = NULL) {
    updateTabsetPanel(session, "main_tabs", selected = main)
    if (!is.null(sub)) updateTabsetPanel(session, "nongen_tabs", selected = sub)
    if (!is.null(sim)) updateTabsetPanel(session, "sim_subtabs", selected = sim)
  }

  observeEvent(input$go_case, go_to("Non-Genetic Evidence", "Case LR"))
  observeEvent(input$go_dbs, go_to("Non-Genetic Evidence", "Database Search"))
  observeEvent(input$go_choose, go_to("Non-Genetic Evidence", "Profile performance"))
  observeEvent(input$go_gen_rates, go_to("Genetic Evidence", sim = "Rates at Threshold"))
  observeEvent(input$go_gen_rates_dbs, go_to("Genetic Evidence", sim = "Rates at Threshold"))
  observeEvent(input$go_gen_optimal, go_to("Genetic Evidence", sim = "Decision Plot"))

  interpretation <- function(LR, evidence) {
    if (LR == 0) return(tr("Impossible under H1 (exclusion)."))
    if (LR == 1) return(tr("Equally probable under H1 and H2."))
    trf(interpretation_templates[[evidence]][if (LR > 1) 1 else 2],
        format(signif(if (LR > 1) LR else 1 / LR, 3), big.mark = big_mark(), scientific = FALSE))
  }

  threshold_lr_text <- function(t) {
    trf("Match if %s", sprintf("log10(LR) \u2265 %.1f (LR \u2265 %s)", t,
                                format(signif(10^t, 4), big.mark = big_mark(), scientific = FALSE)))
  }

  output$dec_threshold_lr <- renderText(threshold_lr_text(input$dec_threshold))
  output$dbs_threshold_lr <- renderText(threshold_lr_text(input$dbs_threshold))

  is_narrow <- function(id) {
    w <- session$clientData[[paste0("output_", id, "_width")]]
    !is.null(w) && w < 600
  }

  plot_size <- function(id, wide) {
    if (is_narrow(id)) 10 else wide
  }

  ng_profile <- reactive({
    validate(
      check_num(input$ng_MPa, "Age", 1, 80, whole = TRUE),
      check_num(input$ng_MPr, "Age uncertainty", 1, 20, whole = TRUE),
      check_num(input$ng_eps, "Sex error rate", 0, 0.5),
      check_num(input$ng_epa, "Age error rate", 0, 0.5),
      check_num(input$ng_epc, "Hair error rate", 0, 0.5),
      check_num(input$ng_propF, "Female proportion", 0.3, 0.7)
    )
    list(MPs = input$ng_MPs, MPa = input$ng_MPa, MPr = input$ng_MPr, MPc = as.numeric(input$ng_MPc),
         eps = input$ng_eps, epa = input$ng_epa, epc = input$ng_epc, propF = input$ng_propF)
  })

  case_parts <- reactive({
    p <- ng_profile()
    aw <- age_window(p$MPa, p$MPr)
    sex_name <- function(s) tr(names(sex_choices)[sex_choices == s])
    hair_name <- function(h) tr(names(hair_choices)[h])
    rows <- list()
    cell <- list()
    if (isTRUE(input$case_use_sex)) {
      obs <- input$case_obsSex
      rows[[length(rows) + 1]] <- data.frame(
        ev = tr("Sex"), mp = sex_name(p$MPs), obs = sex_name(obs),
        h1 = if (obs == p$MPs) 1 - p$eps else p$eps,
        h2 = if (obs == "F") p$propF else 1 - p$propF)
      cell$sex <- obs
    }
    if (isTRUE(input$case_use_age)) {
      validate(check_num(input$case_obsAge, "Observed age", 1, 80, whole = TRUE))
      obs <- input$case_obsAge
      inside <- obs >= aw$min && obs <= aw$max
      rows[[length(rows) + 1]] <- data.frame(
        ev = tr("Age"), mp = sprintf("%d-%d", aw$min, aw$max),
        obs = sprintf("%d (%s)", obs, tr(if (inside) "in range" else "out of range")),
        h1 = if (inside) 1 - p$epa else p$epa,
        h2 = if (inside) aw$T1p else 1 - aw$T1p)
      cell$age <- if (inside) "T1" else "T0"
    }
    if (isTRUE(input$case_use_hair)) {
      obs <- as.numeric(input$case_obsC)
      rows[[length(rows) + 1]] <- data.frame(
        ev = tr("Hair colour"), mp = hair_name(p$MPc), obs = hair_name(obs),
        h1 = if (obs == p$MPc) 1 - p$epc else p$epc / (length(DEFAULT_HAIR_FREQ) - 1),
        h2 = unname(DEFAULT_HAIR_FREQ[obs]))
      cell$hair <- names(DEFAULT_HAIR_FREQ)[obs]
    }
    tab <- do.call(rbind, rows)
    if (!is.null(tab)) tab$lr <- tab$h1 / tab$h2
    full <- length(cell) == 3
    list(table = tab,
         LR = if (is.null(tab)) NA else prod(tab$lr),
         cell = if (full) c(paste(cell$sex, cell$age, sep = "-"), cell$hair) else NULL)
  })

  output$case_lr_value <- renderText({
    LR <- case_parts()$LR
    if (is.na(LR)) tr("No evidence selected") else format_lr(LR)
  })

  output$case_lr_log <- renderText({
    LR <- case_parts()$LR
    if (is.na(LR)) "" else format_log_lr(LR)
  })

  output$case_interp <- renderText({
    cp <- case_parts()
    if (is.na(cp$LR)) return("")
    interpretation(cp$LR, if (nrow(cp$table) == 1) "single" else "combined")
  })

  output$case_table <- renderTable({
    tab <- case_parts()$table
    req(tab)
    num <- function(x) sprintf("%.4g", x)
    out <- data.frame(tab$ev, tab$mp, tab$obs, num(tab$h1), num(tab$h2), num(tab$lr),
                      check.names = FALSE)
    names(out) <- c(tr("Evidence"), tr("MP"), tr("Observed"),
                    "P(obs | H1)", "P(obs | H2)", "LR")
    out
  }, striped = TRUE, align = "lllrrr")

  output$case_cell_note <- renderUI({
    if (is.null(case_parts()$cell)) {
      p(class = "text-muted small",
        tr("With evidence left out, this LR is not a cell of the Profile performance table."))
    } else {
      p(class = "text-muted small",
        tr("This LR is the outlined cell of the table in"), " ",
        actionLink("goto_performance", tr("Profile performance")))
    }
  })

  observeEvent(input$goto_performance, go_to("Non-Genetic Evidence", "Profile performance"))

  cpt_pop <- reactive({
    p <- ng_profile()
    build_cpt_pop(p$propF, p$MPa, p$MPr, DEFAULT_HAIR_FREQ)
  })

  cpt_mp <- reactive({
    p <- ng_profile()
    build_cpt_mp(p$MPs, p$MPc, p$eps, p$epa, p$epc)
  })

  case_log_lr <- reactive({
    cp <- case_parts()
    if (is.null(cp$cell) || cp$LR == 0) NULL else log10(cp$LR)
  })

  output$cpt_case_note <- renderText({
    if (is.null(case_parts()$cell)) "" else tr("Outlined cell: the observation entered in Case LR.")
  })

  output$cpt_plots <- renderPlot({
    create_cpt_plots(cpt_pop(), cpt_mp(), narrow = is_narrow("cpt_plots"),
                     case_cell = case_parts()$cell)
  }, res = 100, height = function() if (is_narrow("cpt_plots")) 1000 else 380,
  alt = function() tr("Heatmaps of P(D | H2), P(D | H1) and log10(LR) for each sex, age class and hair colour of the missing person profile."))

  lr_dists <- reactive({
    build_lr_distributions(cpt_pop(), cpt_mp())
  })

  output$dist_plot <- renderPlot({
    dists <- lr_dists()
    create_dist_plot(dists$H1, dists$H2, base_size = plot_size("dist_plot", 12),
                     threshold = input$dec_threshold, case_log = case_log_lr())
  }, res = 100, alt = function() one_line(dist_title(lr_dists()$H1, lr_dists()$H2)))

  roc_alt <- function() {
    one_line(roc_title(build_roc(lr_dists()$H1, lr_dists()$H2)$AUC))
  }

  dist_stats <- function(d) {
    if (nrow(d) == 0) {
      cat(tr("No data available"), "\n")
    } else {
      p0 <- sum(d$w[d$x == -Inf])
      xf <- d$x[is.finite(d$x)]
      if (p0 > 0) {
        cat(tr("Weighted mean log10(LR): undefined, P(LR = 0) ="), sprintf("%.3f", p0), "\n")
        cat(tr("Range: LR = 0 to log10(LR) ="), sprintf("%.3f", max(xf)), "\n")
      } else {
        cat(tr("Weighted mean log10(LR):"), sprintf("%.3f", sum(d$x * d$w)), "\n")
        cat(tr("Range:"), "[", sprintf("%.3f", min(xf)), ",", sprintf("%.3f", max(xf)), "]\n")
      }
    }
  }

  output$dist_stats_h1 <- renderPrint({
    dist_stats(lr_dists()$H1)
  })

  output$dist_stats_h2 <- renderPrint({
    dist_stats(lr_dists()$H2)
  })

  output$dec_metrics <- renderPrint({
    dists <- lr_dists()
    m <- calc_metrics(dists$H1, dists$H2, input$dec_threshold)

    rows <- c(tr("True Positive Rate:"), tr("False Positive Rate:"),
              tr("False Negative Rate:"), tr("True Negative Rate:"))
    width <- max(nchar(rows))
    cat(sprintf("%s %.3f\n", formatC(rows, width = -width), c(m$TPR, m$FPR, m$FNR, m$TNR)), sep = "")
  })

  output$dec_mcc <- renderText({
    dists <- lr_dists()
    sprintf("%.3f", calc_metrics(dists$H1, dists$H2, input$dec_threshold)$MCC)
  })

  output$dec_roc_plot <- renderPlot({
    dists <- lr_dists()
    roc_obj <- build_roc(dists$H1, dists$H2)
    df_roc <- roc_obj$data
    auc_val <- roc_obj$AUC

    if (nrow(df_roc) == 0) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5, label = tr("Insufficient data")) +
               ggplot2::theme_void())
    }

    m <- calc_metrics(dists$H1, dists$H2, input$dec_threshold)
    current_point <- data.frame(FPR = m$FPR, TPR = m$TPR)

    ggplot2::ggplot(df_roc, ggplot2::aes(x = FPR, y = TPR)) +
      ggplot2::geom_path(color = "grey35", linewidth = 1) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                           color = COL_REFERENCE_LINE, linewidth = 0.8) +
      ggplot2::geom_point(data = current_point, ggplot2::aes(x = FPR, y = TPR),
                         color = COL_THRESHOLD, size = 4.5, shape = 16) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(title = roc_title(auc_val),
                    subtitle = trf("Black point: selected threshold (log10(LR) \u2265 %.1f)",
                                   input$dec_threshold),
                    x = tr("False Positive Rate"), y = tr("True Positive Rate")) +
      ggplot2::theme_minimal(base_size = plot_size("dec_roc_plot", 14)) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }, res = 100, alt = function() paste(roc_alt(),
    trf("Black point: selected threshold (log10(LR) \u2265 %.1f)", input$dec_threshold)))

  sim_families_val <- reactiveVal(NULL)
  sim_parse_error  <- reactiveVal(NULL)

  observeEvent(input$sim_file, {
    result <- tryCatch(parse_sim_upload(input$sim_file$datapath), error = function(e) NULL)
    if (is.null(result)) {
      sim_families_val(NULL)
      sim_parse_error(c("Could not read '%s'. See the accepted formats below.", input$sim_file$name))
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
               alt = names(genetic_examples)[genetic_examples == input$sim_example],
               style = "width: 100%; height: auto; border: 1px solid #ddd; background: white;"),
      div(class = "ped-legend",
        span(class = "ped-item", span(class = "swatch swatch-ref"), tr_span("Genotyped reference (Ref)")),
        span(class = "ped-item", span(class = "swatch swatch-mp"), tr_span("Missing person (labelled MP)")),
        span(class = "ped-item", span(class = "swatch"), tr_span("Not genotyped"))))
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
      h5(tr_span("Select Comparison")),
      selectInput("sim_family", NULL, choices = names(fams), selected = names(fams)[1]),
      p(class = "text-muted small",
        length(fams), tr_span("comparisons detected in this file."))
    )
  })

  output$sim_parse_msg <- renderUI({
    err <- sim_parse_error()
    if (!is.null(err)) {
      err_lines <- strsplit(fmt_msg(err), "\n")[[1]]
      return(div(class = "warning-box", style = "margin-top: 10px; padding: 10px;",
        tags$b(tr("File format error")),
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
      trf("Data loaded: %s simulations (%d comparisons available)",
          format(nrow(df), big.mark = big_mark()), n_fams)
    } else {
      trf("Data loaded: %s simulations", format(nrow(df), big.mark = big_mark()))
    }
  })

  output$sim_n_sims <- renderText({
    df <- sim_data()
    req(df)
    format(nrow(df), big.mark = big_mark())
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

  output$sim_summary_table <- renderTable({
    df <- sim_data()
    req(df)
    stats_col <- function(x) {
      q <- stats::quantile(x, c(0.05, 0.95), names = FALSE)
      c(format_lr_value(stats::median(x)),
        paste(format_lr_value(q[1]), "\u2013", format_lr_value(q[2])),
        sprintf("%.1f%%", 100 * mean(x == 0)))
    }
    out <- data.frame(c(tr("Median LR"), tr("5th\u201395th percentile of LR"), tr("Simulations with LR = 0")),
                      stats_col(df$Related), stats_col(df$Unrelated))
    names(out) <- c("", tr("H1: the person is the MP (related)"), tr("H2: unrelated person"))
    out
  }, striped = TRUE, bordered = TRUE, align = "lrr")

  zero_caption <- function(df) {
    if (any(df$Related == 0) || any(df$Unrelated == 0)) {
      trf("Not shown (LR = 0): %d related and %d unrelated simulations",
          sum(df$Related == 0), sum(df$Unrelated == 0))
    } else {
      NULL
    }
  }

  sim_violin_title <- reactive({
    df <- sim_data()
    req(df)
    trf("Median log10(LR): %s under H1 (related),\n%s under H2 (unrelated)",
        format_log_value(stats::median(log10(df$Related))),
        format_log_value(stats::median(log10(df$Unrelated))))
  })

  sim_density_title <- reactive({
    df <- sim_data()
    req(df)
    trf("LR > 1 in %.1f%% of related and %.1f%% of unrelated simulations",
        100 * mean(df$Related > 1), 100 * mean(df$Unrelated > 1))
  })

  output$sim_violin <- renderPlot({
    df <- sim_data()
    req(df)

    rel_valid <- df$Related[df$Related > 0]
    unrel_valid <- df$Unrelated[df$Unrelated > 0]

    lab_h1 <- tr("H1: the person is the MP (related)")
    lab_h2 <- tr("H2: unrelated person")
    df_long <- data.frame(
      Hypothesis = factor(c(
        rep(lab_h1, length(rel_valid)),
        rep(lab_h2, length(unrel_valid))
      ), levels = c(lab_h1, lab_h2)),
      log10LR = c(log10(rel_valid), log10(unrel_valid))
    )

    ggplot2::ggplot(df_long, ggplot2::aes(x = Hypothesis, y = log10LR, fill = Hypothesis)) +
      ggplot2::geom_violin(alpha = 0.7, color = "white", linewidth = 0.3) +
      ggplot2::geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, color = "grey30",
                            outlier.size = 0.5, outlier.alpha = 0.3) +
      ggplot2::scale_fill_manual(values = setNames(c(COL_H1, COL_H2), c(lab_h1, lab_h2))) +
      ggplot2::geom_hline(yintercept = 0, color = COL_REFERENCE_LINE, linewidth = 0.8) +
      ggplot2::labs(title = sim_violin_title(),
                    subtitle = tr("Grey line: LR = 1 (neutral evidence)"),
                    caption = zero_caption(df),
                    x = "", y = log_axis_label()) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = plot_size("sim_violin", 13)) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank()
      )
  }, res = 100, alt = function() one_line(sim_violin_title()))

  output$sim_density_plot <- renderPlot({
    df <- sim_data()
    req(df)

    rel_valid <- df$Related[df$Related > 0]
    unrel_valid <- df$Unrelated[df$Unrelated > 0]

    if (length(rel_valid) < 2 || length(unrel_valid) < 2) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5, label = tr("Insufficient positive LR values")) +
               ggplot2::theme_void())
    }

    dens_rel <- stats::density(log10(rel_valid))
    dens_unrel <- stats::density(log10(unrel_valid))

    lab_h1 <- tr("H1: the person is the MP (related)")
    lab_h2 <- tr("H2: unrelated person")
    df_dens <- rbind(
      data.frame(x = dens_rel$x, y = dens_rel$y, Hypothesis = lab_h1),
      data.frame(x = dens_unrel$x, y = dens_unrel$y, Hypothesis = lab_h2)
    )

    ggplot2::ggplot(df_dens, ggplot2::aes(x = x, y = y, fill = Hypothesis)) +
      ggplot2::geom_area(alpha = 0.55, position = "identity", color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = setNames(c(COL_H1, COL_H2), c(lab_h1, lab_h2)),
                                 breaks = c(lab_h1, lab_h2)) +
      ggplot2::geom_vline(xintercept = 0, color = COL_REFERENCE_LINE, linewidth = 0.8) +
      ggplot2::labs(title = sim_density_title(),
                    subtitle = tr("Grey line: LR = 1."),
                    caption = zero_caption(df),
                    x = log_axis_label(), y = tr("Density")) +
      ggplot2::theme_minimal(base_size = plot_size("sim_density_plot", 13)) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  }, res = 100, alt = function() sim_density_title())

  sim_log_stats <- function(x) {
    valid <- x[x > 0]
    if (length(valid) == 0) {
      cat(tr("No positive values"), "\n")
    } else {
      lv <- log10(valid)
      rows <- formatC(c(tr("Mean:"), tr("Median:"), tr("SD:"), tr("Range:")),
                      width = -max(nchar(c(tr("Mean:"), tr("Median:"), tr("SD:"), tr("Range:")))))
      cat(rows[1], " ", sprintf("%.3f", mean(lv)), "\n", sep = "")
      cat(rows[2], " ", sprintf("%.3f", stats::median(lv)), "\n", sep = "")
      cat(rows[3], " ", sprintf("%.3f", stats::sd(lv)), "\n", sep = "")
      cat(rows[4], " [ ", sprintf("%.3f", min(lv)), " , ", sprintf("%.3f", max(lv)), " ]\n", sep = "")
    }
    if (any(x == 0)) {
      cat(trf("Excluded (LR = 0): %d of %d", sum(x == 0), length(x)), "\n")
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
        low = "#e0e0e0", mid = "#a6a6a6", high = "#595959",
        midpoint = stats::median(log10(curve_data$LR)),
        name = bquote(log[10](LR) ~ .(tr("threshold")))
      ) +
      ggplot2::geom_point(data = opt_point, ggplot2::aes(x = FNR, y = FPR),
                         color = COL_OPTIMAL, size = 5.5, shape = 18, inherit.aes = FALSE) +
      ggplot2::annotate("text",
        x = opt_point$FNR, y = opt_point$FPR,
        label = trf("Optimal threshold\nLR = %.2f", opt_point$LR),
        hjust = 0.5, vjust = -0.6, size = 3.5, fontface = "bold", color = COL_OPTIMAL) +
      ggplot2::geom_point(data = sel_point, ggplot2::aes(x = FNR, y = FPR),
                         color = COL_THRESHOLD, size = 4.5, shape = 16, inherit.aes = FALSE) +
      ggplot2::annotate("text",
        x = sel_point$FNR, y = sel_point$FPR,
        label = trf("Selected\nLR = %g", input$sim_threshold),
        hjust = 0.5, vjust = 1.5, size = 3.2, color = COL_THRESHOLD) +
      ggplot2::geom_text(data = label_data,
        ggplot2::aes(x = FNR, y = FPR, label = round(LR)),
        size = 2.8, vjust = -1, color = "gray40", inherit.aes = FALSE,
        check_overlap = TRUE) +
      ggplot2::scale_x_continuous(n.breaks = if (is_narrow("sim_decision_plot")) 4 else 6) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.15))) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(
        title = sim_decision_title(),
        subtitle = if (is_narrow("sim_decision_plot")) {
          tr("Green diamond: optimal threshold\nBlack dot: selected threshold")
        } else {
          tr("Green diamond: optimal threshold. Black dot: selected threshold.")
        },
        x = tr("False Negative Rate (FNR)"),
        y = tr("False Positive Rate (FPR)")
      ) +
      ggplot2::theme_minimal(base_size = plot_size("sim_decision_plot", 13)) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = if (is_narrow("sim_decision_plot")) "bottom" else "right",
        panel.grid.minor = ggplot2::element_blank()
      )
  }, res = 100, alt = function() sim_decision_title())

  sim_decision_title <- reactive({
    trf("Optimal decision threshold: LR \u2265 %.2f (weight %s)", sim_opt()$DT, input$sim_weight)
  })

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

  rate_text <- function(count) {
    n <- nrow(sim_data())
    trf("%s (%s of %s simulations)", sprintf("%.1f%%", 100 * count / n),
        format(count, big.mark = big_mark()), format(n, big.mark = big_mark()))
  }

  output$sim_rate_tpr <- renderText(rate_text(sim_rates()$TP_count))
  output$sim_rate_fpr <- renderText(rate_text(sim_rates()$FP_count))
  output$sim_rate_fnr <- renderText(rate_text(sim_rates()$FN_count))
  output$sim_rate_tnr <- renderText(rate_text(sim_rates()$TN_count))

  output$sim_rates_sentence <- renderText({
    r <- sim_rates()
    validate(sim_n_message())
    ppv <- if (is.na(sim_ppv())) tr("undefined") else sprintf("%.1f%%", 100 * sim_ppv())
    trf("LR \u2265 %s, N = %s: %s expected false positives; true relative detected in %s; PPV %s.",
        format(input$sim_threshold, big.mark = big_mark(), scientific = FALSE),
        format(sim_n_candidates(), big.mark = big_mark(), scientific = FALSE),
        sprintf("%.1f", r$FPR * sim_n_candidates()),
        sprintf("%.1f%%", 100 * r$TPR), ppv)
  })

  output$sim_advanced <- renderUI({
    r <- sim_rates()
    validate(sim_n_message())
    item <- function(label, value, help) {
      tagList(p(class = "small", style = "margin: 6px 0 0 0;", strong(label), " ", value),
              p(class = "text-muted small", style = "margin: 0;", help))
    }
    tagList(
      item(tr("PPV:"), if (is.na(sim_ppv())) tr("undefined") else sprintf("%.4f", sim_ppv()),
           tr("Share of reported matches that are the true relative, with one true relative among N.")),
      item(tr("Accuracy:"), sprintf("%.4f", r$Accuracy),
           tr("Share of correct decisions; balanced classes.")),
      item("MCC:", sprintf("%.4f", r$MCC), tr(mcc_help))
    )
  })

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

  output$sim_dt_sentence <- renderUI({
    req(sim_dt_shown())
    result <- sim_opt()
    validate(sim_n_message())
    div(class = "info-box", style = "margin-bottom: 0;",
      p(style = "margin: 0 0 4px 0;",
        trf("Optimal: LR \u2265 %.2f; FPR %s, FNR %s; %s expected false positives per %s.",
            result$DT, sprintf("%.4f", result$FPR), sprintf("%.4f", result$FNR),
            sprintf("%.1f", result$FPR * sim_n_candidates()),
            format(sim_n_candidates(), big.mark = big_mark(), scientific = FALSE))),
      p(class = "text-muted small", style = "margin: 0 0 8px 0;", "D = sqrt(FNR^2 + (w \u00d7 FPR)^2)"),
      actionButton("sim_use_dt", tr("Use this threshold"), class = "btn-default btn-block"))
  })

  observeEvent(input$sim_use_dt, {
    updateNumericInput(session, "sim_threshold", value = sim_opt()$DT)
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
      div(class = "warning-box", style = "padding: 10px;", tags$b(tr("File error: ")), fmt_msg(err))
    } else if (!is.null(dbs_data())) {
      div(class = "success-box", style = "padding: 10px;",
          trf("Loaded %s candidates.", format(nrow(dbs_data()), big.mark = big_mark())))
    }
  })

  dbs_lr <- reactive({
    cands <- dbs_data()
    validate(need(!is.null(cands), tr("Upload a candidate database or load the example.")))
    p <- ng_profile()
    candidate_lrs(cands, p$MPs, p$MPa, p$MPr, p$MPc, p$eps, p$epa, p$epc, p$propF,
                  isTRUE(input$dbs_empirical))
  })

  dbs_above_title <- reactive({
    res <- dbs_lr()$results
    trf("%s of %s candidates reach LR \u2265 %s",
        format(sum(res$log10LR >= input$dbs_threshold), big.mark = big_mark()),
        format(nrow(res), big.mark = big_mark()),
        format(signif(10^input$dbs_threshold, 4), big.mark = big_mark(), scientific = FALSE))
  })

  output$dbs_sentence <- renderText({
    res <- dbs_lr()$results
    above <- sum(res$log10LR >= input$dbs_threshold)
    paste(trf("%s of %s candidates (%s) reach LR \u2265 %s.",
              format(above, big.mark = big_mark()), format(nrow(res), big.mark = big_mark()),
              sprintf("%.1f%%", 100 * above / nrow(res)),
              format(signif(10^input$dbs_threshold, 4), big.mark = big_mark(), scientific = FALSE)),
          tr("Candidates in the same sex, age and hair class share one LR."))
  })

  output$dbs_summary <- renderPrint({
    out <- dbs_lr()
    res <- out$results
    n <- nrow(res)
    cat(tr("Candidates:"), n, "\n")
    cat(tr("LR = 0 (exclusions):"), sum(res$LR == 0), "\n\n")
    cat(tr(if (isTRUE(input$dbs_empirical)) "Frequencies estimated from the database:" else
          "Default population frequencies:"), "\n")
    cat(sprintf("  P(F) = %.4f\n", out$propF))
    cat(trf("  P(age %g-%g | H2) = %.4f", out$age_min, out$age_max, out$T1p), "\n")
    cat(" ", tr("Hair:"), paste(sprintf("%s %.4f", vapply(names(out$hair_freq), tr, character(1)),
                                       out$hair_freq), collapse = ", "), "\n")
    if (length(out$notes) > 0) {
      cat("\n", paste(out$notes, collapse = "\n"), "\n", sep = "")
    }
  })

  output$dbs_hist <- renderPlot({
    res <- dbs_lr()$results
    res$w <- 1 / nrow(res)
    fin <- res[is.finite(res$log10LR), ]
    p0 <- sum(res$LR == 0) / nrow(res)
    binwidth <- 0.25
    validate(need(nrow(fin) > 0, tr("All candidates have LR = 0.")))
    p <- ggplot2::ggplot(fin, ggplot2::aes(x = log10LR, weight = w)) +
      ggplot2::geom_histogram(binwidth = binwidth, fill = COL_NEUTRAL, alpha = 0.8,
                              color = "white", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = input$dbs_threshold, linetype = "dashed",
                          color = COL_THRESHOLD, linewidth = 0.8) +
      ggplot2::labs(title = tr("Candidate log10(LR)"),
                    subtitle = tr("Black dashed line: decision threshold"),
                    x = log_axis_label(), y = tr("Probability")) +
      ggplot2::theme_minimal(base_size = plot_size("dbs_hist", 12)) +
      ggplot2::theme(plot.title.position = "plot",
                     plot.title = ggplot2::element_text(face = "bold"),
                     panel.grid.minor = ggplot2::element_blank())
    if (p0 > 0) {
      zx <- min(fin$log10LR) - 3 * binwidth
      zero <- data.frame(log10LR = zx, w = p0)
      p <- p +
        ggplot2::geom_col(data = zero, ggplot2::aes(x = log10LR, y = w), inherit.aes = FALSE,
                          width = binwidth, fill = COL_NEUTRAL, alpha = 0.8) +
        ggplot2::geom_text(data = zero, ggplot2::aes(x = log10LR, y = w,
                           label = sprintf("LR = 0\nP = %.3f", w)),
                           inherit.aes = FALSE, vjust = -0.3, size = 3.2) +
        ggplot2::geom_vline(xintercept = zx + 1.5 * binwidth, linetype = "dotted", color = "gray50") +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.2))) +
        ggplot2::labs(subtitle = tr("Black dashed line: decision threshold. Left bar: LR = 0 (exclusion)"))
    }
    p
  }, res = 100, alt = function() dbs_above_title())

  output$dbs_table <- renderTable({
    res <- dbs_lr()$results
    validate(check_num(input$dbs_top, "Top candidates shown", 1, 500, whole = TRUE))
    top <- utils::head(res, input$dbs_top)
    out <- data.frame(
      Rank = seq_len(nrow(top)),
      ID = top$ID, Sex = top$Sex, Age = as.integer(top$Age),
      HairColor = vapply(top$HairColor, tr, character(1), USE.NAMES = FALSE),
      LR = sub("LR = ", "", vapply(top$LR, format_lr, character(1)), fixed = TRUE),
      log10LR = ifelse(top$LR == 0, tr("LR = 0 (exclusion)"), sprintf("%.3f", top$log10LR)),
      check.names = FALSE
    )
    names(out) <- c(tr("Rank"), "ID", tr("Sex"), tr("Age"), tr("Hair color"), "LR", "log10(LR)")
    out
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
