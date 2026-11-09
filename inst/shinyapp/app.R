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

hair_error_help <- "Total probability of observing a wrong color, split equally among the other 4 colors."

hair_freq_help <- "Black: 8%, Brown: 45%, Blonde: 30%, Red: 7%, Gray: 10%"

defaults_help <- "Defaults from Marsico & Caridi 2023; change only if your lab has its own estimates."

error_help <- "Probability that the recorded value is wrong, e.g. 0.05 = 1 in 20 records."

interpretation_templates <- list(
  sex = c("The observed sex is %s times more probable if this person is the missing person than if they are an unrelated person from the population.",
          "The observed sex is %s times more probable if they are an unrelated person from the population than if this person is the missing person."),
  age = c("The observed age is %s times more probable if this person is the missing person than if they are an unrelated person from the population.",
          "The observed age is %s times more probable if they are an unrelated person from the population than if this person is the missing person."),
  hair = c("The observed hair colour is %s times more probable if this person is the missing person than if they are an unrelated person from the population.",
           "The observed hair colour is %s times more probable if they are an unrelated person from the population than if this person is the missing person."),
  combined = c("The combined observations are %s times more probable if this person is the missing person than if they are an unrelated person from the population, assuming sex, age and hair colour are independent.",
               "The combined observations are %s times more probable if they are an unrelated person from the population than if this person is the missing person, assuming sex, age and hair colour are independent.")
)

mcc_help <- "Matthews correlation coefficient: agreement between decisions and truth, from -1 to 1, assuming balanced classes."

check_num <- function(x, label, lo, hi, whole = FALSE) {
  if (is.null(x) || length(x) != 1 || !is.numeric(x) || !is.finite(x)) {
    return(trf("%s: please enter a number.", tr(label)))
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

create_dist_plot <- function(distH1, distH2, binwidth = 0.25, base_size = 12) {
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
    tr("Familias example dataset, 10,000 simulations.")
  } else {
    tr("Simulated with forrel 1.9.1: 10,000 replicates under H1 (the person of interest is the missing person) and H2 (unrelated), 24 autosomal STRs, Argentina allele frequencies, no mutation model. References are genotyped relatives of the missing person.")
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
      notes <- c(notes, tr("Only one sex is present in the database: the default female proportion (0.51) is used."))
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
          h2(tr("Evidence Evaluation for Kinship in Forensic Science")),
          p(class = "interpretation-text",
            tr("This application implements a Bayesian framework for evaluating evidence in kinship cases for forensic science. It includes two main sections: non-genetic evidence (sex, age, hair color) with simulation capabilities, and genetic evidence analysis (upload pre-computed simulation results).")
          ),
          hr(),

          fluidRow(
            column(6,
              div(class = "info-box",
                h4(tr("The Likelihood Ratio")),
                p(tr("The LR compares two hypotheses:")),
                tags$ul(
                  tags$li(strong("H1:"), " ", tr("the person is the missing person (MP), related as stated in the pedigree")),
                  tags$li(strong("H2:"), " ", tr("unrelated person"))
                ),
                div(class = "formula-box",
                  tr("LR = P(Evidence | H1) / P(Evidence | H2)")
                ),
                p(tr("Example: LR = 100 means the evidence is 100 times more probable if related than if unrelated.")),
                p(tr("The LR does "), strong(tr("not")), tr(" provide a probability of identity. It quantifies how much the evidence changes the odds."))
              )
            ),
            column(6,
              div(class = "info-box",
                h4(tr("Using the LR")),
                p(tr("The LR quantifies how many times more likely the evidence is under H1 vs H2.")),
                tags$ul(
                  tags$li(tr("LR > 1: Evidence more probable under H1")),
                  tags$li(tr("LR < 1: Evidence more probable under H2")),
                  tags$li(tr("LR = 1: Evidence equally probable under both hypotheses"))
                ),
                p(tr("The numeric LR value should be reported without verbal qualifiers."))
              )
            )
          ),

          h4(tr("Where to start")),
          fluidRow(
            column(4,
              div(class = "well entry-card",
                h5(tr("How strong is an observed sex/age/hair match?")),
                actionLink("go_individual", tr("Non-genetic: Individual Evidence")),
                actionLink("go_combine", tr("Non-genetic: Combine Evidence"))
              )
            ),
            column(4,
              div(class = "well entry-card",
                h5(tr("How many false matches will a database search give?")),
                actionLink("go_gen_rates", tr("Genetic: Rates at Threshold")),
                actionLink("go_dbs", tr("Non-genetic: Database Search"))
              )
            ),
            column(4,
              div(class = "well entry-card",
                h5(tr("Which LR decision threshold should I use?")),
                actionLink("go_choose", tr("Non-genetic: 3. Choose threshold")),
                actionLink("go_gen_optimal", tr("Genetic: optimal decision threshold"))
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
            p("Marsico F.L., Vigeland M.D., et al. (2021). Forensic Science International:
              Genetics 54:102519")
          )
        )
      )
    ),

    tabPanel(
      tr("Non-Genetic Evidence"), value = "Non-Genetic Evidence",
      br(),
      fluidRow(column(10, offset = 1,
        p(class = "interpretation-text",
          tr("Calculate and simulate likelihood ratios from non-genetic evidence types (sex, age, hair color). All analyses in this section can be run directly in the web application."))
      )),
      tabsetPanel(
        id = "nongen_tabs",
        type = "pills",

    tabPanel(
      tr("Individual Evidence"), value = "Individual Evidence",
      fluidRow(
        column(12,
          br(),
          h3(tr("Single Evidence Type Analysis")),
          p(tr("Calculate the likelihood ratio for one type of evidence about the missing person (MP)."))
        )
      ),
      fluidRow(
        column(3,
          wellPanel(
            selectInput("indiv_type", tr("Evidence Type:"),
              choices = c("Biological Sex" = "sex",
                          "Age" = "age",
                          "Hair Color" = "hair"), selectize = FALSE),
            hr(),

            conditionalPanel(
              condition = "input.indiv_type == 'sex'",
              h5(tr("Case Parameters")),
              selectInput("indiv_MPs", tr("MP reported sex:"), sex_choices, selectize = FALSE),
              selectInput("indiv_obsSex", tr("Observed sex:"), sex_choices, selectize = FALSE),
              hr(),
              p(class = "text-muted small", tr(defaults_help)),
              h5(tr("Model Parameters")),
              sliderInput("indiv_eps_sex", tr("Observation error rate:"),
                0, 0.3, 0.05, 0.01),
              sliderInput("indiv_propF", tr("Female proportion in population:"),
                0.3, 0.7, 0.51, 0.01)
            ),

            conditionalPanel(
              condition = "input.indiv_type == 'age'",
              h5(tr("Case Parameters")),
              numericInput("indiv_MPa", tr("MP estimated age (years):"), 35, 1, 80),
              numericInput("indiv_MPr", tr("Age uncertainty (\u00b1 years):"), 5, 1, 20),
              numericInput("indiv_obsAge", tr("Observed age (years):"), 37, 1, 80),
              hr(),
              p(class = "text-muted small", tr(defaults_help)),
              h5(tr("Model Parameters")),
              sliderInput("indiv_epa", tr("Age category error rate:"),
                0, 0.3, 0.05, 0.01)
            ),

            conditionalPanel(
              condition = "input.indiv_type == 'hair'",
              h5(tr("Case Parameters")),
              selectInput("indiv_MPc", tr("MP hair color:"), hair_choices, selected = 2, selectize = FALSE),
              selectInput("indiv_obsC", tr("Observed hair color:"), hair_choices, selected = 2, selectize = FALSE),
              hr(),
              p(class = "text-muted small", tr(defaults_help)),
              h5(tr("Model Parameters")),
              sliderInput("indiv_epc", tr("Color misclassification rate (total):"),
                0, 0.3, 0.05, 0.01),
              p(class = "text-muted small", tr(hair_error_help)),
              p(class = "text-muted small", tr(hair_freq_help))
            )
          )
        ),
        column(9,
          div(class = "result-box",
            h4(tr("Likelihood Ratio")),
            h2(textOutput("indiv_lr_value")),
            p(textOutput("indiv_lr_log")),
            p(class = "plain-sentence", textOutput("indiv_interp", inline = TRUE))
          ),
          div(class = "formula-box",
            h5(style = "margin-top: 0;", tr("Calculation Details")),
            verbatimTextOutput("indiv_formula")
          )
        )
      )
    ),

    tabPanel(
      tr("Combine Evidence"), value = "Combine Evidence",
      fluidRow(
        column(12,
          br(),
          h3(tr("Evidence Combination")),
          p(tr("Under conditional independence, individual LRs for the missing person (MP) multiply to give a combined LR."))
        )
      ),
      fluidRow(
        column(12,
          div(class = "warning-box",
            strong(tr("Independence Assumption: ")),
            tr("LR multiplication assumes conditional independence of evidence sources.")
          )
        )
      ),
      fluidRow(
        column(3,
          wellPanel(
            h5(tr("Sex Evidence")),
            checkboxInput("comb_use_sex", tr("Include"), TRUE),
            conditionalPanel(
              condition = "input.comb_use_sex",
              selectInput("comb_MPs", tr("MP:"), sex_choices, selectize = FALSE),
              selectInput("comb_obsSex", tr("Observed:"), sex_choices, selectize = FALSE),
              numericInput("comb_eps", tr("Observation error rate (0\u20131):"), 0.05, 0, 0.3, 0.01),
              p(class = "text-muted small", tr(error_help)),
              sliderInput("comb_propF", tr("Female proportion in population:"),
                0.3, 0.7, 0.51, 0.01)
            ),
            hr(),
            strong(textOutput("comb_lr_sex"))
          )
        ),
        column(3,
          wellPanel(
            h5(tr("Age Evidence")),
            checkboxInput("comb_use_age", tr("Include"), TRUE),
            conditionalPanel(
              condition = "input.comb_use_age",
              numericInput("comb_MPa", tr("MP age (years):"), 35, 1, 80),
              numericInput("comb_MPr", tr("Age uncertainty (\u00b1 years):"), 5, 1, 20),
              numericInput("comb_obsAge", tr("Observed age (years):"), 37, 1, 80),
              numericInput("comb_epa", tr("Observation error rate (0\u20131):"), 0.05, 0, 0.3, 0.01),
              p(class = "text-muted small", tr(error_help))
            ),
            hr(),
            strong(textOutput("comb_lr_age"))
          )
        ),
        column(3,
          wellPanel(
            h5(tr("Hair Color Evidence")),
            checkboxInput("comb_use_hair", tr("Include"), TRUE),
            conditionalPanel(
              condition = "input.comb_use_hair",
              selectInput("comb_MPc", tr("MP:"), hair_choices, selected = 2, selectize = FALSE),
              selectInput("comb_obsC", tr("Observed:"), hair_choices, selected = 2, selectize = FALSE),
              numericInput("comb_epc", tr("Observation error rate, total (0\u20131):"), 0.05, 0, 0.3, 0.01),
              p(class = "text-muted small", tr(hair_error_help))
            ),
            hr(),
            strong(textOutput("comb_lr_hair"))
          )
        ),
        column(3,
          div(class = "result-box",
            h4(tr("Combined Result")),
            h2(textOutput("comb_lr_total")),
            p(textOutput("comb_lr_log")),
            p(class = "plain-sentence", style = "text-align: left;", textOutput("comb_interp", inline = TRUE)),
            hr(),
            verbatimTextOutput("comb_formula")
          )
        )
      )
    ),

    tabPanel(
      tr("1. Define profile"), value = "Probability Tables",
      fluidRow(
        column(12,
          br(),
          h3(tr("Conditional Probability Tables")),
          p(tr("Define the missing person (MP) profile and error rates. 2. LR distribution and 3. Choose threshold use these settings."))
        )
      ),
      fluidRow(
        column(2,
          wellPanel(
            h5(tr("Missing person (MP)")),
            selectInput("cpt_MPs", tr("Sex:"), sex_choices, selectize = FALSE),
            numericInput("cpt_MPa", tr("Age (years):"), 40, 1, 80),
            numericInput("cpt_MPr", tr("Age uncertainty (\u00b1 years):"), 6, 1, 20),
            selectInput("cpt_MPc", tr("Hair color:"), hair_choices, selected = 2, selectize = FALSE)
          )
        ),
        column(2,
          wellPanel(
            h5(tr("Observation error rates (0\u20131)")),
            p(class = "text-muted small", tr(error_help)),
            numericInput("cpt_eps", tr("Sex:"), 0.05, 0, 0.5, 0.01),
            numericInput("cpt_epa", tr("Age:"), 0.05, 0, 0.5, 0.01),
            numericInput("cpt_epc", tr("Hair (total):"), 0.05, 0, 0.5, 0.01),
            p(class = "text-muted small", tr(hair_error_help))
          ),
          wellPanel(
            h5(tr("Population")),
            sliderInput("cpt_propF", tr("Female proportion (0\u20131):"), 0.3, 0.7, 0.51, 0.01)
          )
        ),
        column(8,
          plotOutput("cpt_plots", height = "auto"),
          div(class = "info-box",
            fluidRow(
              column(4,
                strong("P(D|H2):"), " ", tr("Phenotype probability in reference population")
              ),
              column(4,
                strong("P(D|H1):"), " ", tr("Probability under H1: the person is the MP (related)")
              ),
              column(4,
                strong("log10(LR):"), " ",
                tr("Blue = supports H1, orange = supports H2, white cell = LR = 0")
              )
            ),
            p(class = "text-muted small", style = "margin: 8px 0 0 0;",
              tr("Rows: F/M = sex; T1 = age inside the missing person's range, T0 = outside."))
          )
        )
      )
    ),

    tabPanel(
      tr("2. LR distribution"), value = "LR Distribution",
      fluidRow(
        column(12,
          br(),
          h3(tr("Likelihood Ratio Distributions")),
          p(tr("Comparison of LR distributions under H1 and H2.")),
          uiOutput("profile_banner_dist")
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
          div(class = "success-box h1-box",
            h5(tr("H1: the person is the MP (related)")),
            verbatimTextOutput("dist_stats_h1")
          )
        ),
        column(6,
          div(class = "warning-box h2-box",
            h5(tr("H2: unrelated person")),
            verbatimTextOutput("dist_stats_h2")
          )
        )
      )
    ),

    tabPanel(
      tr("3. Choose threshold"), value = "Decision Analysis",
      fluidRow(
        column(12,
          br(),
          h3(tr("Decision Threshold Selection")),
          p(tr("Evaluate the trade-off between false positive and false negative rates.")),
          uiOutput("profile_banner_dec")
        )
      ),
      fluidRow(
        column(4,
          wellPanel(
            h5(tr("Decision threshold")),
            sliderInput("dec_threshold", tr("log10(LR) decision threshold:"),
                               min = -3, max = 3, value = 0, step = 0.1),
            p(class = "text-muted small", textOutput("dec_threshold_lr", inline = TRUE)),
            hr(),
            h5(tr("Performance Metrics")),
            verbatimTextOutput("dec_metrics"),
            tags$details(
              tags$summary(tr("Advanced metrics")),
              p(class = "small", strong("MCC: "), textOutput("dec_mcc", inline = TRUE)),
              p(class = "text-muted small", tr(mcc_help))
            )
          ),
          div(class = "info-box",
            h6(tr("Metric Definitions")),
            tags$ul(class = "small",
              tags$li(strong("TPR:"), " ", tr("True Positive Rate (sensitivity)")),
              tags$li(strong("FPR:"), " ", tr("False Positive Rate (1 - specificity)")),
              tags$li(strong(tr("Rule:")), " ", tr("a person is declared a match when log10(LR) >= threshold"))
            )
          )
        ),
        column(8,
          plotOutput("dec_roc_plot", height = "400px")
        )
      )
    ),

    tabPanel(
      tr("Database Search"), value = "Database Search",
      fluidRow(
        column(12,
          br(),
          h3(tr("Database Search")),
          p(tr("Compute the non-genetic LR of every candidate in a database against the missing person (MP) profile, and count candidates at or above a decision threshold (expected false positives in a database search, Marsico & Caridi 2023).")),
          div(class = "profile-banner",
            tr("For genetic false-positive rates, use Genetic Evidence > Rates at Threshold."), " ",
            actionLink("go_gen_rates_dbs", tr("Go there.")))
        )
      ),
      fluidRow(
        column(3,
          wellPanel(
            h5(tr("Candidate Database")),
            fileInput("dbs_file", NULL, accept = c(".csv", ".tsv", ".txt", "text/csv", "text/plain"),
              buttonLabel = tr("Browse...")),
            p(class = "text-muted small",
              tr("CSV with columns ID, Sex (F/M), Age (whole number 1-80) and HairColor (Black, Brown, Blonde, Red, Gray). Comma, semicolon or tab separated.")),
            actionButton("dbs_load_example", tr("Load example"), class = "btn-primary btn-block"),
            downloadButton("dbs_download_example", tr("Download example"),
              class = "btn-default btn-block", style = "margin-top: 6px;"),
            p(class = "text-muted small", style = "margin-top: 6px;",
              tr("Synthetic example: 2,001 candidates. C0001 matches the default MP profile, but every candidate in the same sex, age-range and hair-colour class gets the same LR: non-genetic evidence ranks classes, it does not identify individuals.")),
            uiOutput("dbs_msg"),
            hr(),
            h5(tr("Missing person (MP)")),
            selectInput("dbs_MPs", tr("Sex:"), sex_choices, selectize = FALSE),
            numericInput("dbs_MPa", tr("Age (years):"), 40, 1, 80),
            numericInput("dbs_MPr", tr("Age uncertainty (\u00b1 years):"), 6, 1, 20),
            selectInput("dbs_MPc", tr("Hair color:"), hair_choices, selected = 2, selectize = FALSE),
            hr(),
            h5(tr("Observation error rates (0\u20131)")),
            p(class = "text-muted small", tr(error_help)),
            numericInput("dbs_eps", tr("Sex:"), 0.05, 0, 0.5, 0.01),
            numericInput("dbs_epa", tr("Age:"), 0.05, 0, 0.5, 0.01),
            numericInput("dbs_epc", tr("Hair (total):"), 0.05, 0, 0.5, 0.01),
            p(class = "text-muted small", tr(hair_error_help)),
            hr(),
            checkboxInput("dbs_empirical", tr("Estimate population frequencies from the database"), TRUE),
            p(class = "text-muted small",
              tr("Unchecked: female proportion 0.51, default hair color frequencies and the uniform 1-80 age model."))
          )
        ),
        column(9,
          wellPanel(
            sliderInput("dbs_threshold", tr("log10(LR) decision threshold:"),
              min = -3, max = 3, value = 0, step = 0.1),
            p(class = "text-muted small", textOutput("dbs_threshold_lr", inline = TRUE)),
            p(class = "plain-sentence", textOutput("dbs_sentence", inline = TRUE)),
            br(),
            verbatimTextOutput("dbs_summary")
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
      fluidRow(column(12,
        tags$details(class = "info-box",
          tags$summary(tr("About Genetic Simulations")),
          p(tr("The web version of mispitools cannot perform genetic simulations due to computational requirements. Simulations must be run using dedicated software and then uploaded here for analysis:")),
          tags$ul(
            tags$li(strong("Familias"), " - ", tr("standalone software for forensic genetics")),
            tags$li(strong("forrel"), " - ", tr("R package for forensic pedigree analysis")),
            tags$li(strong("mispitools"), " - ", tr("R package (command-line version)")),
            tags$li(strong("FamLink"), " - ", tr("software for kinship analysis"))
          ),
          p(tr("Once you have simulation results, upload them below to analyze statistical power, LR distributions, and optimize decision thresholds."))
        )
      )),
      fluidRow(
        column(12,
          h3(tr("Simulation Results Analysis")),
          p(class = "interpretation-text",
            tr("Upload simulation data from forensic genetics software (Familias, Forrel, FamLink, mispitools) to evaluate statistical power and optimize LR decision thresholds.")
          )
        )
      ),
      fluidRow(
        column(4,
          wellPanel(
            div(class = "upload-zone",
              h5(tr("Upload Simulation Data")),
              fileInput("sim_file", NULL,
                accept = c(".csv", ".tsv", ".txt",
                           "text/csv", "text/tab-separated-values", "text/plain"),
                buttonLabel = tr("Browse...")),
              p(class = "text-muted small",
                tr("Accepted: Familias, FamLink, forrel, mispitools, or any CSV/TSV with Unrelated & Related columns.")),
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
              h5(tr("LR decision threshold")),
              numericInput("sim_threshold", tr("LR threshold value (match if LR >= threshold):"),
                value = 100, min = 0, max = 1000000, step = 1),
              numericInput("sim_db_size", tr("N = number of unrelated candidates compared:"),
                value = 10000, min = 1, step = 100),
              hr(),
              h5(tr("Optimal decision threshold")),
              sliderInput("sim_weight", tr("How much worse is a false match than a missed relative? (\u00d7)"),
                min = 1, max = 50, value = 10, step = 1),
              p(class = "text-muted small",
                tr("Weight > 1 penalizes false positives more. Default: 10 (FP is 10x worse than FN).")),
              actionButton("sim_calc_dt", tr("Compute optimal threshold"),
                class = "btn-primary btn-block"),
              uiOutput("sim_dt_sentence"),
              br(),
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
                        span(class = "metric-label", tr("Proportion of relatives with LR \u2265 10,000 (fixed reference)")),
                        div(class = "metric-value", textOutput("sim_ip"))
                      )
                    )
                  ),
                  div(style = "overflow-x: auto;", tableOutput("sim_summary_table")),
                  tags$details(
                    tags$summary(tr("Advanced metrics")),
                    p(class = "small", strong(tr("Exclusion power: ")), textOutput("sim_ep", inline = TRUE)),
                    p(class = "text-muted small", tr("Proportion of unrelated simulations with LR = 0, i.e. excluded with certainty."))
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
                  div(class = "info-box",
                    p(tr("Each point represents a different LR threshold. The curve shows the trade-off between false negative rate (missing a true match) and false positive rate (incorrectly flagging a non-match). Points closer to the origin (0,0) indicate better thresholds. Grey labels show the LR value at selected points."))
                  ),
                  plotOutput("sim_decision_plot", height = "450px")
                ),

                tabPanel(
                  tr("Rates at Threshold"), value = "Rates at Threshold",
                  br(),
                  div(class = "info-box",
                    p(class = "plain-sentence", style = "margin: 0;", textOutput("sim_rates_sentence", inline = TRUE))
                  ),
                  p(class = "text-muted small",
                    tr("A comparison is declared a match when LR >= threshold. Set the decision threshold and N in the left panel.")),
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
                  fluidRow(
                    column(4, offset = 4,
                      div(class = "metric-card",
                        span(class = "metric-label", tr("Expected false positives (FPR x N)")),
                        div(class = "metric-value", textOutput("sim_expected_fp"))
                      )
                    )
                  ),
                  tableOutput("sim_rates_table"),
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
            div(style = "padding: 10px 0; color: #555;",
              h4(tr("Upload simulation data to begin analysis"))
            ),
            div(class = "info-box",
              h5(tr("Accepted file formats:")),
              tags$ol(
                tags$li(tags$b(tr("Familias (single family):")), " ",
                  tr("Tab-separated .txt file with 2 header rows (case name + pedigree label), then numeric LR data. First column = Unrelated (Reference Pedigree), second = Related.")),
                tags$li(tags$b(tr("Familias (multi-family):")), " ",
                  tr("Same structure but with column pairs per family. Row 1 repeats the family name, row 2 has 'Reference pedigree' for H2 and a relationship label (Father, Brother, etc.) for H1. Can contain one or more families.")),
                tags$li(tags$b("FamLink:"), " ",
                  tr("Starts with 'Raw data LR:s...', then 'Pedigree (Numerator)', 'Versus (Denominator)' and 'True:' header rows. Groups are separated by '|'. The 'True' row indicates which column is Related (H1). Can contain one or more pedigree comparisons.")),
                tags$li(tags$b(tr("Generic CSV/TSV:")), " ",
                  tr("Any .csv, .tsv, or .txt file with columns named 'Unrelated' and 'Related', or two numeric columns (first = Unrelated, second = Related). Supports comma, tab, or semicolon separators."))
              ),
              hr(),
              h5(tr("Example format (Familias, tab-separated):")),
              tags$pre(style = "background: #f7f7f7; padding: 10px; font-size: 12px; overflow-x: auto;",
                "FamilyName\\tFamilyName\nReference pedigree\\tFather\n3.33333\\t5\n2.5\\t3.33333\n1.25\\t2.08333"
              ),
              h5(tr("Example format (Generic CSV):")),
              tags$pre(style = "background: #f7f7f7; padding: 10px; font-size: 12px; overflow-x: auto;",
                "Unrelated,Related\n0.5,3.2\n0.8,5.1\n0.3,2.7"
              ),
              p(class = "text-muted small",
                tr("Note: Files must be text-based (.txt, .csv, .tsv). Excel files (.xlsx) are not supported. Simulations can be generated using Familias, forrel, mispitools (R version), or FamLink."))
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
    indiv_type = c("Biological Sex" = "sex", "Age" = "age", "Hair Color" = "hair"),
    indiv_MPs = sex_choices, indiv_obsSex = sex_choices, cpt_MPs = sex_choices,
    comb_MPs = sex_choices, comb_obsSex = sex_choices, dbs_MPs = sex_choices,
    indiv_MPc = hair_choices, indiv_obsC = hair_choices, cpt_MPc = hair_choices,
    comb_MPc = hair_choices, comb_obsC = hair_choices, dbs_MPc = hair_choices,
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

  observeEvent(input$go_individual, go_to("Non-Genetic Evidence", "Individual Evidence"))
  observeEvent(input$go_combine, go_to("Non-Genetic Evidence", "Combine Evidence"))
  observeEvent(input$go_dbs, go_to("Non-Genetic Evidence", "Database Search"))
  observeEvent(input$go_choose, go_to("Non-Genetic Evidence", "Decision Analysis"))
  observeEvent(input$go_gen_rates, go_to("Genetic Evidence", sim = "Rates at Threshold"))
  observeEvent(input$go_gen_rates_dbs, go_to("Genetic Evidence", sim = "Rates at Threshold"))
  observeEvent(input$go_gen_optimal, go_to("Genetic Evidence", sim = "Decision Plot"))
  observeEvent(input$goto_profile_dist, go_to("Non-Genetic Evidence", "Probability Tables"))
  observeEvent(input$goto_profile_dec, go_to("Non-Genetic Evidence", "Probability Tables"))

  profile_banner <- function(link_id) {
    div(class = "profile-banner",
      trf("These results use the missing-person profile and error rates set in 1. Define profile (Sex: %s, Age %s \u00b1 %s, Hair: %s).",
          tr(names(sex_choices)[sex_choices == input$cpt_MPs]),
          input$cpt_MPa, input$cpt_MPr,
          tr(names(hair_choices)[hair_choices == input$cpt_MPc])),
      " ",
      actionLink(link_id, tr("Change them there.")))
  }

  output$profile_banner_dist <- renderUI(profile_banner("goto_profile_dist"))
  output$profile_banner_dec <- renderUI(profile_banner("goto_profile_dec"))

  interpretation <- function(LR, evidence) {
    if (LR == 0) return(tr("The observation is impossible under the stated model (exclusion)."))
    if (LR == 1) return(tr("The observation is equally probable under both hypotheses."))
    paste(trf(interpretation_templates[[evidence]][if (LR > 1) 1 else 2],
              sprintf("%.4g", if (LR > 1) LR else 1 / LR)),
          tr("This is not the probability that they are the missing person."))
  }

  output$indiv_interp <- renderText({
    interpretation(indiv_lr(), input$indiv_type)
  })

  output$comb_interp <- renderText({
    interpretation(comb_total(), "combined")
  })

  threshold_lr_text <- function(t) {
    sprintf("log10(LR) \u2265 %.1f (LR \u2265 %s)", t, format(signif(10^t, 4), big.mark = big_mark(), scientific = FALSE))
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
        check_num(input$indiv_MPr, "Age uncertainty", 1, 20, whole = TRUE),
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
      concordant <- input$indiv_MPs == input$indiv_obsSex
      cat(trf("Sex Evidence (%s)", tr(if (concordant) "concordant" else "discordant")), "\n", sep = "")
      cat("----------------------------------------\n")
      if (concordant) {
        cat(tr("P(observed | H1) = 1 - eps ="), 1 - input$indiv_eps_sex, "\n")
      } else {
        cat(tr("P(observed | H1) = eps ="), input$indiv_eps_sex, "\n")
      }
      obs_freq <- if (input$indiv_obsSex == "F") input$indiv_propF else (1 - input$indiv_propF)
      cat(tr("P(observed | H2) ="), round(obs_freq, 3), "\n")

    } else if (input$indiv_type == "age") {
      aw <- age_window(input$indiv_MPa, input$indiv_MPr)
      in_range <- input$indiv_obsAge >= aw$min && input$indiv_obsAge <= aw$max
      cat(trf("Age Evidence (%s)", tr(if (in_range) "in range" else "out of range")), "\n", sep = "")
      cat("----------------------------------------\n")
      cat(tr("MP age:"), input$indiv_MPa, "+/-", input$indiv_MPr, tr("years"), "\n")
      cat(tr("Range used (ages 1-80):"), aw$min, tr("to"), aw$max,
          "=", aw$max - aw$min + 1, tr("integer ages"), "\n")
      cat(tr("Observed:"), input$indiv_obsAge, tr("years"), "\n")
      if (in_range) {
        cat(tr("P(in range | H1) = 1 - eps ="), 1 - input$indiv_epa, "\n")
        cat(tr("P(in range | H2) ="), aw$max - aw$min + 1, "/ 80 =", round(aw$T1p, 4), "\n")
      } else {
        cat(tr("P(out of range | H1) = eps ="), input$indiv_epa, "\n")
        cat(tr("P(out of range | H2) = 1 -"), aw$max - aw$min + 1, "/ 80 =",
            round(1 - aw$T1p, 4), "\n")
      }

    } else {
      concordant <- input$indiv_MPc == input$indiv_obsC
      cat(trf("Hair Color Evidence (%s)", tr(if (concordant) "concordant" else "discordant")), "\n", sep = "")
      cat("----------------------------------------\n")
      cat(tr("MP:"), tr(hair_labels[[input$indiv_MPc]]), "\n")
      cat(tr("Observed:"), tr(hair_labels[[input$indiv_obsC]]), "\n")
      if (concordant) {
        cat(tr("P(observed | H1) = 1 - eps ="), 1 - input$indiv_epc, "\n")
      } else {
        cat(tr("P(observed | H1) = eps / 4 ="), input$indiv_epc / 4, "\n")
      }
      cat(tr("P(observed | H2) = population frequency ="),
          DEFAULT_HAIR_FREQ[as.numeric(input$indiv_obsC)], "\n")
    }
    cat(if (LR == 0) tr("LR = 0 (exclusion)") else paste("LR =", round(LR, 4)), "\n")
  })

  cpt_pop <- reactive({
    validate(
      check_num(input$cpt_MPa, "Age", 1, 80, whole = TRUE),
      check_num(input$cpt_MPr, "Age uncertainty", 1, 20, whole = TRUE),
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
  }, res = 100, height = function() if (is_narrow("cpt_plots")) 1000 else 380,
  alt = function() tr("Heatmaps of P(D | H2), P(D | H1) and log10(LR) for each sex, age class and hair colour of the missing person profile."))

  lr_dists <- reactive({
    build_lr_distributions(cpt_pop(), cpt_mp())
  })

  output$dist_plot <- renderPlot({
    dists <- lr_dists()
    create_dist_plot(dists$H1, dists$H2, base_size = plot_size("dist_plot", 12))
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

  output$dist_roc <- renderPlot({
    dists <- lr_dists()
    roc_obj <- build_roc(dists$H1, dists$H2)
    df_roc <- roc_obj$data
    auc_val <- roc_obj$AUC

    if (nrow(df_roc) == 0) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5, label = tr("Insufficient data")) +
               ggplot2::theme_void())
    }

    ggplot2::ggplot(df_roc, ggplot2::aes(x = FPR, y = TPR)) +
      ggplot2::geom_path(color = "grey35", linewidth = 1) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                           color = COL_REFERENCE_LINE, linewidth = 0.8) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::labs(title = roc_title(auc_val),
                    x = tr("False Positive Rate (1 - Specificity)"),
                    y = tr("True Positive Rate\n(Sensitivity)")) +
      ggplot2::theme_minimal(base_size = plot_size("dist_roc", 12)) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }, res = 100, alt = roc_alt)

  comb_lr_sex_val <- reactive({
    if (input$comb_use_sex) {
      validate(
        check_num(input$comb_eps, "Sex observation error rate", 0, 0.3),
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
        check_num(input$comb_MPr, "Age uncertainty", 1, 20, whole = TRUE),
        check_num(input$comb_obsAge, "Observed age", 1, 80, whole = TRUE),
        check_num(input$comb_epa, "Age observation error rate", 0, 0.3)
      )
      calc_lr_age(input$comb_MPa, input$comb_MPr, input$comb_obsAge, input$comb_epa)
    } else {
      1
    }
  })

  comb_lr_hair_val <- reactive({
    if (input$comb_use_hair) {
      validate(check_num(input$comb_epc, "Hair observation error rate", 0, 0.3))
      calc_lr_hair(as.numeric(input$comb_MPc), as.numeric(input$comb_obsC),
                   input$comb_epc, DEFAULT_HAIR_FREQ)
    } else {
      1
    }
  })

  output$comb_lr_sex <- renderText({
    if (input$comb_use_sex) format_lr(comb_lr_sex_val()) else tr("Not included")
  })

  output$comb_lr_age <- renderText({
    if (input$comb_use_age) format_lr(comb_lr_age_val()) else tr("Not included")
  })

  output$comb_lr_hair <- renderText({
    if (input$comb_use_hair) format_lr(comb_lr_hair_val()) else tr("Not included")
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
      cat(tr("No evidence selected"), "\n")
    } else {
      cat(paste(sprintf("%.4g", vals), collapse = " x "), "\n")
      cat("=", sprintf("%.4g", prod(vals)), "\n")
    }
  })

  output$dec_metrics <- renderPrint({
    dists <- lr_dists()
    m <- calc_metrics(dists$H1, dists$H2, input$dec_threshold)

    rows <- c(tr("True Positive Rate:"), tr("False Positive Rate:"),
              tr("False Negative Rate:"), tr("True Negative Rate:"))
    width <- max(nchar(rows))
    cat(tr("At decision threshold log10(LR) ="), input$dec_threshold, "\n")
    cat(tr("Match declared when log10(LR) >= threshold"), "\n")
    cat("========================================\n\n")
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
                    subtitle = trf("Black point: selected decision threshold (log10(LR) >= %.1f)",
                                   input$dec_threshold),
                    x = tr("False Positive Rate"), y = tr("True Positive Rate")) +
      ggplot2::theme_minimal(base_size = plot_size("dec_roc_plot", 14)) +
      ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }, res = 100, alt = function() paste(roc_alt(),
    trf("Black point: selected decision threshold (log10(LR) >= %.1f)", input$dec_threshold)))

  sim_families_val <- reactiveVal(NULL)
  sim_parse_error  <- reactiveVal(NULL)

  observeEvent(input$sim_file, {
    result <- tryCatch(parse_sim_upload(input$sim_file$datapath), error = function(e) NULL)
    if (is.null(result)) {
      sim_families_val(NULL)
      sim_parse_error(c(paste0(
        "Could not parse the uploaded file ('%s'). ",
        "Please check the accepted formats below.\n\n",
        "Common issues:\n",
        "- The file must be a text file (.txt, .csv, .tsv), not Excel (.xlsx)\n",
        "- For Familias format: Row 1 = family names, Row 2 = pedigree labels ",
        "(e.g., 'Reference pedigree', 'Father'), then numeric LR values\n",
        "- For FamLink format: must start with 'Raw data LR' and have ",
        "'Pedigree (Numerator)', 'Versus (Denominator)' and 'True:' header rows\n",
        "- For generic CSV/TSV: need at least 2 numeric columns\n",
        "- Check that the separator is tab, comma, or semicolon"
      ), input$sim_file$name))
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
                    subtitle = tr("Vertical grey line: LR = 1 (neutral). Less overlap = better discrimination."),
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
          tr("Green diamond: optimal decision threshold\nBlack dot: selected decision threshold\n(match if LR >= threshold)")
        } else {
          tr("Green diamond: optimal decision threshold | Black dot: selected decision threshold\n(match if LR >= threshold)")
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
    out <- trf("With LR \u2265 %s, searching %s unrelated people you would expect about %s false matches, and you would detect the true relative %s of the time.",
               format(input$sim_threshold, big.mark = big_mark(), scientific = FALSE),
               format(sim_n_candidates(), big.mark = big_mark(), scientific = FALSE),
               sprintf("%.1f", r$FPR * sim_n_candidates()),
               sprintf("%.1f%%", 100 * r$TPR))
    if (!is.na(sim_ppv())) {
      out <- paste(out, trf("If one true relative is among them, a reported match is the true one %s of the time.",
                            sprintf("%.1f%%", 100 * sim_ppv())))
    }
    out
  })

  output$sim_advanced <- renderUI({
    r <- sim_rates()
    validate(sim_n_message())
    item <- function(label, value, help) {
      tagList(p(class = "small", style = "margin: 6px 0 0 0;", strong(label), " ", value),
              p(class = "text-muted small", style = "margin: 0;", help))
    }
    tagList(
      item(tr("PPV (1 true match among N candidates)"),
           if (is.na(sim_ppv())) tr("undefined") else sprintf("%.4f", sim_ppv()),
           tr("Probability that a reported match is the true relative, assuming one true relative among N unrelated candidates.")),
      item(tr("Accuracy (balanced classes)"), sprintf("%.4f", r$Accuracy),
           tr("Proportion of correct decisions when related and unrelated comparisons are equally frequent.")),
      item(tr("Matthews Correlation Coefficient (MCC, balanced classes)"), sprintf("%.4f", r$MCC),
           tr(mcc_help))
    )
  })

  output$sim_expected_fp <- renderText({
    validate(sim_n_message())
    sprintf("%.1f", sim_rates()$FPR * sim_n_candidates())
  })

  output$sim_rates_table <- renderTable({
    r <- sim_rates()
    out <- data.frame(
      Metric = c(tr("True Positive Rate (TPR)"),
                 tr("False Positive Rate (FPR)"),
                 tr("True Negative Rate (TNR)"),
                 tr("False Negative Rate (FNR)"),
                 trf("Expected false positives (FPR x N, N = %s)",
                     format(sim_n_candidates(), big.mark = big_mark()))),
      Value = c(sprintf("%.6f", r$TPR),
                sprintf("%.6f", r$FPR),
                sprintf("%.6f", r$TNR),
                sprintf("%.6f", r$FNR),
                sprintf("%.1f", r$FPR * sim_n_candidates()))
    )
    names(out) <- c(tr("Metric"), tr("Value"))
    out
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

  output$sim_dt_sentence <- renderUI({
    req(sim_dt_shown())
    result <- sim_opt()
    validate(sim_n_message())
    div(class = "info-box", style = "margin-bottom: 0;",
      p(style = "margin: 0 0 8px 0;",
        trf("With weight w = %s, the best threshold is LR \u2265 %.2f: about %s false matches per %s people compared, and %s of true relatives detected.",
            input$sim_weight, result$DT, sprintf("%.1f", result$FPR * sim_n_candidates()),
            format(sim_n_candidates(), big.mark = big_mark(), scientific = FALSE),
            sprintf("%.1f%%", 100 * (1 - result$FNR)))),
      actionButton("sim_use_dt", tr("Use this threshold"), class = "btn-default btn-block"))
  })

  observeEvent(input$sim_use_dt, {
    updateNumericInput(session, "sim_threshold", value = sim_opt()$DT)
  })

  output$sim_dt_result <- renderPrint({
    if (!sim_dt_shown()) {
      cat(tr("Click 'Compute optimal threshold' to calculate."), "\n")
    } else {
      result <- sim_opt()
      rows <- c(tr("LR threshold:"), tr("FPR at optimal threshold:"), tr("FNR at optimal threshold:"), tr("Weight used:"))
      rows <- formatC(rows, width = -max(nchar(rows)))
      cat(tr("Optimal decision threshold"), "\n")
      cat("========================================\n\n")
      cat("  ", rows[1], " ", sprintf("%.2f", result$DT), "\n", sep = "")
      cat("  ", rows[2], " ", sprintf("%.6f", result$FPR), "\n", sep = "")
      cat("  ", rows[3], " ", sprintf("%.6f", result$FNR), "\n", sep = "")
      cat("  ", rows[4], " ", input$sim_weight, "\n", sep = "")
      cat("\n  ", tr("Formula:"), " D = sqrt(FNR^2 + (w*FPR)^2)\n", sep = "")
      cat("  ", tr("Match declared when LR >= threshold"), "\n", sep = "")
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
      div(class = "warning-box", style = "padding: 10px;", tags$b(tr("File error: ")), fmt_msg(err))
    } else if (!is.null(dbs_data())) {
      div(class = "success-box", style = "padding: 10px;",
          trf("Loaded %s candidates.", format(nrow(dbs_data()), big.mark = big_mark())))
    }
  })

  dbs_lr <- reactive({
    cands <- dbs_data()
    validate(need(!is.null(cands), tr("Upload a candidate database or load the example.")))
    validate(
      check_num(input$dbs_MPa, "Age", 1, 80, whole = TRUE),
      check_num(input$dbs_MPr, "Age uncertainty", 1, 20, whole = TRUE),
      check_num(input$dbs_eps, "Sex error rate", 0, 0.5),
      check_num(input$dbs_epa, "Age error rate", 0, 0.5),
      check_num(input$dbs_epc, "Hair error rate", 0, 0.5)
    )
    candidate_lrs(cands, input$dbs_MPs, input$dbs_MPa, input$dbs_MPr, as.numeric(input$dbs_MPc),
                  input$dbs_eps, input$dbs_epa, input$dbs_epc, isTRUE(input$dbs_empirical))
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
          tr("All candidates in the same sex, age and hair class get the same LR, so non-genetic evidence narrows the list but cannot single out one person."))
  })

  output$dbs_summary <- renderPrint({
    out <- dbs_lr()
    res <- out$results
    n <- nrow(res)
    above <- sum(res$log10LR >= input$dbs_threshold)
    cat(tr("Candidates:"), n, "\n")
    cat(trf("LR >= %s (log10(LR) >= %.1f): %d (%.4f of candidates)",
            format(signif(10^input$dbs_threshold, 4)), input$dbs_threshold, above, above / n), "\n")
    cat(tr("LR = 0 (exclusions):"), sum(res$LR == 0), "\n\n")
    cat(tr(if (isTRUE(input$dbs_empirical)) "Frequencies estimated from the database:" else
          "Default population frequencies:"), "\n")
    cat(sprintf("  P(F) = %.4f\n", out$propF))
    cat(trf("  P(age %g-%g | H2) = %.4f", out$age_min, out$age_max, out$T1p), "\n")
    cat(" ", tr("Hair:"), paste(sprintf("%s %.4f", vapply(names(out$hair_freq), tr, character(1)),
                                       out$hair_freq), collapse = ", "), "\n")
    if (length(out$notes) > 0) {
      cat("\n", paste(tr("Note:"), out$notes, collapse = "\n"), "\n", sep = "")
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
      ggplot2::labs(title = dbs_above_title(),
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
