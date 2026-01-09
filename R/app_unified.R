#' Comprehensive Shiny App for Missing Person Identification
#'
#' @description
#' Launches a comprehensive interactive Shiny application for calculating
#' likelihood ratios (LRs) from non-genetic evidence in missing person cases.
#' This unified app integrates all evidence types (sex, age, hair color,
#' birthdate) with tutorials, visualizations, and decision analysis tools.
#'
#' @return A Shiny app object. When run interactively, launches a multi-tab
#'   web interface with:
#'   \itemize{
#'     \item \strong{Welcome}: Introduction to LR concepts
#'     \item \strong{Individual Evidence}: Calculate LR for each evidence type
#'     \item \strong{CPT Analysis}: Visualize conditional probability tables
#'     \item \strong{Distribution}: Simulate and visualize LR distributions
#'     \item \strong{Combine Evidence}: Combine multiple evidence types
#'     \item \strong{Decision Analysis}: Threshold selection and error metrics
#'     \item \strong{Tutorial}: Step-by-step educational content
#'   }
#'
#' @details
#' This app provides a complete workflow for forensic identification using
#' non-genetic evidence. It implements the Bayesian framework where:
#' \itemize{
#'   \item H1: The unidentified person IS the missing person
#'   \item H2: The unidentified person is NOT the missing person
#'   \item LR = P(Evidence | H1) / P(Evidence | H2)
#' }
#'
#' Evidence types supported:
#' \itemize{
#'   \item Biological sex (male/female)
#'   \item Age (within expected range)
#'   \item Hair color (5 categories)
#'   \item Birth date (discrepancy analysis)
#' }
#'
#' @seealso
#' \code{\link{lr_sex}}, \code{\link{lr_age}}, \code{\link{lr_hair_color}},
#' \code{\link{lr_birthdate}} for individual LR calculations,
#' \code{\link{lr_combine}} for combining evidence,
#' \code{\link{decision_threshold}}, \code{\link{threshold_rates}} for
#' decision analysis.
#'
#' @references
#' Marsico FL, Caridi I (2023). "Incorporating non-genetic evidence in large
#' scale missing person searches: A general approach beyond filtering."
#' \emph{Forensic Science International: Genetics}, 66, 102891.
#' \doi{10.1016/j.fsigen.2023.102891}
#'
#' Marsico FL, Vigeland MD, et al. (2021). "Making decisions in missing person
#' identification cases with low statistical power."
#' \emph{Forensic Science International: Genetics}, 52, 102519.
#' \doi{10.1016/j.fsigen.2021.102519}
#'
#' @import shiny
#' @import shinythemes
#' @import ggplot2
#' @import patchwork
#' @import pROC
#' @import reshape2
#' @import dplyr
#' @export
#' @examples
#' if (interactive()) {
#'   mispitools_app()
#' }
mispitools_app <- function() {

  # Suppress R CMD check NOTE for ggplot2 aesthetics
  Hypothesis <- logLR <- w <- Var1 <- Var2 <- value <- TPR <- FPR <- NULL

  # ============================================================================
  # CONSTANTS: Realistic forensic parameters based on published data

  # ============================================================================

  # Realistic population hair color frequencies (European population)
  # Source: Walsh et al. (2017) Forensic Sci Int Genet
  DEFAULT_HAIR_FREQ <- c(
    Black = 0.08,
    Brown = 0.45,
    Blonde = 0.30,
    Red = 0.07,
    Gray = 0.10
  )

  # Hair color labels
  hair_labels <- c("1" = "Black", "2" = "Brown", "3" = "Blonde",
                   "4" = "Red", "5" = "Gray/White")

  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================

  # Build CPT under H2 (Population)
  build_cpt_pop <- function(propF, MPa, MPr, propC) {
    Age <- seq(1, 80)
    MPmin <- max(1, MPa - MPr)
    MPmax <- min(80, MPa + MPr)
    T1p <- (MPmax - MPmin) / length(Age)
    T1p <- max(0.001, min(0.999, T1p))
    T0p <- 1 - T1p
    propS <- c(propF, 1 - propF)
    jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")
    jointprob <- c(propS[1]*T1p, propS[1]*T0p, propS[2]*T1p, propS[2]*T0p)
    names(jointprob) <- jointname
    CPT_POP <- outer(jointprob, propC)
    rownames(CPT_POP) <- jointname
    colnames(CPT_POP) <- names(propC)
    return(CPT_POP)
  }

  # Build CPT under H1 (MP)
  build_cpt_mp <- function(MPs, MPc, eps, epa, epc) {
    jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")
    if (MPs == "F") {
      jointprob <- c((1-eps)*(1-epa), (1-eps)*epa, eps*(1-epa), eps*epa)
    } else {
      jointprob <- c(eps*(1-epa), eps*epa, (1-eps)*(1-epa), (1-eps)*epa)
    }
    names(jointprob) <- jointname

    # Build error matrix for hair color (uniform error model)
    n_colors <- 5
    errorMat <- matrix(epc / (n_colors - 1), nrow = n_colors, ncol = n_colors)
    diag(errorMat) <- 1 - epc
    # Normalize rows
    errorMat <- errorMat / rowSums(errorMat)

    probC <- errorMat[MPc, ]
    CPT_MP <- outer(jointprob, probC)
    rownames(CPT_MP) <- jointname
    colnames(CPT_MP) <- names(DEFAULT_HAIR_FREQ)
    return(CPT_MP)
  }

  # Create 3-panel heatmap with improved aesthetics
  create_cpt_plots <- function(CPT_POP, CPT_MP) {
    dfPOP <- reshape2::melt(CPT_POP)
    dfMP <- reshape2::melt(CPT_MP)
    LR_mat <- CPT_MP / CPT_POP
    LR_mat[CPT_POP == 0] <- NA
    dfLR <- reshape2::melt(log10(LR_mat))

    base_theme <- ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        plot.subtitle = ggplot2::element_text(size = 9, color = "gray40"),
        axis.text = ggplot2::element_text(size = 9),
        legend.position = "bottom",
        legend.key.width = ggplot2::unit(1.5, "cm")
      )

    p1 <- ggplot2::ggplot(dfPOP, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c",
                                   limits = c(0, max(dfPOP$value) * 1.1)) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)),
                         size = 2.8, color = "black") +
      ggplot2::labs(title = "P(D | H2)", subtitle = "Reference population",
                    x = "Hair Color", y = "", fill = "Probability") +
      base_theme

    p2 <- ggplot2::ggplot(dfMP, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#08519c",
                                   limits = c(0, max(dfMP$value) * 1.1)) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)),
                         size = 2.8, color = "black") +
      ggplot2::labs(title = "P(D | H1)", subtitle = "If person IS the MP",
                    x = "Hair Color", y = "", fill = "Probability") +
      base_theme

    # Handle NA values in LR
    dfLR$value[is.na(dfLR$value)] <- 0
    dfLR$value[is.infinite(dfLR$value)] <- max(dfLR$value[is.finite(dfLR$value)], na.rm = TRUE)

    p3 <- ggplot2::ggplot(dfLR, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_gradient2(low = "#d73027", mid = "#ffffbf", high = "#1a9850",
                                    midpoint = 0, na.value = "gray80") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", value)),
                         size = 2.8, color = "black") +
      ggplot2::labs(title = expression(log[10](LR)), subtitle = "Evidence weight",
                    x = "Hair Color", y = "", fill = expression(log[10](LR))) +
      base_theme

    (p1 | p2 | p3) + patchwork::plot_layout(ncol = 3)
  }

  # Build LR distributions from CPTs
  build_lr_distributions <- function(CPT_POP, CPT_MP) {
    pop_vec <- as.vector(CPT_POP)
    mp_vec <- as.vector(CPT_MP)
    lr_vec <- mp_vec / pop_vec
    lr_vec[pop_vec == 0] <- NA
    log_lr <- log10(lr_vec)

    df_H1 <- data.frame(x = log_lr, w = mp_vec)
    df_H2 <- data.frame(x = log_lr, w = pop_vec)
    df_H1 <- df_H1[!is.na(df_H1$x) & is.finite(df_H1$x), ]
    df_H2 <- df_H2[!is.na(df_H2$x) & is.finite(df_H2$x), ]

    if (nrow(df_H1) > 0) df_H1$w <- df_H1$w / sum(df_H1$w)
    if (nrow(df_H2) > 0) df_H2$w <- df_H2$w / sum(df_H2$w)

    list(H1 = df_H1, H2 = df_H2)
  }

  # Create distribution histogram with improved aesthetics
  create_dist_plot <- function(distH1, distH2, binwidth = 0.25) {
    if (nrow(distH1) == 0 || nrow(distH2) == 0) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") +
               ggplot2::theme_void())
    }

    dH1 <- data.frame(logLR = distH1$x, w = distH1$w, Hypothesis = "H1: Person IS the MP")
    dH2 <- data.frame(logLR = distH2$x, w = distH2$w, Hypothesis = "H2: Person is NOT the MP")
    dAll <- rbind(dH1, dH2)

    ggplot2::ggplot(dAll, ggplot2::aes(x = logLR, weight = w, fill = Hypothesis)) +
      ggplot2::geom_histogram(binwidth = binwidth, position = "identity",
                              alpha = 0.65, color = "white", linewidth = 0.3) +
      ggplot2::scale_fill_manual(values = c("H1: Person IS the MP" = "#2171b5",
                                            "H2: Person is NOT the MP" = "#cb181d")) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
      ggplot2::labs(title = expression("Distribution of " * log[10](LR)),
                    subtitle = "Vertical line: LR = 1 (uninformative)",
                    x = expression(log[10](LR)), y = "Density") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "top",
        legend.title = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  # Calculate metrics at threshold
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

  # Build ROC curve
  build_roc <- function(distH1, distH2) {
    if (nrow(distH1) == 0 || nrow(distH2) == 0) {
      return(list(data = data.frame(th = numeric(0), TPR = numeric(0), FPR = numeric(0)),
                  AUC = NA))
    }

    all_vals <- sort(unique(c(distH1$x, distH2$x)))
    if (length(all_vals) < 2) {
      return(list(data = data.frame(th = 0, TPR = 1, FPR = 1), AUC = 0.5))
    }

    minx <- min(all_vals) - 0.1
    maxx <- max(all_vals) + 0.1
    thresholds <- seq(minx, maxx, length.out = 100)

    TPR_vec <- numeric(length(thresholds))
    FPR_vec <- numeric(length(thresholds))

    for (i in seq_along(thresholds)) {
      m <- calc_metrics(distH1, distH2, thresholds[i])
      TPR_vec[i] <- m$TPR
      FPR_vec[i] <- m$FPR
    }

    # Calculate AUC using trapezoidal rule
    ord <- order(FPR_vec)
    FPR_sorted <- FPR_vec[ord]
    TPR_sorted <- TPR_vec[ord]
    AUC_value <- sum(diff(FPR_sorted) * (TPR_sorted[-1] + TPR_sorted[-length(TPR_sorted)]) / 2)
    AUC_value <- abs(AUC_value)

    list(data = data.frame(th = thresholds, TPR = TPR_vec, FPR = FPR_vec),
         AUC = AUC_value)
  }

  # Calculate LR for sex
  calc_lr_sex <- function(MPs, obsSex, eps, propF) {
    propS <- c(F = propF, M = 1 - propF)
    if (MPs == obsSex) {
      LR <- (1 - eps) / propS[obsSex]
    } else {
      LR <- eps / propS[obsSex]
    }
    LR
  }

  # Calculate LR for age
  calc_lr_age <- function(MPa, MPr, obsAge, epa) {
    MPmin <- max(1, MPa - MPr)
    MPmax <- min(80, MPa + MPr)
    T1p <- (MPmax - MPmin) / 80
    T1p <- max(0.001, min(0.999, T1p))
    T0p <- 1 - T1p

    if (obsAge >= MPmin && obsAge <= MPmax) {
      LR <- (1 - epa) / T1p
    } else {
      LR <- epa / T0p
    }
    LR
  }

  # Calculate LR for hair color
  calc_lr_hair <- function(MPc, obsC, epc, propC) {
    # Uniform error model
    n_colors <- length(propC)
    if (MPc == obsC) {
      p_obs_h1 <- 1 - epc
    } else {
      p_obs_h1 <- epc / (n_colors - 1)
    }
    LR <- p_obs_h1 / propC[obsC]
    LR
  }

  # ============================================================================
  # UI DEFINITION
  # ============================================================================

  ui <- shiny::fluidPage(
    theme = shinythemes::shinytheme("flatly"),

    # Custom CSS for professional appearance
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        @import url('https://fonts.googleapis.com/css2?family=Source+Sans+Pro:wght@400;600;700&display=swap');

        body {
          font-family: 'Source Sans Pro', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
          background-color: #f8f9fa;
        }

        .main-header {
          background: linear-gradient(135deg, #1e3a5f 0%, #2c5282 100%);
          color: white;
          padding: 20px 30px;
          margin: -15px -15px 25px -15px;
          display: flex;
          align-items: center;
          box-shadow: 0 2px 10px rgba(0,0,0,0.15);
        }

        .main-header img {
          height: 60px;
          margin-right: 20px;
        }

        .main-header h1 {
          margin: 0;
          font-size: 28px;
          font-weight: 700;
          letter-spacing: 0.5px;
        }

        .main-header p {
          margin: 5px 0 0 0;
          opacity: 0.9;
          font-size: 14px;
        }

        .info-box {
          background: #ffffff;
          border-left: 4px solid #2c5282;
          padding: 15px 20px;
          margin: 15px 0;
          border-radius: 0 8px 8px 0;
          box-shadow: 0 1px 3px rgba(0,0,0,0.08);
        }

        .result-box {
          background: linear-gradient(135deg, #ebf8ff 0%, #bee3f8 100%);
          border: 2px solid #2c5282;
          padding: 25px;
          margin: 15px 0;
          border-radius: 12px;
          text-align: center;
        }

        .result-box h2 {
          margin: 10px 0;
          color: #1e3a5f;
          font-size: 32px;
          font-weight: 700;
        }

        .result-box h4 {
          margin: 0;
          color: #4a5568;
          font-weight: 600;
        }

        .warning-box {
          background: #fffaf0;
          border-left: 4px solid #dd6b20;
          padding: 15px 20px;
          margin: 15px 0;
          border-radius: 0 8px 8px 0;
        }

        .success-box {
          background: #f0fff4;
          border-left: 4px solid #38a169;
          padding: 15px 20px;
          margin: 15px 0;
          border-radius: 0 8px 8px 0;
        }

        .formula-box {
          background: #2d3748;
          color: #e2e8f0;
          padding: 15px 20px;
          font-family: 'Fira Code', 'Consolas', monospace;
          border-radius: 8px;
          margin: 15px 0;
          font-size: 14px;
        }

        .tutorial-section {
          background: white;
          padding: 25px;
          margin: 15px 0;
          border-radius: 12px;
          box-shadow: 0 2px 8px rgba(0,0,0,0.08);
        }

        .tutorial-section h4 {
          color: #1e3a5f;
          border-bottom: 2px solid #e2e8f0;
          padding-bottom: 10px;
          margin-bottom: 15px;
        }

        .nav-tabs {
          border-bottom: 2px solid #e2e8f0;
        }

        .nav-tabs > li > a {
          font-weight: 600;
          color: #4a5568;
          border: none;
          padding: 12px 20px;
        }

        .nav-tabs > li.active > a {
          color: #2c5282;
          border-bottom: 3px solid #2c5282;
          background: transparent;
        }

        .nav-tabs > li > a:hover {
          background: #f7fafc;
          border: none;
        }

        .well {
          background: white;
          border: 1px solid #e2e8f0;
          border-radius: 12px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        }

        .interpretation-text {
          font-size: 16px;
          line-height: 1.6;
          color: #2d3748;
        }

      "))
    ),

    # Header with logo
    shiny::div(class = "main-header",
      shiny::tags$img(src = "www/MispiIcon.png", alt = "mispitools logo"),
      shiny::div(
        shiny::h1("mispitools"),
        shiny::p("Likelihood Ratio Analysis for Missing Person Identification")
      )
    ),

    # Main tabs
    shiny::tabsetPanel(
      id = "main_tabs",
      type = "tabs",

      # ==========================================================================
      # TAB 1: WELCOME
      # ==========================================================================
      shiny::tabPanel(
        "Overview",
        shiny::fluidRow(
          shiny::column(10, offset = 1,
            shiny::br(),
            shiny::h2("Evidence Evaluation in Missing Person Cases"),
            shiny::p(class = "interpretation-text",
              "This application implements a Bayesian framework for evaluating non-genetic
              evidence in missing person identification. It calculates likelihood ratios (LRs)
              that quantify the weight of evidence, following established methodology in
              forensic genetics as implemented in software such as Familias (Kling et al., 2014)."
            ),
            shiny::hr(),

            shiny::fluidRow(
              shiny::column(6,
                shiny::div(class = "info-box",
                  shiny::h4("The Likelihood Ratio"),
                  shiny::p("The LR compares two hypotheses:"),
                  shiny::tags$ul(
                    shiny::tags$li(shiny::strong("H1:"), " The unidentified person is the missing person"),
                    shiny::tags$li(shiny::strong("H2:"), " The unidentified person is not the missing person")
                  ),
                  shiny::div(class = "formula-box",
                    "LR = P(Evidence | H1) / P(Evidence | H2)"
                  ),
                  shiny::p("The LR does ", shiny::strong("not"), " provide a probability of identity.
                           It quantifies how much the evidence changes the odds.")
                )
              ),
              shiny::column(6,
                shiny::div(class = "info-box",
                  shiny::h4("Using the LR"),
                  shiny::p("The LR quantifies how many times more likely the evidence is under H1 vs H2."),
                  shiny::tags$ul(
                    shiny::tags$li("LR > 1: Evidence more probable under H1"),
                    shiny::tags$li("LR < 1: Evidence more probable under H2"),
                    shiny::tags$li("LR = 1: Evidence equally probable under both hypotheses")
                  ),
                  shiny::p("The numeric LR value should be reported without verbal qualifiers.")
                )
              )
            ),

            shiny::h4("Available Evidence Types"),
            shiny::fluidRow(
              shiny::column(3,
                shiny::div(class = "well", style = "text-align: center; padding: 20px;",
                  shiny::h5("Biological Sex"),
                  shiny::p(class = "text-muted", "Male/female classification with observation error")
                )
              ),
              shiny::column(3,
                shiny::div(class = "well", style = "text-align: center; padding: 20px;",
                  shiny::h5("Age"),
                  shiny::p(class = "text-muted", "Age range matching with estimation uncertainty")
                )
              ),
              shiny::column(3,
                shiny::div(class = "well", style = "text-align: center; padding: 20px;",
                  shiny::h5("Hair Color"),
                  shiny::p(class = "text-muted", "Five-category classification with error matrix")
                )
              ),
              shiny::column(3,
                shiny::div(class = "well", style = "text-align: center; padding: 20px;",
                  shiny::h5("Birth Date"),
                  shiny::p(class = "text-muted", "Date discrepancy using Dirichlet model")
                )
              )
            ),

            shiny::hr(),
            shiny::div(style = "text-align: center; color: #718096; padding: 20px;",
              shiny::p(shiny::strong("References")),
              shiny::p("Marsico et al. (2023). FSI:Genetics 66:102891 |
                       Marsico et al. (2021). FSI:Genetics 52:102519"),
              shiny::p("Kling et al. (2014). Familias 3. FSI:Genetics 13:121-127")
            )
          )
        )
      ),

      # ==========================================================================
      # TAB 2: INDIVIDUAL EVIDENCE
      # ==========================================================================
      shiny::tabPanel(
        "Individual Evidence",
        shiny::fluidRow(
          shiny::column(12,
            shiny::br(),
            shiny::h3("Single Evidence Type Analysis"),
            shiny::p("Calculate the likelihood ratio for one type of evidence. All parameters
                     can be adjusted based on case-specific information.")
          )
        ),
        shiny::fluidRow(
          shiny::column(3,
            shiny::wellPanel(
              shiny::selectInput("indiv_type", "Evidence Type:",
                choices = c("Biological Sex" = "sex",
                            "Age" = "age",
                            "Hair Color" = "hair")),
              shiny::hr(),

              # Sex parameters
              shiny::conditionalPanel(
                condition = "input.indiv_type == 'sex'",
                shiny::h5("Case Parameters"),
                shiny::selectInput("indiv_MPs", "MP reported sex:",
                  c("Female" = "F", "Male" = "M")),
                shiny::selectInput("indiv_obsSex", "Observed sex:",
                  c("Female" = "F", "Male" = "M")),
                shiny::hr(),
                shiny::h5("Model Parameters"),
                shiny::sliderInput("indiv_eps_sex", "Observation error rate:",
                  0, 0.3, 0.05, 0.01),
                shiny::sliderInput("indiv_propF", "Female proportion in population:",
                  0.3, 0.7, 0.51, 0.01),
                shiny::p(class = "text-muted small",
                  "Error rate: probability of recording wrong sex.
                   Population proportion based on demographic data.")
              ),

              # Age parameters
              shiny::conditionalPanel(
                condition = "input.indiv_type == 'age'",
                shiny::h5("Case Parameters"),
                shiny::numericInput("indiv_MPa", "MP estimated age (years):", 35, 1, 90),
                shiny::numericInput("indiv_MPr", "Estimation uncertainty (years):", 5, 1, 20),
                shiny::numericInput("indiv_obsAge", "Observed age:", 37, 1, 90),
                shiny::hr(),
                shiny::h5("Model Parameters"),
                shiny::sliderInput("indiv_epa", "Age category error rate:",
                  0, 0.3, 0.05, 0.01),
                shiny::p(class = "text-muted small",
                  "The matching range is MP age +/- uncertainty.
                   Reference assumes uniform age distribution (1-80 years).")
              ),

              # Hair color parameters
              shiny::conditionalPanel(
                condition = "input.indiv_type == 'hair'",
                shiny::h5("Case Parameters"),
                shiny::selectInput("indiv_MPc", "MP hair color:",
                  choices = c("Black" = 1, "Brown" = 2, "Blonde" = 3,
                              "Red" = 4, "Gray/White" = 5), selected = 2),
                shiny::selectInput("indiv_obsC", "Observed hair color:",
                  choices = c("Black" = 1, "Brown" = 2, "Blonde" = 3,
                              "Red" = 4, "Gray/White" = 5), selected = 2),
                shiny::hr(),
                shiny::h5("Model Parameters"),
                shiny::sliderInput("indiv_epc", "Color misclassification rate:",
                  0, 0.3, 0.05, 0.01),
                shiny::p(class = "text-muted small", "Population frequencies (European):"),
                shiny::p(class = "text-muted small",
                  "Black: 8%, Brown: 45%, Blonde: 30%, Red: 7%, Gray: 10%")
              )
            )
          ),
          shiny::column(9,
            shiny::div(class = "result-box",
              shiny::h4("Likelihood Ratio"),
              shiny::h2(shiny::textOutput("indiv_lr_value")),
              shiny::p(shiny::textOutput("indiv_lr_log"))
            ),
            shiny::div(class = "formula-box",
              shiny::h5(style = "color: #a0aec0; margin-top: 0;", "Calculation Details"),
              shiny::verbatimTextOutput("indiv_formula")
            )
          )
        )
      ),

      # ==========================================================================
      # TAB 3: CPT ANALYSIS
      # ==========================================================================
      shiny::tabPanel(
        "Probability Tables",
        shiny::fluidRow(
          shiny::column(12,
            shiny::br(),
            shiny::h3("Conditional Probability Tables"),
            shiny::p("Visualization of evidence probabilities under each hypothesis.
                     The tables show joint probabilities for sex, age category, and hair color.")
          )
        ),
        shiny::fluidRow(
          shiny::column(2,
            shiny::wellPanel(
              shiny::h5("Missing Person"),
              shiny::selectInput("cpt_MPs", "Sex:", c("Female" = "F", "Male" = "M")),
              shiny::numericInput("cpt_MPa", "Age:", 40, 1, 90),
              shiny::numericInput("cpt_MPr", "Age range:", 6, 1, 20),
              shiny::selectInput("cpt_MPc", "Hair color:", 1:5, selected = 2)
            )
          ),
          shiny::column(2,
            shiny::wellPanel(
              shiny::h5("Error Rates"),
              shiny::numericInput("cpt_eps", "Sex:", 0.05, 0, 0.5, 0.01),
              shiny::numericInput("cpt_epa", "Age:", 0.05, 0, 0.5, 0.01),
              shiny::numericInput("cpt_epc", "Hair:", 0.05, 0, 0.5, 0.01)
            ),
            shiny::wellPanel(
              shiny::h5("Population"),
              shiny::sliderInput("cpt_propF", "Female %:", 0.3, 0.7, 0.51, 0.01)
            )
          ),
          shiny::column(8,
            shiny::plotOutput("cpt_plots", height = "380px"),
            shiny::div(class = "info-box",
              shiny::fluidRow(
                shiny::column(4,
                  shiny::strong("P(D|H2):"), " Phenotype probability in reference population"
                ),
                shiny::column(4,
                  shiny::strong("P(D|H1):"), " Probability if person IS the MP"
                ),
                shiny::column(4,
                  shiny::strong("log10(LR):"),
                  shiny::span(style = "color: #38a169;", " Green = supports H1,"),
                  shiny::span(style = "color: #e53e3e;", " Red = supports H2")
                )
              )
            )
          )
        )
      ),

      # ==========================================================================
      # TAB 4: DISTRIBUTION
      # ==========================================================================
      shiny::tabPanel(
        "LR Distribution",
        shiny::fluidRow(
          shiny::column(12,
            shiny::br(),
            shiny::h3("Likelihood Ratio Distributions"),
            shiny::p("Comparison of LR distributions under H1 and H2. The separation between
                     distributions indicates discriminating power.")
          )
        ),
        shiny::fluidRow(
          shiny::column(6,
            shiny::plotOutput("dist_plot", height = "350px")
          ),
          shiny::column(6,
            shiny::plotOutput("dist_roc", height = "350px")
          )
        ),
        shiny::fluidRow(
          shiny::column(6,
            shiny::div(class = "success-box",
              shiny::h5("Under H1 (Person IS the MP)"),
              shiny::verbatimTextOutput("dist_stats_h1")
            )
          ),
          shiny::column(6,
            shiny::div(class = "warning-box",
              shiny::h5("Under H2 (Person is NOT the MP)"),
              shiny::verbatimTextOutput("dist_stats_h2")
            )
          )
        )
      ),

      # ==========================================================================
      # TAB 5: COMBINE EVIDENCE
      # ==========================================================================
      shiny::tabPanel(
        "Combine Evidence",
        shiny::fluidRow(
          shiny::column(12,
            shiny::br(),
            shiny::h3("Evidence Combination"),
            shiny::p("Under conditional independence, individual LRs multiply to give a combined LR.
                     This assumption requires that evidence types are independent given identity.")
          )
        ),
        shiny::fluidRow(
          shiny::column(12,
            shiny::div(class = "warning-box",
              shiny::strong("Independence Assumption: "),
              "LR multiplication assumes conditional independence of evidence sources.
               Correlated evidence (e.g., hair and skin color) requires joint modeling."
            )
          )
        ),
        shiny::fluidRow(
          shiny::column(3,
            shiny::wellPanel(
              shiny::h5("Sex Evidence"),
              shiny::checkboxInput("comb_use_sex", "Include", TRUE),
              shiny::conditionalPanel(
                condition = "input.comb_use_sex",
                shiny::selectInput("comb_MPs", "MP:", c("F", "M")),
                shiny::selectInput("comb_obsSex", "Observed:", c("F", "M")),
                shiny::numericInput("comb_eps", "Error:", 0.05, 0, 0.3, 0.01)
              ),
              shiny::hr(),
              shiny::strong(shiny::textOutput("comb_lr_sex"))
            )
          ),
          shiny::column(3,
            shiny::wellPanel(
              shiny::h5("Age Evidence"),
              shiny::checkboxInput("comb_use_age", "Include", TRUE),
              shiny::conditionalPanel(
                condition = "input.comb_use_age",
                shiny::numericInput("comb_MPa", "MP age:", 35, 1, 90),
                shiny::numericInput("comb_MPr", "Range:", 5, 1, 20),
                shiny::numericInput("comb_obsAge", "Observed:", 37, 1, 90),
                shiny::numericInput("comb_epa", "Error:", 0.05, 0, 0.3, 0.01)
              ),
              shiny::hr(),
              shiny::strong(shiny::textOutput("comb_lr_age"))
            )
          ),
          shiny::column(3,
            shiny::wellPanel(
              shiny::h5("Hair Color Evidence"),
              shiny::checkboxInput("comb_use_hair", "Include", TRUE),
              shiny::conditionalPanel(
                condition = "input.comb_use_hair",
                shiny::selectInput("comb_MPc", "MP:", 1:5, selected = 2),
                shiny::selectInput("comb_obsC", "Observed:", 1:5, selected = 2),
                shiny::numericInput("comb_epc", "Error:", 0.05, 0, 0.3, 0.01)
              ),
              shiny::hr(),
              shiny::strong(shiny::textOutput("comb_lr_hair"))
            )
          ),
          shiny::column(3,
            shiny::div(class = "result-box", style = "min-height: 320px;",
              shiny::h4("Combined Result"),
              shiny::h2(shiny::textOutput("comb_lr_total")),
              shiny::p(shiny::textOutput("comb_lr_log")),
              shiny::hr(),
              shiny::verbatimTextOutput("comb_formula")
            )
          )
        )
      ),

      # ==========================================================================
      # TAB 6: DECISION ANALYSIS
      # ==========================================================================
      shiny::tabPanel(
        "Decision Analysis",
        shiny::fluidRow(
          shiny::column(12,
            shiny::br(),
            shiny::h3("Decision Threshold Selection"),
            shiny::p("Evaluate the trade-off between false positive and false negative rates
                     at different LR thresholds.")
          )
        ),
        shiny::fluidRow(
          shiny::column(4,
            shiny::wellPanel(
              shiny::h5("Threshold Selection"),
              shiny::sliderInput("dec_threshold", HTML("log<sub>10</sub>(LR) Threshold:"),
                                 min = -3, max = 3, value = 0, step = 0.1),
              shiny::hr(),
              shiny::h5("Performance Metrics"),
              shiny::verbatimTextOutput("dec_metrics")
            ),
            shiny::div(class = "info-box",
              shiny::h6("Metric Definitions"),
              shiny::tags$ul(class = "small",
                shiny::tags$li(shiny::strong("TPR:"), " True Positive Rate (sensitivity)"),
                shiny::tags$li(shiny::strong("FPR:"), " False Positive Rate (1 - specificity)"),
                shiny::tags$li(shiny::strong("MCC:"), " Matthews Correlation Coefficient")
              )
            )
          ),
          shiny::column(8,
            shiny::plotOutput("dec_roc_plot", height = "400px")
          )
        )
      ),

      # ==========================================================================
      # TAB 7: TUTORIAL
      # ==========================================================================
      shiny::tabPanel(
        "Tutorial",
        shiny::fluidRow(
          shiny::column(10, offset = 1,
            shiny::br(),
            shiny::h2("Understanding Likelihood Ratios in Forensic Science"),

            shiny::div(class = "tutorial-section",
              shiny::h4("1. The Bayesian Framework"),
              shiny::p("In forensic identification, we evaluate evidence by comparing
                       two competing propositions:"),
              shiny::tags$ul(
                shiny::tags$li(shiny::strong("H1 (Prosecution):"),
                  " The unidentified remains belong to the missing person"),
                shiny::tags$li(shiny::strong("H2 (Defense):"),
                  " The unidentified remains belong to someone else")
              ),
              shiny::p("The likelihood ratio quantifies the relative probability of
                       the evidence under each hypothesis:"),
              shiny::div(class = "formula-box",
                "LR = P(Evidence | H1) / P(Evidence | H2)"
              ),
              shiny::p(shiny::strong("Important:"), " The LR is NOT the probability
                       that the person is the missing individual. It tells us how
                       much the evidence changes our belief.")
            ),

            shiny::div(class = "tutorial-section",
              shiny::h4("2. Worked Example: Sex Evidence"),
              shiny::p("Consider a case where:"),
              shiny::tags$ul(
                shiny::tags$li("The missing person is reported as female"),
                shiny::tags$li("The unidentified remains are determined to be female"),
                shiny::tags$li("Observation error rate: 5% (probability of misclassification)"),
                shiny::tags$li("Population: 51% female")
              ),
              shiny::div(class = "formula-box",
                "P(observed female | H1) = 1 - 0.05 = 0.95\n",
                "P(observed female | H2) = 0.51\n",
                "LR = 0.95 / 0.51 = 1.86"
              ),
              shiny::p("The evidence is 1.86 times more probable under H1 than H2.")
            ),

            shiny::div(class = "tutorial-section",
              shiny::h4("3. Worked Example: Age Evidence"),
              shiny::p("Consider a case where:"),
              shiny::tags$ul(
                shiny::tags$li("Missing person's estimated age: 35 years"),
                shiny::tags$li("Age estimation uncertainty: +/- 5 years (range: 30-40)"),
                shiny::tags$li("Observed age from remains: 37 years (within range)"),
                shiny::tags$li("Reference population: uniform distribution 1-80 years")
              ),
              shiny::div(class = "formula-box",
                "P(age in range | H1) = 1 - error = 0.95\n",
                "P(age in range | H2) = 10/80 = 0.125\n",
                "LR = 0.95 / 0.125 = 7.6"
              ),
              shiny::p("The evidence is 7.6 times more probable under H1 than H2.")
            ),

            shiny::div(class = "tutorial-section",
              shiny::h4("4. Combining Evidence"),
              shiny::p("When evidence sources are conditionally independent given
                       identity, their LRs multiply:"),
              shiny::div(class = "formula-box",
                "LR_combined = LR_sex x LR_age x LR_hair\n",
                "Example: 1.86 x 7.6 x 3.2 = 45.2"
              ),
              shiny::div(class = "warning-box",
                shiny::strong("Caution:"), " The independence assumption may be
                violated for correlated traits. For example, hair color and skin
                color are genetically correlated. In such cases, joint modeling
                is required."
              )
            ),

            shiny::div(class = "tutorial-section",
              shiny::h4("5. References"),
              shiny::tags$ul(
                shiny::tags$li("Marsico FL et al. (2023). Likelihood ratios for
                  non-genetic evidence in missing person cases. ",
                  shiny::em("FSI:Genetics"), " 66:102891"),
                shiny::tags$li("Marsico FL et al. (2021). Making decisions in
                  missing person identification cases with low statistical power. ",
                  shiny::em("FSI:Genetics"), " 52:102519"),
                shiny::tags$li("Kling D, Tillmar AO, Egeland T (2014). Familias 3:
                  Extensions and new functionality. ",
                  shiny::em("FSI:Genetics"), " 13:121-127"),
                shiny::tags$li("Balding DJ, Steele CD (2015). ",
                  shiny::em("Weight-of-Evidence for Forensic DNA Profiles"),
                  ". 2nd ed. Wiley")
              )
            )
          )
        )
      )
    )
  )

  # ============================================================================
  # SERVER LOGIC
  # ============================================================================

  server <- function(input, output, session) {

    # --------------------------------------------------------------------------
    # INDIVIDUAL EVIDENCE TAB
    # --------------------------------------------------------------------------

    indiv_lr <- shiny::reactive({
      if (input$indiv_type == "sex") {
        calc_lr_sex(input$indiv_MPs, input$indiv_obsSex,
                    input$indiv_eps_sex, input$indiv_propF)
      } else if (input$indiv_type == "age") {
        calc_lr_age(input$indiv_MPa, input$indiv_MPr,
                    input$indiv_obsAge, input$indiv_epa)
      } else if (input$indiv_type == "hair") {
        propC <- DEFAULT_HAIR_FREQ
        calc_lr_hair(as.numeric(input$indiv_MPc), as.numeric(input$indiv_obsC),
                     input$indiv_epc, propC)
      } else {
        1
      }
    })

    output$indiv_lr_value <- shiny::renderText({
      sprintf("LR = %.3f", indiv_lr())
    })

    output$indiv_lr_log <- shiny::renderText({
      sprintf("log10(LR) = %.3f", log10(indiv_lr()))
    })

    output$indiv_formula <- shiny::renderPrint({
      if (input$indiv_type == "sex") {
        match_status <- if (input$indiv_MPs == input$indiv_obsSex) "match" else "mismatch"
        cat("Sex Evidence (", match_status, ")\n", sep = "")
        cat("----------------------------------------\n")
        if (match_status == "match") {
          cat("P(observed | H1) = 1 - eps =", 1 - input$indiv_eps_sex, "\n")
        } else {
          cat("P(observed | H1) = eps =", input$indiv_eps_sex, "\n")
        }
        obs_freq <- if (input$indiv_obsSex == "F") input$indiv_propF else (1 - input$indiv_propF)
        cat("P(observed | H2) =", round(obs_freq, 3), "\n")
        cat("LR =", round(indiv_lr(), 4), "\n")

      } else if (input$indiv_type == "age") {
        in_range <- input$indiv_obsAge >= (input$indiv_MPa - input$indiv_MPr) &&
                    input$indiv_obsAge <= (input$indiv_MPa + input$indiv_MPr)
        cat("Age Evidence (", if (in_range) "in range" else "out of range", ")\n", sep = "")
        cat("----------------------------------------\n")
        cat("MP age:", input$indiv_MPa, "+/-", input$indiv_MPr, "years\n")
        cat("Observed:", input$indiv_obsAge, "years\n")
        range_prob <- (2 * input$indiv_MPr) / 80
        if (in_range) {
          cat("P(in range | H1) = 1 - eps =", 1 - input$indiv_epa, "\n")
          cat("P(in range | H2) =", round(range_prob, 4), "\n")
        } else {
          cat("P(out of range | H1) = eps =", input$indiv_epa, "\n")
          cat("P(out of range | H2) =", round(1 - range_prob, 4), "\n")
        }
        cat("LR =", round(indiv_lr(), 4), "\n")

      } else if (input$indiv_type == "hair") {
        match_status <- if (input$indiv_MPc == input$indiv_obsC) "match" else "mismatch"
        cat("Hair Color Evidence (", match_status, ")\n", sep = "")
        cat("----------------------------------------\n")
        cat("MP:", hair_labels[input$indiv_MPc], "\n")
        cat("Observed:", hair_labels[input$indiv_obsC], "\n")
        cat("Population frequency:", DEFAULT_HAIR_FREQ[as.numeric(input$indiv_obsC)], "\n")
        cat("LR =", round(indiv_lr(), 4), "\n")
      }
    })

    # --------------------------------------------------------------------------
    # CPT ANALYSIS TAB
    # --------------------------------------------------------------------------

    cpt_pop <- shiny::reactive({
      build_cpt_pop(input$cpt_propF, input$cpt_MPa, input$cpt_MPr, DEFAULT_HAIR_FREQ)
    })

    cpt_mp <- shiny::reactive({
      build_cpt_mp(input$cpt_MPs, as.numeric(input$cpt_MPc),
                   input$cpt_eps, input$cpt_epa, input$cpt_epc)
    })

    output$cpt_plots <- shiny::renderPlot({
      create_cpt_plots(cpt_pop(), cpt_mp())
    }, res = 100)

    # --------------------------------------------------------------------------
    # DISTRIBUTION TAB
    # --------------------------------------------------------------------------

    lr_dists <- shiny::reactive({
      build_lr_distributions(cpt_pop(), cpt_mp())
    })

    output$dist_plot <- shiny::renderPlot({
      dists <- lr_dists()
      create_dist_plot(dists$H1, dists$H2)
    }, res = 100)

    output$dist_stats_h1 <- shiny::renderPrint({
      dists <- lr_dists()
      if (nrow(dists$H1) == 0) {
        cat("No data available\n")
        return()
      }
      x <- dists$H1$x
      w <- dists$H1$w
      cat("Weighted mean log10(LR):", sprintf("%.3f", sum(x * w)), "\n")
      cat("Range: [", sprintf("%.3f", min(x)), ",", sprintf("%.3f", max(x)), "]\n")
      cat("\nExpected evidence weight favors H1\n")
    })

    output$dist_stats_h2 <- shiny::renderPrint({
      dists <- lr_dists()
      if (nrow(dists$H2) == 0) {
        cat("No data available\n")
        return()
      }
      x <- dists$H2$x
      w <- dists$H2$w
      cat("Weighted mean log10(LR):", sprintf("%.3f", sum(x * w)), "\n")
      cat("Range: [", sprintf("%.3f", min(x)), ",", sprintf("%.3f", max(x)), "]\n")
      cat("\nExpected evidence weight under H2\n")
    })

    output$dist_roc <- shiny::renderPlot({
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
        ggplot2::geom_line(color = "#2c5282", linewidth = 1.2) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             color = "gray50", linewidth = 0.8) +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(title = sprintf("ROC Curve (AUC = %.3f)", auc_val),
                      subtitle = "Discrimination between H1 and H2",
                      x = "False Positive Rate (1 - Specificity)",
                      y = "True Positive Rate (Sensitivity)") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          panel.grid.minor = ggplot2::element_blank()
        )
    }, res = 100)

    # --------------------------------------------------------------------------
    # COMBINE EVIDENCE TAB
    # --------------------------------------------------------------------------

    comb_lr_sex_val <- shiny::reactive({
      if (input$comb_use_sex) {
        calc_lr_sex(input$comb_MPs, input$comb_obsSex, input$comb_eps, 0.51)
      } else {
        1
      }
    })

    comb_lr_age_val <- shiny::reactive({
      if (input$comb_use_age) {
        calc_lr_age(input$comb_MPa, input$comb_MPr, input$comb_obsAge, input$comb_epa)
      } else {
        1
      }
    })

    comb_lr_hair_val <- shiny::reactive({
      if (input$comb_use_hair) {
        calc_lr_hair(as.numeric(input$comb_MPc), as.numeric(input$comb_obsC),
                     input$comb_epc, DEFAULT_HAIR_FREQ)
      } else {
        1
      }
    })

    output$comb_lr_sex <- shiny::renderText({
      if (input$comb_use_sex) sprintf("LR = %.2f", comb_lr_sex_val()) else "Not included"
    })

    output$comb_lr_age <- shiny::renderText({
      if (input$comb_use_age) sprintf("LR = %.2f", comb_lr_age_val()) else "Not included"
    })

    output$comb_lr_hair <- shiny::renderText({
      if (input$comb_use_hair) sprintf("LR = %.2f", comb_lr_hair_val()) else "Not included"
    })

    comb_total <- shiny::reactive({
      comb_lr_sex_val() * comb_lr_age_val() * comb_lr_hair_val()
    })

    output$comb_lr_total <- shiny::renderText({
      sprintf("LR = %.2f", comb_total())
    })

    output$comb_lr_log <- shiny::renderText({
      sprintf("log10(LR) = %.3f", log10(comb_total()))
    })

    output$comb_formula <- shiny::renderPrint({
      parts <- c()
      vals <- c()
      if (input$comb_use_sex) { parts <- c(parts, "Sex"); vals <- c(vals, comb_lr_sex_val()) }
      if (input$comb_use_age) { parts <- c(parts, "Age"); vals <- c(vals, comb_lr_age_val()) }
      if (input$comb_use_hair) { parts <- c(parts, "Hair"); vals <- c(vals, comb_lr_hair_val()) }

      if (length(parts) == 0) {
        cat("No evidence selected\n")
        return()
      }

      cat(paste(sprintf("%.2f", vals), collapse = " x "), "\n")
      cat("=", sprintf("%.2f", prod(vals)), "\n")
    })

    # --------------------------------------------------------------------------
    # DECISION ANALYSIS TAB
    # --------------------------------------------------------------------------

    output$dec_metrics <- shiny::renderPrint({
      dists <- lr_dists()
      m <- calc_metrics(dists$H1, dists$H2, input$dec_threshold)

      cat("At threshold log10(LR) =", input$dec_threshold, "\n")
      cat("========================================\n\n")
      cat("True Positive Rate:  ", sprintf("%.3f", m$TPR), "\n")
      cat("False Positive Rate: ", sprintf("%.3f", m$FPR), "\n")
      cat("False Negative Rate: ", sprintf("%.3f", m$FNR), "\n")
      cat("True Negative Rate:  ", sprintf("%.3f", m$TNR), "\n")
      cat("\nMCC: ", sprintf("%.3f", m$MCC), "\n")
    })

    output$dec_roc_plot <- shiny::renderPlot({
      dists <- lr_dists()
      roc_obj <- build_roc(dists$H1, dists$H2)
      df_roc <- roc_obj$data
      auc_val <- roc_obj$AUC

      if (nrow(df_roc) == 0) {
        return(ggplot2::ggplot() +
                 ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Insufficient data") +
                 ggplot2::theme_void())
      }

      # Find current threshold point
      th_idx <- which.min(abs(df_roc$th - input$dec_threshold))
      current_point <- df_roc[th_idx, ]

      ggplot2::ggplot(df_roc, ggplot2::aes(x = FPR, y = TPR)) +
        ggplot2::geom_line(color = "#2c5282", linewidth = 1.5) +
        ggplot2::geom_point(data = current_point, ggplot2::aes(x = FPR, y = TPR),
                           color = "#e53e3e", size = 6, shape = 16) +
        ggplot2::geom_point(data = current_point, ggplot2::aes(x = FPR, y = TPR),
                           color = "white", size = 3, shape = 16) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             color = "gray50", linewidth = 0.8) +
        ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(title = sprintf("ROC Curve (AUC = %.3f)", auc_val),
                      subtitle = sprintf("Red point: current threshold (log10(LR) = %.1f)",
                                        input$dec_threshold),
                      x = "False Positive Rate", y = "True Positive Rate") +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          panel.grid.minor = ggplot2::element_blank()
        )
    }, res = 100)
  }

  # Add resource path for logo
  shiny::addResourcePath("www", system.file("www", package = "mispitools"))

  # Launch the app
  shiny::shinyApp(ui = ui, server = server)
}
