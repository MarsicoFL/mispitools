#' Missing person shiny app (versión mejorada)
# Minor formatting update for documentation clarity

#' @import shiny
#' @import shinythemes
#' @import ggplot2
#' @export
#' 
#' @return Interfaz con pestañas (tabs) para computar LRs y tablas de probabilidad condicionadas (CPT),
#' incluyendo una sección descriptiva, más plots y estética mejorada.
#' 
#' @examples
#' CPT_MP_improved()

library(shiny)
library(shinythemes)
library(ggplot2)
install.packages("reshape2")
install.packages("pROC")  # para calcular AUC fácilmente

#' @param None
#' Context: Launches the interactive Shiny app for LR comparison and ROC visualization.
#' Users can adjust missing-person parameters and explore calibration metrics.

#' @import shiny
#' @import shinythemes
#' @import ggplot2
#' @import reshape2
#' @import patchwork
#' @import pROC
#' @export
lrComparisonApp <- function() {
  
  #---------------------------------------------------
  # 1) Auxiliary functions
  #---------------------------------------------------
  
  library(shiny)
  library(shinythemes)
  library(ggplot2)
  library(reshape2)
  library(patchwork)
  library(pROC)
  
  # (a) Build CPT under H2 (Population)
  cpt_pop_function <- function(propF, MPa, MPr, propC) {
    Age <- seq(1:80)
    MPmin <- MPa - MPr
    MPmax <- MPa + MPr
    
    # Probability age in [MPmin..MPmax]
    T1p <- (MPmax - MPmin) / length(Age)
    if(T1p < 0)  T1p <- 0
    if(T1p > 1)  T1p <- 1
    T0p <- 1 - T1p
    
    # Probability sex = F or M
    propS <- c(propF, 1 - propF)
    
    # 4 possible combos: F-T1, F-T0, M-T1, M-T0
    jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")
    jointprob <- c(propS[1]*T1p, propS[1]*T0p, propS[2]*T1p, propS[2]*T0p)
    names(jointprob) <- jointname
    
    # Expand to hair color
    CPT_POP <- outer(jointprob, propC)
    rownames(CPT_POP) <- jointname
    colnames(CPT_POP) <- paste0("C", seq_along(propC))
    
    return(CPT_POP)
  }
  
  # (b) Build CPT under H1 (MP)
  cpt_mp_function <- function(MPs, MPc, eps, epa, epc) {
    jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")
    if (MPs == "F") {
      jointprob <- c((1-eps)*(1-epa), (1-eps)*epa, eps*(1-epa), eps*epa)
    } else {
      jointprob <- c(eps*(1-epa), eps*epa, (1-eps)*(1-epa), (1-eps)*epa)
    }
    names(jointprob) <- jointname
    
    # Error matrix for hair color (simplified model)
    ep12 <- ep13 <- ep14 <- ep15 <- ep23 <- ep24 <- ep25 <- ep34 <- ep35 <- ep45 <- epc
    l1 = 1/(1+ep12+ep13+ep14+ep15)
    l2 = 1/(1+ep12+ep23+ep24+ep25)
    l3 = 1/(1+ep13+ep23+ep34+ep35)
    l4 = 1/(1+ep14+ep24+ep34+ep45)
    l5 = 1/(1+ep15+ep25+ep35+ep45)
    
    errorMat <- rbind(
      c(l1,       l1*ep12,  l1*ep13,  l1*ep14,  l1*ep15),
      c(l2*ep12,  l2,       l2*ep13,  l2*ep14,  l2*ep15),
      c(l3*ep13,  l3*ep23,  l3,       l3*ep34,  l3*ep35),
      c(l4*ep14,  l4*ep24,  l4*ep34,  l4,       l4*ep45),
      c(l5*ep15,  l5*ep25,  l5*ep35,  l5*ep45,  l5)
    )
    rownames(errorMat) <- paste0("C", 1:5)
    colnames(errorMat) <- paste0("C", 1:5)
    
    # Probability of color = row MPc
    probC <- errorMat[MPc, ]
    
    # Expand
    CPT_MP <- outer(jointprob, probC)
    rownames(CPT_MP) <- jointname
    colnames(CPT_MP) <- paste0("C", 1:5)
    
    return(CPT_MP)
  }
  
  # (c) Main 3 plots: P(D|H2), P(D|H1), log10(LR)
  main_3_plots <- function(CPT_POP, CPT_MP) {
    dfPOP <- melt(CPT_POP)
    dfMP  <- melt(CPT_MP)
    dfLR  <- melt(log10(CPT_MP / CPT_POP))
    
    base_theme <- theme_minimal(base_size = 16) +
      theme(panel.grid = element_blank())
    
    p1 <- ggplot(dfPOP, aes(x=Var2, y=Var1, fill=value)) +
      geom_raster() +
      scale_fill_gradient(low="white", high="blue", limits=c(0,1)) +
      geom_text(aes(label=formatC(value, digits=3, format="f")), size=4) +
      labs(title="P(D|H2)", x="Hair colour", y="Sex - Age", fill="Prob") +
      base_theme
    
    p2 <- ggplot(dfMP, aes(x=Var2, y=Var1, fill=value)) +
      geom_raster() +
      scale_fill_gradient(low="white", high="blue", limits=c(0,1)) +
      geom_text(aes(label=formatC(value, digits=3, format="f")), size=4) +
      labs(title="P(D|H1)", x="Hair colour", y="Sex - Age", fill="Prob") +
      base_theme
    
    p3 <- ggplot(dfLR, aes(x=Var2, y=Var1, fill=value)) +
      geom_raster() +
      scale_fill_gradient(low="yellow", high="red") +
      geom_text(aes(label=formatC(value, digits=2, format="f")), size=4) +
      labs(title="log10(LR)", x="Hair colour", y="Sex - Age", fill="log10(LR)") +
      base_theme
    
    # Layout horizontally
    out <- (p1 + p2 + p3) + plot_layout(ncol=3)
    return(out)
  }
  
  # (d) Build distributions (logLR, with H1/H2 weights)
  build_lr_distributions <- function(CPT_POP, CPT_MP) {
    pop_vec <- as.vector(CPT_POP)
    mp_vec  <- as.vector(CPT_MP)
    lr_vec  <- mp_vec / pop_vec
    lr_vec[pop_vec == 0] <- NA
    log_lr  <- log10(lr_vec)
    
    df_H1 <- data.frame(x=log_lr, w=mp_vec)
    df_H2 <- data.frame(x=log_lr, w=pop_vec)
    df_H1 <- df_H1[!is.na(df_H1$x), ]
    df_H2 <- df_H2[!is.na(df_H2$x), ]
    
    df_H1$w <- df_H1$w / sum(df_H1$w)
    df_H2$w <- df_H2$w / sum(df_H2$w)
    
    return(list(H1=df_H1, H2=df_H2))
  }
  
  # (e) Overlapped histogram for H1 & H2
  distribution_barplots <- function(distH1, distH2, binwidth=0.5) {
    dH1 <- data.frame(logLR=distH1$x, w=distH1$w, Condition="H1")
    dH2 <- data.frame(logLR=distH2$x, w=distH2$w, Condition="H2")
    dAll <- rbind(dH1, dH2)
    
    base_theme <- theme_minimal(base_size = 16) +
      theme(panel.grid = element_blank())
    
    p <- ggplot(dAll, aes(x=logLR, weight=w, fill=Condition)) +
      geom_histogram(binwidth=binwidth, position="identity", alpha=0.5, color="black") +
      scale_fill_manual(values=c("H1"="blue", "H2"="red")) +
      labs(
        title="Distribution of log10(LR) - Overlapped (H1 & H2)",
        x="log10(LR)",
        y="Probability"
      ) +
      base_theme
    return(p)
  }
  
  # (f) Threshold-based metrics (TPR, FPR, TNR, FNR, MCC)
  calc_metrics_threshold <- function(distH1, distH2, threshold) {
    x1 <- distH1$x
    w1 <- distH1$w
    x2 <- distH2$x
    w2 <- distH2$w
    
    TP <- sum(w1[x1 >= threshold])
    FN <- sum(w1[x1 < threshold])
    FP <- sum(w2[x2 >= threshold])
    TN <- sum(w2[x2 < threshold])
    
    denom <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    if(denom == 0) {
      MCC <- NA
    } else {
      MCC <- (TP*TN - FP*FN) / denom
    }
    out <- list(TPR=TP, FPR=FP, TNR=TN, FNR=FN, MCC=MCC)
    return(out)
  }
  
  # (g) Build ROC + AUC
  build_roc <- function(distH1, distH2) {
    all_vals <- sort(unique(c(distH1$x, distH2$x)))
    if(length(all_vals) < 2) {
      # Edge case: everything might be the same. Just return a dummy.
      return(list(data=data.frame(th=numeric(0), TPR=numeric(0), FPR=numeric(0)), AUC=NA))
    }
    minx <- min(all_vals) - 0.01
    maxx <- max(all_vals) + 0.01
    thresholds <- seq(minx, maxx, length.out=200)
    
    TPR_vec <- numeric(length(thresholds))
    FPR_vec <- numeric(length(thresholds))
    
    for (i in seq_along(thresholds)) {
      m <- calc_metrics_threshold(distH1, distH2, thresholds[i])
      TPR_vec[i] <- m$TPR
      FPR_vec[i] <- m$FPR
    }
    AUC_value <- pROC::auc(FPR_vec, TPR_vec, type="linear")
    df_roc <- data.frame(th=thresholds, TPR=TPR_vec, FPR=FPR_vec)
    return(list(data=df_roc, AUC=AUC_value))
  }
  
  # (h) Helper to detect T1/T0 from an observed age wrt MPa ± MPr
  range_label <- function(obsAge, MPa, MPr) {
    mn <- MPa - MPr
    mx <- MPa + MPr
    if(obsAge >= mn && obsAge <= mx) return("T1")
    return("T0")
  }
  
  
  #---------------------------------------------------
  # 2) UI
  #---------------------------------------------------
  ui <- fluidPage(
    theme = shinytheme("yeti"),
    
    # Slightly larger text for labels/inputs
    tags$style(HTML("
      .shiny-input-panel { font-size: 16px; }
      label { font-size: 16px; }
    ")),
    
    titlePanel("LR Comparison for a Missing Person (MP)"),
    
    fluidRow(
      column(
        width = 12,
        h4("Welcome to the LR comparison app"),
        p("This application evaluates the likelihood of certain traits (age, sex, hair color) 
           in a missing person (MP) under two hypotheses: 
           H1 = 'It is the sought individual' vs. H2 = 'It is not.'  
           We compute Probability Tables (CPT) and the Likelihood Ratio (LR), 
           along with its distribution and calibration metrics."),
        hr()
      )
    ),
    
    tabsetPanel(
      
      #======== 1) Parameters & Main Plots ==============
      tabPanel(
        "Parameters & Main Plots",
        fluidRow(
          # Left panel (narrow)
          column(
            width = 2,
            wellPanel(
              h4("MP Parameters"),
              numericInput("MPa", "MP age:", 40, min=0),
              numericInput("MPr", "Age error range:", 6, min=0),
              numericInput("MPc", "MP hair color (1..5):", 1, min=1, max=5),
              selectInput("MPs", "MP sex:", choices=c("F","M")),
              hr(),
              h4("Error rates"),
              numericInput("eps", "Sex error:", 0.05, 0, 1),
              numericInput("epa", "Age error:", 0.05, 0, 1),
              numericInput("epc", "Hair color error:", 0.02, 0, 1)
            )
          ),
          
          # Center panel for the main plots
          column(
            width = 8,
            h4("Comparing Features"),
            p("We visualize P(D|H2), P(D|H1), and log10(LR)."),
            plotOutput("main_plots", height="400px")
          ),
          
          # Right panel (narrow)
          column(
            width = 2,
            wellPanel(
              h4("Population Parameters"),
              sliderInput("propF", "Prop. F:", 0, 1, 0.5, step=0.1),
              sliderInput("propC1", "Prop. C1:", 0, 1, 0.3, step=0.05),
              sliderInput("propC2", "Prop. C2:", 0, 1, 0.2, step=0.05),
              sliderInput("propC3", "Prop. C3:", 0, 1, 0.25, step=0.05),
              sliderInput("propC4", "Prop. C4:", 0, 1, 0.15, step=0.05),
              sliderInput("propC5", "Prop. C5:", 0, 1, 0.1, step=0.05)
            )
          )
        )
      ),
      
      #======== 2) Distribution & Calibration ============
      tabPanel(
        "Distribution & Calibration",
        fluidRow(
          column(
            width=12,
            br(),
            p("Below is the overlapped histogram of log10(LR) for H1 (blue) 
               and H2 (red). Then, a set of calibration metrics (TPR, FPR, etc.) 
               using a chosen threshold, and finally a ROC curve with AUC."),
            plotOutput("distPlot", height="350px")
          )
        ),
        fluidRow(
          column(
            width = 4,
            wellPanel(
              h4("Weighted Summary (H1 & H2)"),
              verbatimTextOutput("summary_metrics"),
              hr(),
              sliderInput("thresholdLR", "Threshold (log10(LR)):", min=-3, max=3, value=0, step=0.1),
              verbatimTextOutput("confusion_metrics")
            )
          ),
          column(
            width = 8,
            h4("ROC Curve"),
            plotOutput("rocPlot", height="350px"),
            p("AUC is shown in the plot title.")
          )
        )
      ),
      
      #======== 3) Expected LR ===========================
      tabPanel(
        "Expected LR",
        fluidRow(
          column(
            width=4,
            wellPanel(
              h4("Observed Person"),
              numericInput("obsAge", "Observed age:", 40, min=0),
              selectInput("obsSex", "Observed sex:", choices=c("F","M")),
              numericInput("obsColor", "Observed hair color (1..5):", 2, min=1, max=5),
              hr(),
              verbatimTextOutput("lr_expected")
            )
          ),
          column(
            width=8,
            br(),
            p("We compute the LR for the specific cell: [obsSex - T1/T0, obsColor]. 
               The T1/T0 label depends on whether 'observed age' lies within the 
               MP's age range (MPa ± MPr). If the population probability is zero 
               in that cell, the LR is undefined (NA).")
          )
        )
      )
    )
  )
  
  #---------------------------------------------------
  # 3) SERVER
  #---------------------------------------------------
  server <- function(input, output, session) {
    
    # Reactive: CPT_POP & CPT_MP
    cpt_pop_react <- reactive({
      cpt_pop_function(
        propF = input$propF,
        MPa   = input$MPa,
        MPr   = input$MPr,
        propC = c(input$propC1, input$propC2, input$propC3, input$propC4, input$propC5)
      )
    })
    
    cpt_mp_react <- reactive({
      cpt_mp_function(
        MPs = input$MPs,
        MPc = input$MPc,
        eps = input$eps,
        epa = input$epa,
        epc = input$epc
      )
    })
    
    # Main plots
    output$main_plots <- renderPlot({
      pop <- cpt_pop_react()
      mp  <- cpt_mp_react()
      main_3_plots(pop, mp)
    })
    
    # Distribution plot
    output$distPlot <- renderPlot({
      pop <- cpt_pop_react()
      mp  <- cpt_mp_react()
      dists <- build_lr_distributions(pop, mp)
      distribution_barplots(dists$H1, dists$H2, binwidth=0.5)
    })
    
    # For metrics & ROC, build distributions reactive
    lrDist <- reactive({
      pop <- cpt_pop_react()
      mp  <- cpt_mp_react()
      build_lr_distributions(pop, mp)
    })
    
    # Weighted summary stats
    output$summary_metrics <- renderPrint({
      d <- lrDist()
      xH1 <- d$H1$x
      wH1 <- d$H1$w
      xH2 <- d$H2$x
      wH2 <- d$H2$w
      
      meanH1 <- sum(xH1 * wH1)
      meanH2 <- sum(xH2 * wH2)
      
      weighted_median <- function(x, w) {
        o <- order(x)
        xx <- x[o]
        ww <- w[o]
        cum <- cumsum(ww)
        idx <- which(cum >= 0.5)[1]
        xx[idx]
      }
      
      medH1 <- weighted_median(xH1, wH1)
      medH2 <- weighted_median(xH2, wH2)
      
      cat("Under H1:\n")
      cat("  Weighted mean of log10(LR):", round(meanH1, 3), "\n")
      cat("  Weighted median of log10(LR):", round(medH1, 3), "\n\n")
      cat("Under H2:\n")
      cat("  Weighted mean of log10(LR):", round(meanH2, 3), "\n")
      cat("  Weighted median of log10(LR):", round(medH2, 3), "\n")
    })
    
    # Metrics at a chosen threshold
    output$confusion_metrics <- renderPrint({
      d <- lrDist()
      thr <- input$thresholdLR
      m   <- calc_metrics_threshold(d$H1, d$H2, thr)
      
      cat("Threshold (log10(LR)):", thr, "\n")
      cat("TPR =", round(m$TPR, 3), "\n")
      cat("FPR =", round(m$FPR, 3), "\n")
      cat("TNR =", round(m$TNR, 3), "\n")
      cat("FNR =", round(m$FNR, 3), "\n")
      cat("MCC =", round(m$MCC, 3), "\n")
    })
    
    # ROC Plot
    output$rocPlot <- renderPlot({
      d <- lrDist()
      roc_obj <- build_roc(d$H1, d$H2)
      df_roc  <- roc_obj$data
      auc_val <- roc_obj$AUC
      
      base_theme <- theme_minimal(base_size=16) +
        theme(panel.grid=element_blank())
      
      ggplot(df_roc, aes(x=FPR, y=TPR)) +
        geom_line(color="blue", size=1.2) +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="gray50") +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
        labs(
          title = paste0("ROC Curve - AUC: ", round(auc_val, 3)),
          x = "False Positive Rate",
          y = "True Positive Rate"
        ) +
        base_theme
    })
    
    # Expected LR
    output$lr_expected <- renderPrint({
      pop <- cpt_pop_react()
      mp  <- cpt_mp_react()
      
      # T1 or T0
      tObs <- range_label(input$obsAge, input$MPa, input$MPr)
      row_name <- paste0(input$obsSex, "-", tObs)
      col_name <- paste0("C", input$obsColor)
      
      if(! row_name %in% rownames(pop) || ! col_name %in% colnames(pop)) {
        cat("Combination out of range.\n")
        return(NULL)
      }
      pop_val <- pop[row_name, col_name]
      mp_val  <- mp[row_name, col_name]
      
      if(is.na(pop_val) || pop_val == 0) {
        cat("Population probability is 0 => LR is undefined (NA)\n")
      } else {
        LR <- mp_val / pop_val
        cat("LR =", round(LR, 5), "\n")
        cat("log10(LR) =", round(log10(LR), 3), "\n")
      }
    })
  }
  
  shinyApp(ui=ui, server=server)
}

lrComparisonApp()
