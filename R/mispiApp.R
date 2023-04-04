#' Missing person shiny app
#'
#' @importFrom shiny sidebarLayout sidebarPanel numericInput selectInput sliderInput renderPlot shinyApp
#' @import ggplot2
#' @export
#' @return An user interface for computing non-genetic LRs and conditioned probability tables.
#' @examples
#' CPT_MP()

mispiApp <- function() {
# UI
###spinner:
options(spinner.color = "#5661f4", spinner.type = 6, spinner.color.background = "#ffffff", spinner.size = 0.5)


# User interface
ui = shiny::fluidPage(

    shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::numericInput("MPa", "MP age:", 40, min = 0),
      shiny::numericInput("MPr", "MP age error range:", 6, min = 0),
      shiny::numericInput("MPc", "MP hair color (from 1 to 5):", 1, min = 1, max = 5),
      shiny::numericInput("eps", "Epsilon (error rate) for sex:", 0.05, min = 0, max = 1),
      shiny::numericInput("epa", "Epsilon (error rate) for age:", 0.05, min = 0, max = 1),
      shiny::numericInput("epc", "Epsilon (error rate) for hair color:", 0.02, min = 0, max = 1),
      shiny::selectInput("MPs", "MP sex:",
                  choices = c("F", "M")),
      shiny::sliderInput("propF", "Female proportion:",
                  min = 0, max = 1, value = 0.5, step = 0.1),
      shiny::sliderInput("propC1", "Hair color 1 proportion:",
                  min = 0, max = 1, value = 0.3, step = 0.05),
      shiny::sliderInput("propC2", "Hair color 2 proportion:",
                  min = 0, max = 1, value = 0.2, step = 0.05),
      shiny::sliderInput("propC3", "Hair color 3 proportion:",
                  min = 0, max = 1, value = 0.25, step = 0.05),
      shiny::sliderInput("propC4", "Hair color 4 proportion:",
                  min = 0, max = 1, value = 0.15, step = 0.05),
      shiny::sliderInput("propC5", "Hair color 5 proportion:",
                  min = 0, max = 1, value = 0.1, step = 0.05)
    ),
    shiny::mainPanel(
      shiny::plotOutput("myplot")
    )
  )
)

# Server
server <- function(input, output) {

  MainPlot <- function(propF = 0.5, MPa = 40, MPr = 6, propC = c(0.3,0.2, 0.25, 0.15,0.1), MPs = "F", MPc = 1, eps = 0.05, epa = 0.05, epc = 0.02){

    #CPT POP
    Age <- seq(1:80)
    MPmin <- MPa - MPr
    MPmax <- MPa + MPr
    T1p <- (MPmax-MPmin)/length(Age)  # Para una uniforme
    T0p <-  1-T1p
    propS = c(propF, 1-propF)

    jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")
    jointprob <- c(propS[1]*T1p, propS[1]*T0p, propS[2]*T1p, propS[2]*T0p)
    names(jointprob) <- jointname

    CPT_POP <- outer(jointprob,propC)


    # CVModel (Only Equal)
    ep12 <- ep13 <- ep14 <- ep15 <- ep23 <- ep24 <- ep25 <- ep34 <- ep35 <- ep45 <- ep <- epc

    l1 = 1/(1+ep12+ep13+ep14+ep15)
    l2 = 1/(1+ep12+ep23+ep24+ep25)
    l3 = 1/(1+ep13+ep23+ep34+ep35)
    l4 = 1/(1+ep14+ep24+ep34+ep45)
    l5 = 1/(1+ep15+ep25+ep35+ep45)

    errorMat <- rbind(c(l1, l1*ep12, l1*ep13, l1*ep14, l1*ep15), c(l2*ep12,l2,l2*ep13,l2*ep14,l2*ep15), c(l3*ep13, l3*ep23, l3, l3*ep34, l3*ep35), c(l4*ep14, l4*ep24, l4*ep34, l4, l4*ep45), c(l5*ep15, l5*ep25, l5*ep35, l5*ep45, l5))

    #CPT MP

    jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")

    if (MPs == "F") {jointprob <- c((1-eps)*(1-epa), (1-eps)*epa, eps*(1-epa), eps*epa)}
    else if  (MPs == "M") {jointprob <- c(eps*(1-epa), eps*epa,(1-eps)*(1-epa), (1-eps)*epa)}

    names(jointprob) <- jointname

    Col <- c(1,2,3,4,5)
    probC = errorMat[MPc,]
    names(probC) <- Col

    CPT_MP <- outer(jointprob,probC)

    graphics::par(mfrow = c(2, 1), mar = c(2, 4, 4, 2))
    POP <- reshape2::melt(CPT_POP)
    Var1 <- Var2  <- value <- NULL

    p1 <- ggplot2::ggplot(POP, aes(x = Var2, y = Var1)) +
      ggplot2::geom_raster(aes(fill=value)) +
      scale_fill_gradient(low="grey90", high="blue", limits = c(0,1)) +
      labs(x="Hair colour (C)", y="Biologocial sex-Age", title="P(D|H2)", limits = c(0,1)) +
      theme_bw() + theme(axis.text.x=element_text(size=13, angle=0, vjust=0.3),
                         axis.text.y=element_text(size=13),
                         plot.title=element_text(size=13)) +
      geom_text(aes(label = as.character(format(round(value, 3), nsmall = 3))))
    MP <- reshape2::melt(CPT_MP)
    p2 <- ggplot2::ggplot(MP, aes(x = Var2, y = Var1)) +
      ggplot2::geom_raster(aes(fill=value)) +
      scale_fill_gradient(low="grey90", high="blue", limits = c(0,1)) +
      labs(x="Hair colour (C)", y="Biologocial sex-Age", title="P(D|H1)") +
      theme_bw() + theme(axis.text.x=element_text(size=13, angle=0, vjust=0.3),
                         axis.text.y=element_text(size=13),
                         plot.title=element_text(size=13)) +
      geom_text(aes(label = as.character(format(round(value, 3), nsmall = 3))))
    LRtable <- log10(CPT_MP/CPT_POP)
    LR <- reshape2::melt(LRtable)
    p3 <- ggplot2::ggplot(LR, aes(x = Var2, y = Var1)) +
      ggplot2::geom_raster(aes(fill=value)) +
      scale_fill_gradient(low="grey90", high="blue") +
      labs(x="Hair colour (C)", y="Biologocial sex-Age", title="Log10(LR)") +
      theme_bw() + theme(axis.text.x=element_text(size=13, angle=0, vjust=0.3),
                         axis.text.y=element_text(size=13),
                         plot.title=element_text(size=13)) +
      geom_text(aes(label = as.character(format(round(value, 2), nsmall = 2))))

    p <- ((p1+p2+p3)) + patchwork::plot_annotation(tag_levels = 'A')
    p
    return(p)}


  output$myplot <- shiny::renderPlot({
    p <- MainPlot(propF = input$propF,
                  MPa = input$MPa,
                  MPr = input$MPr,
                  propC = c(input$propC1, input$propC2, input$propC3, input$propC4, input$propC5),
                  MPs = input$MPs,
                  MPc = input$MPc,
                  eps = input$eps,
                  epa = input$epa,
                  epc = input$epc)
    print(p)
  })
}

# Run the app
shiny::shinyApp(ui = ui, server = server)
}
