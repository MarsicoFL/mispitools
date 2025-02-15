#' postSim: A function for simulating posterior odds 
#'
#' @param datasim Output from simLRgen function.
#' @param PriorModel Prior odds model: "prelim" is based on preliminary data, and "uniform" uses only the prior probability of H1
#' @param Prior Prior probability for H1
#' @param eps epsilon parameter sex
#' @param erRs error parameter sex
#' @param epc epsilon parameter hair color
#' @param erRc error parameter hair color
#' @param MPc Missing person hair color
#' @param epa epsilon parameter age
#' @param erRa error parameter age
#' @param MPa Missing person age
#' @param MPr Missing person age error range
#' @export
#' @return A value of posterior odds.
#' @examples
#' library(forrel)
#' x = linearPed(2)
#' plot(x)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' datasim = simLRgen(x, missing = 5, 10, 123)
#' postSim(datasim)
postSim <- function(datasim, Prior = 0.01, PriorModel = c("prelim","uniform")[1], 
                   eps = 0.05, erRs = 0.01, epc = Cmodel(), erRc = Cmodel(), 
                   MPc = 1, epa = 0.05, erRa = 0.01, MPa = 10, MPr = 2) {
  
  unrelated_values <- vector()
  related_values <- vector()
  list_length <- length(datasim[["Unrelated"]])
  
  for (i in 1:list_length) {
    unrelated_value <- datasim[["Unrelated"]][[i]][["LRtotal"]][["H1:H2"]]
    related_value <- datasim[["Related"]][[i]][["LRtotal"]][["H1:H2"]]
    
    unrelated_values <- c(unrelated_values, unrelated_value)
    related_values <- c(related_values, related_value)
  }
  
  results_df <- data.frame(Unrelated = unrelated_values, Related = related_values)
  datasim <- results_df
  nsims = nrow(datasim)
  
  if (PriorModel == "prelim") {
    LRgenH2 <- datasim$Unrelated
    
    # simulaciones para H2
    LRs2 <- LRsex(LR = TRUE, H = 2, nsims = nsims, eps = eps, erRs = erRs)
    LRc2 <- LRcol(LR = TRUE, H = 2, nsims = nsims, epc = epc, erRc = erRc, MPc = MPc)
    LRa2 <- LRage(LR = TRUE, H = 2, nsims = nsims, epa = epa, erRa = erRa, MPa = MPa, MPr = MPr)
    
    datasim2 <- data.frame(
      LRs = LRs2$LRs,
      LRc = LRc2$LRc,
      LRa = LRa2$LRa
    )
    
    datasim2$LRtot <- datasim2$LRa * datasim2$LRc * datasim2$LRs
    datasim2$Label <- "H2"
    datasim2$LRgenH2 <- LRgenH2
    datasim2$PriorOdds <- Prior/(1-Prior)
    datasim2$PostOdds <- datasim2$PriorOdds * datasim2$LRtot
    datasim2$PostOddsGen <- datasim2$PostOdds * datasim2$LRgenH2
    
    # simulaciones para H1
    LRgenH1 <- datasim$Related
    
    LRs1 <- LRsex(LR = TRUE, H = 1, nsims = nsims, eps = eps, erRs = erRs)
    LRc1 <- LRcol(LR = TRUE, H = 1, nsims = nsims, epc = epc, erRc = erRc, MPc = MPc)
    LRa1 <- LRage(LR = TRUE, H = 1, nsims = nsims, epa = epa, erRa = erRa, MPa = MPa, MPr = MPr)
    
    datasim1 <- data.frame(
      LRs = LRs1$LRs,
      LRc = LRc1$LRc,
      LRa = LRa1$LRa
    )
    
    datasim1$LRtot <- datasim1$LRa * datasim1$LRc * datasim1$LRs
    datasim1$Label <- "H1"
    datasim1$LRgenH1 <- LRgenH1
    datasim1$PriorOdds <- Prior/(1-Prior)
    datasim1$PostOdds <- datasim1$PriorOdds * datasim1$LRtot
    datasim1$PostOddsGen <- datasim1$PostOdds * datasim1$LRgenH1
    
  } else {
    LRgenH2 <- datasim$Unrelated
    datasim2 <- data.frame(Unrelated = LRgenH2)
    datasim2$PriorOdds <- Prior/(1-Prior)
    datasim2$PostOdds <- datasim2$PriorOdds * datasim2$Unrelated
    
    LRgenH1 <- datasim$Related
    datasim1 <- data.frame(Related = LRgenH1)
    datasim1$PriorOdds <- Prior/(1-Prior)
    datasim1$PostOdds <- datasim1$PriorOdds * datasim1$Related
  }
  
  data <- data.frame(
    Unrelated = datasim2$PostOddsGen,
    Related = datasim1$PostOddsGen
  )
  
  return(data)
}
