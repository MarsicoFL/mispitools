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
#' @import tidyverse
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

postSim <- function(datasim, Prior = 0.01, PriorModel = c("prelim","uniform")[1] ,eps = 0.05, erRs = 0.01, 
                    epc = Cmodel(), erRc = Cmodel(), MPc = 1, 
                    epa = 0.05, erRa = 0.01, MPa = 10, MPr = 2) {
  
  nsims = nrow(datasim)
  
  LRa <- LRc <- LRs <- LRtot <- PostOdds <- PriorOdds <- NULL

  if (PriorModel == "prelim") {
    LRgenH2 <- datasim$Unrelated
    #simulaciones para H2
    LRs2 <- LRsex(LR = TRUE, H = 2, nsims = nsims, eps = eps, erRs = erRs )
    LRc2 <- LRcol( LR = TRUE, H = 2,  nsims = nsims, epc = epc, erRc = erRc, MPc = MPc)
    LRa2 <- LRage(LR = TRUE, H = 2,  nsims = nsims, epa = epa, erRa = erRa, MPa = MPa, MPr = MPr )
    datasim2 <- as.data.frame(cbind(LRs2$LRs, LRc2$LRc, LRa2$LRa))
    names(datasim2) <- c("LRs", "LRc", "LRa")
    datasim2 <- mutate(datasim2, LRtot = LRa*LRc*LRs)
    datasim2 <- mutate(datasim2, Label = "H2")
    datasim2 <- cbind(datasim2, LRgenH2)
    datasim2 <- mutate(datasim2, PriorOdds = Prior/(1-Prior))
    datasim2 <- mutate(datasim2, PostOdds = PriorOdds*LRtot)
    datasim2 <- mutate(datasim2, PostOddsGen = PostOdds*LRgenH2)
    
    #simulaciones para H1
    LRgenH1 <- datasim$Related
    
    LRs1 <- LRsex(LR = TRUE, H = 1, nsims = nsims, eps = eps, erRs = erRs )
    LRc1 <- LRcol(LR = TRUE, H = 1,  nsims = nsims, epc = epc, erRc = erRc, MPc = MPc)
    LRa1 <- LRage(LR = TRUE, H = 1,  nsims = nsims, epa = epa, erRa = erRa, MPa = MPa, MPr = MPr )
    datasim1 <- as.data.frame(cbind(LRs1$LRs, LRc1$LRc, LRa1$LRa))
    names(datasim1) <- c("LRs", "LRc", "LRa")
    datasim1 <- mutate(datasim1, LRtot = LRa*LRc*LRs)
    datasim1 <- mutate(datasim1, Label = "H1")
    datasim1 <- cbind(datasim1, LRgenH1)
    datasim1 <- mutate(datasim1, PriorOdds = Prior/(1-Prior))
    datasim1 <- mutate(datasim1, PostOdds = PriorOdds*LRtot)
    datasim1 <- mutate(datasim1, PostOddsGen = PostOdds*LRgenH1)} else {
      LRgenH2 <- datasim$Unrelated
      datasim2 <- mutate(LRgenH2, PriorOdds = Prior/(1-Prior))
      datasim2 <- mutate(datasim2, PostOdds = PriorOdds*Unrelated)
      
      LRgenH1 <- datasim$Related
      datasim1 <- mutate(LRgenH1, PriorOdds = Prior/(1-Prior))
      datasim1 <- mutate(datasim1, PostOdds = PriorOdds*Related)
    }
  
  Related <- datasim1$PostOddsGen
  Unrelated <- datasim2$PostOddsGen
  data <- cbind(Unrelated, Related)
  
  return(as.data.frame(data))}
