#' Likelihood ratio for age variable
#'
#' @param MPa Missing person age
#' @param MPr Missing person age range.
#' @param UHRr Unidentified person range
#' @param gam Simulation parameter for UHR ages. 
#' @param nsims number of simulations.
#' @param erRa error rate in the database.
#' @param epa epsilon age
#' @param H hipothesis tested, H1: UHR is MP, H2: UHR is not MP.
#' @param modelA reference database probabilities, uniform assumes equally probable ages. Custom needs a vector with ages frequencies.
#' @param LR compute LR values
#' @export
#' @return A value of Likelihood ratio based on preliminary investigation data. In this case, Age.

LRage <- function(MPa = 40, MPr = 6, UHRr = 1, gam = 0.07, nsims = 1000, epa = 0.05, erRa = epa, H=1, modelA = c("uniform", "custom")[1], LR = FALSE) {

sims = list()  
Age = seq(1:80)
MPmin = MPa - MPr
MPmax = MPa + MPr

if (modelA == "uniform"){

T1p <- (MPmax-MPmin)/length(Age)  # Para una uniforme
T0p <-  1-T1p
LR1 <- (1 - epa)/T1p
LR0 <- epa/T0p}

T1a <- Age[Age < MPmax & Age > MPmin]
T0a <- Age[-T1a]

if(H == 1) {
  group= unlist(sample(c("T1", "T0"), size = nsims, prob = c(1 - erRa,  erRa), replace = TRUE))
  ages = unlist(lapply(group, function(x) ifelse(x=="T1",  sample(T1a, 1), sample(T0a, 1))))
  
  sims = as.data.frame(cbind(group,(ages)))
  names(sims) <- c("group", "age")

  sims <- mutate(sims, UHRmin = as.numeric(ages) - gam*as.numeric(ages) - UHRr) 
  sims <- mutate(sims, UHRmax = as.numeric(ages) + gam*as.numeric(ages) + UHRr)# queda el gamma asociado a la edad, modificar en futuras versiones?
  }

else if (H == 2) {
  ages = unlist(sample(Age, nsims, replace = TRUE))
  group= unlist(lapply(ages, function(x) ifelse(x > MPmin & x < MPmax,  "T1", "T0")))
  
  sims = as.data.frame(cbind(group,ages))
  names(sims) <- c("group", "age")
  sims <- mutate(sims, UHRmin = as.numeric(ages) - gam*as.numeric(ages)- UHRr) 
  sims <- mutate(sims, UHRmax = as.numeric(ages) + gam*as.numeric(ages)+ UHRr)# queda el gamma asociado a la edad, modificar en futuras versiones?
}

if (modelA == "custom"){
  T1p <- length(subset(sims$group, sims$group == "T1"))/length(sims$group)
  T0p <- 1 - T1p
  LR1 <- (1 - epa)/T1p
  LR0 <- epa/T0p}

if (LR == TRUE) { 
  sims <- mutate(sims, ifelse(group == "T1", LRa <- LR1, LRa <- LR0))
  names(sims) <- c("group", "Age", "UHRmin", "UHRmax", "LRa")
  return(sims)}

sims2 <- select(sims, group)
return(as.data.frame(sims))}

