#' Likelihood ratio for age variable
#'
#' @param eps epsilon paramenter.
#' @param erRs error rate in the database.
#' @param nsims number of simulations performed.
#' @param Ps Sex probabilities in the population. 
#' @param H hypothesis tested, H1: UHR is MP, H2: UHR is no MP
#' @param MPs MP sex
#' @param LR compute LR values
#' @export
#' @return A value of Likelihood ratio based on preliminary investigation data. In this case, sex.
#' @examples
#' LRsex() 


LRsex <- function(MPs = "F", eps = 0.05, erRs = eps, nsims = 1000, Ps = c(0.5,0.5), H = 1,  LR = FALSE) {
#Seteo general
sims <- list()  
S = c("F", "M")

#Interpretar sexo de MP
MPss <- which(S == MPs)
NoMPsn <- S[-MPss]
noMPs <- which(S == NoMPsn)

#Simular cuando H1 es cierta
if(H == 1) {

  x = c(S[MPss], S[noMPs])
  sims=as.data.frame(sample(x, size = nsims, prob = c(1-erRs, erRs), replace = TRUE))
names(sims) <- "Sexo"}

#Simular cuando H2 es cierta
else if (H == 2) {

  x = c(S[MPss], S[noMPs])
  sims=as.data.frame(sample(x, size = nsims, prob = c(Ps[1], Ps[2]), replace = TRUE))
  names(sims) <- "Sexo"
}

if (LR == TRUE) { 
  LRmatch = (1 - eps)/Ps[MPss]
  LRnomatch =  eps/Ps[noMPs]
  LRs <- lapply(sims, function(x) ifelse(x==MPs,  LRmatch, LRnomatch))
  sims <- cbind(sims, LRs)
  names(sims) <- c("Sexo", "LRs")
  return(sims)}

#Devolver simulaciones crudas o proporciones
return(sims)
}
