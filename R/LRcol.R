#' Likelihood ratio for age variable
#'
#' @param epc epsilon paramenter.
#' @param erRc error rate in the database.
#' @param nsims number of simulations performed.
#' @param Pc hair color probabilities. 
#' @param H hypothesis tested, H1: UHR is MP, H2: UHR is no MP
#' @param Qprop Query color tested.
#' @param MPc MP hair color
#' @param LR compute LR values
#' @export
#' @return A value of Likelihood ratio based on preliminary investigation data. In this case, hair color.
#' @examples
#' LRcol() 


LRcol <- function(MPc = 1, epc = Cmodel(),  erRc = epc, nsims = 1000, Pc = c(0.3,0.2, 0.25, 0.15,0.1), H=1,  Qprop = MPc, LR = FALSE) {

sims <- list()  
Col <- c(1,2,3,4,5)


if(H == 1) {

  x = Col
  sims=as.data.frame(sample(x, size = nsims, prob = erRc[MPc,], replace = TRUE))
names(sims) <- "Col"}

else if (H == 2) {

  x = Col
  sims=as.data.frame(sample(x, size = nsims, prob = Pc, replace = TRUE))
  names(sims) <- "Col"
}

if (LR == TRUE) { 
  LRs <- lapply(sims, function(x) ifelse(x==MPc,  epc[MPc,x]/Pc[x], epc[MPc,x]/Pc[x]))
  sims <- cbind(sims, LRs)
  names(sims) <- c("Col", "LRc")
  return(sims)}
else {return(sims)
}}

