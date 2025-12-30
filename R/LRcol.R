#' Likelihood ratio for hair colour variable
#'
#' Simulates a distribution of observed hair colours under either H1 (UHR is MP)
#' or H2 (UHR is not MP), and optionally returns LR values per simulated observation.
#'
#' @param epc Error matrix for hair colour observations (rows = true MP colour).
#' @param erRc Error matrix used for simulation under H1; defaults to `epc`.
#' @param nsims Number of simulations performed.
#' @param Pc Population hair colour probabilities.
#' @param H Hypothesis tested: 1 for H1 (UHR is MP), 2 for H2 (UHR is not MP).
#' @param Qprop Query colour tested.
#' @param MPc MP hair colour.
#' @param seed For reproducible simulations.
#' @param LR If `TRUE`, return LR values alongside simulated colours.
#' @export
#' @return A data.frame of simulated hair colours, and if `LR = TRUE`, a `LRc` column.
#' @examples
#' LRcol()
#' LRcol(MPc = 2, H = 2, nsims = 100, LR = TRUE)


LRcol <- function(MPc = 1, epc = Cmodel(),  erRc = epc, nsims = 1000, Pc = c(0.3,0.2, 0.25, 0.15,0.1), H=1,  Qprop = MPc, LR = FALSE, seed = 1234) {

sims <- list()  
Col <- c(1,2,3,4,5)

set.seed(seed)
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
