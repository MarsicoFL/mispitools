#' Population based conditioned probability
#'
#' @param MPa Missing person sex 
#' @param MPr Missing person hair color
#' @param propC sex epsilon
#' @param propS age epsilon - Age is not specified in this first version, because it asumes uniformity.
#' @export
#' @return A value of Likelihood ratio based on preliminary investigation data. In this case, sex.
#' @examples
#' CPT_POP()


CPT_POP <- function(propS = c(0.5,0.5), MPa = 40, MPr = 6, propC = c(0.3,0.2, 0.25, 0.15,0.1)) {
  Age <- seq(1:80)
  MPmin <- MPa - MPr
  MPmax <- MPa + MPr
  T1p <- (MPmax-MPmin)/length(Age)  # Para una uniforme
  T0p <-  1-T1p
  
  jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")
  jointprob <- c(propS[1]*T1p, propS[1]*T0p, propS[2]*T1p, propS[2]*T0p)
  names(jointprob) <- jointname
  
  CPTpop <- outer(jointprob,propC)
  return(CPTpop)
} 
