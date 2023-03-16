#' Missing person based conditioned probability
#'
#' @param MPs Missing person sex 
#' @param MPc Missing person hair color
#' @param eps sex epsilon
#' @param epa age epsilon - Age is not specified in this first version, because it asumes uniformity.
#' @param epc color model
#' @export
#' @return A value of Likelihood ratio based on preliminary investigation data. In this case, sex.
#' @examples
#' CPT_MP()


CPT_MP <- function(MPs = "F", MPc = 1, eps = 0.05, epa = 0.05, epc = Cmodel()) {
  jointname <- c("F-T1", "F-T0", "M-T1", "M-T0")
  jointprob <- c((1-eps)*(1-epa), (1-eps)*epa, eps*(1-epa), eps*epa)
  names(jointprob) <- jointname

  Col <- c(1,2,3,4,5)
  probC = epc[MPc,]
  names(probC) <- Col

  CPTmp <- outer(jointprob,probC)
  return(CPTmp)
}

