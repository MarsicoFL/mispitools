#' Epsilon hair color matrix
#'
#' @param errorModel custom allows selecting a specfic epsilon for each MP-UHR pair, uniform use ep for all. 
#' @param ep epsilon
#' @param ep12 epsilon
#' @param ep13 epsilon
#' @param ep14 epsilon
#' @param ep15 epsilon
#' @param ep23 epsilon
#' @param ep24 epsilon
#' @param ep25 epsilon
#' @param ep34 epsilon
#' @param ep35 epsilon
#' @param ep45 epsilon
#' @export
#' @return A value of Likelihood ratio based on preliminary investigation data. In this case, sex.
#' @examples
#' Cmodel()


Cmodel <- function(errorModel = c("custom","uniform")[1], ep = 0.01, ep12 = 0.01, ep13 = 0.005, ep14 = 0.01, ep15 = 0.003, ep23 =0.01, ep24 = 0.003, ep25 =0.01, ep34 = 0.003, ep35 = 0.003, ep45 = 0.01){
if(errorModel == "uniform") {
  ep12 <- ep13 <- ep14 <- ep15 <- ep23 <- ep25 <- ep25 <- ep34 <- ep35 <- ep45 <- ep
}
l1 = 1/(1+ep12+ep13+ep14+ep15)
l2 = 1/(1+ep12+ep23+ep24+ep25)
l3 = 1/(1+ep13+ep23+ep34+ep35)
l4 = 1/(1+ep14+ep24+ep34+ep45)
l5 = 1/(1+ep15+ep25+ep35+ep45)

errorMat <- rbind(c(l1, l1*ep12, l1*ep13, l1*ep14, l1*ep15), c(l2*ep12,l2,l2*ep13,l2*ep14,l2*ep15), c(l3*ep13, l3*ep23, l3, l3*ep34, l3*ep35), c(l4*ep14, l4*ep24, l4*ep34, l4, l4*ep45), c(l5*ep15, l5*ep25, l5*ep35, l5*ep45, l5))
return(errorMat)} 
