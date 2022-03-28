#' Simulate likelihoods ratio (LRs) based on preliminary investigation data: a function for obtaining expected LRs under relatedness and unrelatedness kinship hypothesis.
#'
#' @param vartype Indicates type of preliminary investigation variable. Options are: sex, region, age, birthDate and height.
#' @param ErrorRate Error rate for sex, region, age and Height LR calculations. 
#' @param alphaBdate Vector containing alpha parameters for Dirichlet distribution.
#' @param numReg Number of regions present in the case.
#' @param seed Seed for simulations.
#' @param numsims Number of simulations performed.
#'
#' @return An object of class data.frame with LRs obtained for both hypothesis, Unrelated where POI/UHR is not MP or Related where POI/UHR is MP.
#' @export
#' @import DirichletReg
#' @import purrr
#' @examples
#' library(mispitools) 
#' simLRprelim("sex")




simLRprelim = function(vartype, numsims = 1000, seed = 123, ErrorRate = 0.05, alphaBdate = c(1, 4, 60, 11, 6, 4, 4), numReg = 6) {

set.seed(seed)

if (vartype == "sex") {
	sexLRvalues <- c((1-ErrorRate)/0.5, ErrorRate/0.5)
	a <- sample(sexLRvalues, numsims, replace = TRUE, prob = c(1-ErrorRate, ErrorRate))
        b <- sample(sexLRvalues, numsims, replace = TRUE, prob = c(0.5, 0.5))}
else if (vartype == "region") {
	H2 = 1/numReg
	RegLRvalues <- c((1-ErrorRate)/(1/H2), ErrorRate/H2)
        a <- sample(RegLRvalues, numsims, replace = TRUE, prob = c(1-ErrorRate, ErrorRate))
        b <- sample(RegLRvalues, numsims, replace = TRUE)}
else if (vartype == "age") {
	AgeLRvalues <- c((1-ErrorRate)/0.5, ErrorRate/0.5)
        a <- sample(AgeLRvalues, numsims, replace = TRUE, prob = c(1-ErrorRate, ErrorRate))
        b <- sample(AgeLRvalues, numsims, replace = TRUE, prob = c(0.5, 0.5))}
else if (vartype == "height") {
	heightLRvalues <- c((1-ErrorRate)/0.5, ErrorRate/0.5)
        a <- sample(heightLRvalues, numsims, replace = TRUE, prob = c(1-ErrorRate, ErrorRate))
        b <- sample(heightLRvalues, numsims, replace = TRUE, prob = c(0.5, 0.5))}
else if (vartype == "birthdate") {
	x = DirichletReg::rdirichlet(500, alphaBdate)
	temp <- dim(x); n <- temp[1]; m <- temp[2]
	lpb <- apply(log(x),2,mean)
	mom <- apply(x,2,mean)*(mean(x[,1])-mean(x[,1]^2))/(mean(x[,1]^2) - ((mean(x[,1]))^2))
	fit2= as.list(mom/sum(mom))
	bdateLRvalues <- purrr::map2(fit2, (1/length(fit2)), ~ .x / .y)
	a <- sample(bdateLRvalues, numsims, replace = TRUE, prob = fit2)
        b <- sample(bdateLRvalues, numsims, replace = TRUE)}

LRsimulated <- base::cbind(a,b)
base::colnames(LRsimulated) <- c("Unrelated", "Related")
base::structure(base::as.data.frame(LRsimulated))
}
