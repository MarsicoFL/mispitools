#' Simulate likelihoods ratio (LRs) based on preliminary investigation data: a function for obtaining expected LRs under relatedness and unrelatedness kinship hypothesis.
#'
#' @param vartype Indicates type of preliminary investigation variable. Options are: sex, region, age, birthDate and height.
#' @param ErrorRate Error rate for sex, region, age and Height LR calculations. 
#' @param alphaBdate Vector containing alpha parameters for Dirichlet distribution. Usually they are the frequencies of the solved cases in each category.
#' @param numReg Number of regions present in the case.
#' @param seed Seed for simulations.
#' @param numsims Number of simulations performed.
#' @param MP Introduce the preliminary data of the selected variable (vartype) of the MP. If it is null, open search is carried out. If it is not NULL, close search LR is computed. Variables values must be named as those presented in makePOIprelim function.
#' @param database It is used when the close search (MP not NULL), is carried out. It could be the output from makePOIprelim or a database with the same structure.
#' @param cuts Value of differences between DBD and ABD used for category definition. They must be the same as the ones selected for alphaBdate vector.
#' @param int Interval parameter, used for height and age vartypes. It defines the estimation range, for example, if MP age is 55, and int is 10, the estimated age range will be between 45 and 65.
#'
#' @return An object of class data.frame with LRs obtained for both hypothesis, Unrelated where POI/UHR is not MP or Related where POI/UHR is MP.
#' @export
#' @import DirichletReg
#' @import purrr
#' @examples
#' library(mispitools) 
#' simLRprelim("sex")




simLRprelim = function(vartype, numsims = 1000, seed = 123, int = 5, ErrorRate = 0.05, alphaBdate = c(1, 4, 60, 11, 6, 4, 4), numReg = 6, MP = NULL, database, cuts = c(-120, -30, 30, 120, 240, 360)) {

set.seed(seed)
if (is.null(MP)){
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
}

#Laburando el escenario cerrado
else if (!is.null(MP)){
if (vartype == "sex") {
	fs = sum(database$Sex == MP)/length(database$Sex)
        sexLRvalues <- c((1-ErrorRate)/fs, ErrorRate/(1-fs))
        a <- sample(sexLRvalues, numsims, replace = TRUE, prob = c(1-ErrorRate, ErrorRate))
        b <- sample(sexLRvalues, numsims, replace = TRUE, prob = c(fs, 1-fs))}
else if (vartype == "region") {
	fr = sum(database$Region == MP)/length(database$Region)
        RegLRvalues <- c((1-ErrorRate)/fr, ErrorRate/(1-fr))
        a <- sample(RegLRvalues, numsims, replace = TRUE, prob = c(1-ErrorRate, ErrorRate))
        b <- sample(RegLRvalues, numsims, replace = TRUE, prob = c(fr, 1-fr))}
else if (vartype == "age") {
	fa = sum(database$Age < (MP+int) & database$Age  > (MP-int))/length(database$Age)
	AgeLRvalues <- c((1-ErrorRate)/fa, ErrorRate/(1-fa))
        a <- sample(AgeLRvalues, numsims, replace = TRUE, prob = c(1-ErrorRate, ErrorRate))
        b <- sample(AgeLRvalues, numsims, replace = TRUE, prob = c(fa, 1-fa))}
else if (vartype == "height") {
        fh = sum(database$Height < (MP+int) & database$Height  > (MP-int))/length(database$Height)
	heightLRvalues <- c((1-ErrorRate)/fh, ErrorRate/(1-fh))
        a <- sample(heightLRvalues, numsims, replace = TRUE, prob = c(1-ErrorRate, ErrorRate))
        b <- sample(heightLRvalues, numsims, replace = TRUE, prob = c(fh, 1-fh))}
else if (vartype == "birthdate") {
	MP = as.Date(MP)
	ABD <- NULL
	DBD <- NULL
	Freq <- NULL
	PrelimData <- mutate(database, ABD = MP)
	PrelimData <- mutate(PrelimData, Dis = as.numeric(DBD - ABD))
	PrelimData <- as.data.frame(table(cut(PrelimData$Dis, breaks = c(-Inf,cuts, Inf))))
	PrelimData <- mutate(PrelimData, relfreq = Freq/sum(PrelimData$Freq))
	alpha2 <- as.vector(PrelimData$relfreq)
	x = DirichletReg::rdirichlet(500, alphaBdate)
	temp <- dim(x); n <- temp[1]; m <- temp[2]
	lpb <- apply(log(x),2,mean)
	mom <- apply(x,2,mean)*(mean(x[,1])-mean(x[,1]^2))/(mean(x[,1]^2) - ((mean(x[,1]))^2))
	fit2= as.list(mom/sum(mom))
	bdateLRvalues <- purrr::map2(fit2, alpha2, ~ .x / .y)
	a <- sample(bdateLRvalues, numsims, replace = TRUE, prob = fit2)
	b <- sample(bdateLRvalues, numsims, replace = TRUE, prob = alpha2)}}
LRsimulated <- base::cbind(b,a)
base::colnames(LRsimulated) <- c("Unrelated", "Related")
base::structure(base::as.data.frame(LRsimulated))
}
