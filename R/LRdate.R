#' Likelihood ratio for birth date in missing person searches
#'
#' @param ABD Actual birth date of the missing person.
#' @param DBD Declared birth date of the person of interest.
#' @param alpha A vector containing the alpha values for the dirichlet. It should contain the number of categories of differences between DBD and ABD.
#' @param cuts Value of differences between DBD and ABD used for category definition.
#' @param type Type of scenario, type 1 is an "open search", where it is unknown if the missing person is in the database. Type 2 refers to a scenario where the missing person is in the database.
#' @param PrelimData Used when type = 2, is the dataframe with the DBD of the persons of interest in the database.
#' @param draw Number of simulations for Dirichlet distribution computation.
#' @param seed Seed for simulations.
#' @import DirichletReg
#' @export
#' @return A value of Likelihood ratio based on preliminary investigation data. In this case, birth date.
#' @examples
#' library(DirichletReg)
#' LRdate(ABD = "1976-05-31", DBD = "1976-07-15", 
#'	PrelimData, alpha = c(1, 4, 60, 11, 6, 4, 4), 
#'	cuts = c(-120, -30, 30, 120, 240, 360), 
#'	type = 1, seed = 123)



LRdate = function(ABD = "1976-05-31", DBD = "1976-07-15", PrelimData, alpha = c(1, 4, 60, 11, 6, 4, 4), cuts = c(-120, -30, 30, 120, 240, 360), draw = 500, type = 1, seed = 123) {
seed = seed
x = DirichletReg::rdirichlet(draw, alpha)
temp <- dim(x); n <- temp[1]; m <- temp[2]
lpb <- apply(log(x),2,mean)
mom <- apply(x,2,mean)*(mean(x[,1])-mean(x[,1]^2))/(mean(x[,1]^2) - ((mean(x[,1]))^2))
fit2= as.list(mom/sum(mom))


DBD = as.Date(DBD)
ABD = as.Date(ABD)

Dis0 = julian(DBD,ABD)
Dist = Dis0[1]
i = 0   
w = length(alpha) - 1
if(Dist <= cuts[1]) {H1 = fit2[1]} 
if(Dist >= cuts[w]) {H1 = fit2[w]}
if(Dist >= cuts[1] & Dist <= cuts[w]) {
for (i in 1:w) {
    up = i + 1
    if (Dist > cuts[i]  & Dist < cuts[up]){H1 = fit2[i]}}}
                
if(type == 1) {
  H2 = 1/length(alpha)
  LR = as.numeric(H1)/H2
  }

if(type == 2) {
  PrelimData <- mutate(PrelimData, Dis = julian(DBD,ABD))
  PrelimData <- as.data.frame(table(cut(PrelimData$Dis, breaks = c(-Inf,cuts, Inf))))
  alpha2 <- as.vector(PrelimData$Freq)
  x2 = DirichletReg::rdirichlet(draw, alpha2)
  temp2 <- dim(x2); n <- temp2[1]; m2 <- temp2[2]
  lpb2 <- apply(log(x2),2,mean)
  mom2 <- apply(x2,2,mean)*(mean(x2[,1])-mean(x2[,1]^2))/(mean(x2[,1]^2) - ((mean(x2[,1]))^2))
  fit3= as.list(mom2/sum(mom2))
  if(Dist <= cuts[1]){H2 = fit3[1]}
  if(Dist >= cuts[w]){H2 = fit3[w]}
  if(Dist >= cuts[1] & Dist <= cuts[w]) {
  for (i in 1:w) {
  up = i + 1
  if (Dist > cuts[i]  & Dist < cuts[up]){H2 = fit3[i]}}}

  LR = as.numeric(H1)/as.numeric(H2)
}
print(LR)
}
