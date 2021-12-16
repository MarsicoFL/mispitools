#' Likelihood ratio for preliminary investigation data: a function for computing likelihood ratio based on preliminary investigation data.
#'
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#' @param ABD Actual birth date of the missing person.
#' @param DBD Declared birth date of the person of interest.
#' @param alpha A vector containing the alpha values for the dirichlet. It should contain the number of categories of differences between DBD and ABD.
#' @param cuts Value of differences between DBD and ABD used for category definition.
#' @param type Type of scenario, type 1 is an "open search", where it is unknown if the missing person is in the database. Type 2 refers to a scenario where the missing person is in the database.
#' @param PrelimData Used when type = 2, is the dataframe with the DBD of the persons of interest in the database.
#' @param Draw Number of simulations for Dirichlet distribution computation.
#' @import dirichlet
#' @export
#' @return A value of Likelihood ratio based on preliminary investigation data. In this case, birth date.
#' @examples
#' library(dirichlet)
#' LRprelim(ABD = "1976-05-31", DBD = "1976-07-15", PrelimData, alpha = c(1, 4, 60, 11, 6, 4, 4), cuts = c(-120, -30, 30, 120, 240, 360), type = 1)



LRprelim = function(ABD = "1976-05-31", DBD = "1976-07-15", PrelimData, alpha = c(1, 4, 60, 11, 6, 4, 4), cuts = c(-120, -30, 30, 120, 240, 360), draw = 500, type = 1) {

x = dirichlet::rdirichlet(draw, alpha)
fit = data.frame(dirichlet::fit.genDirichlet(x, "ml"))
fit2 = as.list(fit$p)
DBD = as.Date(DBD)
ABD = as.Date(ABD)

Dis0 = julian(DBD,ABD)
Dis = Dis0[1]
i = 0   
w = length(alpha) - 1
if(Dis <= cuts[1]) {H1 = fit2[1]} 
if(Dis >= cuts[w]) {H1 = fit2[w]}
if(Dis >= cuts[1] & Dis <= cuts[w]) {
for (i in 1:w) {
    up = i + 1
    if (Dis > cuts[i]  & Dis < cuts[up]){H1 = fit2[i]}}}
                
if(type == 1) {
  H2 = 1/length(alpha)
  LR = as.numeric(H1)/H2
  }

if(type == 2) {
  PrelimData <- mutate(PrelimData, Dis = julian(DBD,ABD))
  PrelimData <- as.data.frame(table(cut(PrelimData$Dis, breaks = c(-Inf,cuts, Inf))))
  alpha2 <- as.vector(PrelimData$Freq)
  x2 = dirichlet::rdirichlet(draw, alpha2)
  fit3 = data.frame(fit.genDirichlet(x2, "ml"))
  
  if(Dis <= cuts[1]){H2 = fit3[1]}
  if(Dis >= cuts[w]){H2 = fit3[w]}
  if(Dis >= cuts[1] & Dis <= fit3[w]) {
  for (i in 1:w) {
  up = i + 1
  if (Dis > cuts[i]  & Dis < cuts[up]){H2 = fit3[i]}}}

  LR = as.numeric(H1)/as.numeric(H2)
}
print(LR)
}
