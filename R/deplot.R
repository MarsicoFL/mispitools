
#' Decision making plot: a function for plotting false positive and false negative rates for each LR threshold.
#'
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#'
#' @export
#'
#' @examples
#' #library(forrel)
#' #x = linearPed(2)
#' #x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' #x = profileSim(x, N = 1, ids = 2)[[1]]
#' #datasim = makeLRsims(x, missing = 5, 1000)
#' #deplot(datasim)
#' @import plotly


deplot = function(datasim) {

TPED = datasim$Related# We define the obtained values for H1
RPED = datasim$Unrelated# We define the obtained values for H2
nsimul = nrow(datasim) # Number of simulations performed

## False positives and false negatives calculations
ValoresLR = seq(1, nsimul, length.out= nsimul)
FPs = 0
FNs = 0
for(i in 1:nsimul) { FPs[i] = sum(RPED > ValoresLR[i]); #False positives.
FNs[i] = sum(TPED < ValoresLR[i])} #False negatives.

Datos = data.frame(x = ValoresLR, y= FPs/nsimul, z= FNs/nsimul, w=FPs)
p <- plotly::plot_ly(Datos, x = Datos$z, y = Datos$y,
  # Hover text:
  text = ~paste("LR threshold: ", Datos$x,
                "<br>expected false positives:", Datos$y*nsimul),
  type   = 'scatter', 
  mode   = 'markers',
  color = Datos$z, size = 1) %>%
  layout(xaxis = list(autotypenumbers = 'strict', title = 'False negative rate'),
         yaxis = list(title = 'False positive rate'))
p

}

