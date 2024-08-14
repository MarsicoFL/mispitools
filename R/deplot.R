
#' Decision making plot: a function for plotting false positive and false negative rates for each LR threshold.
#'
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#' @param LRmax Maximum LR value used as a threshold. 1000 setted by default.
#'
#' @export
#' @return A plot showing false positive and false negative rates for each likelihood ratio threshold.
#' @examples
#' library(forrel)
#' x = linearPed(2)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' datasim = simLRgen(x, missing = 5, 10, 123)
#' deplot(datasim)
#' @import ggplot2


deplot = function(datasim, LRmax = 1000) {
if (!is.data.frame(datasim)) {
   datasim <- simLR2dataframe(datasim)
 }

x <- y <- z <- NULL

TPED = datasim$Related
RPED = datasim$Unrelated
nsimul = nrow(datasim) 

ValoresLR = seq(1, LRmax, length.out= LRmax)
FPs = 0
FNs = 0
for(i in 1:nsimul) { FPs[i] = sum(RPED > ValoresLR[i]); 
FNs[i] = sum(TPED < ValoresLR[i])} 

Datos = data.frame(x = ValoresLR, y= FPs/nsimul, z= FNs/nsimul, w=FPs)
primer_punto <- 1
ultimo_punto <- nrow(Datos)

p <- ggplot2::ggplot(Datos, ggplot2::aes(x = z, y = y)) +
  ggplot2::geom_point(size = 2, color = "black") +  
  
  ggplot2::geom_text(data = Datos[c(primer_punto, ultimo_punto), ],
                     ggplot2::aes(label = x),
                     hjust = 1.1, vjust = -0.5, size = 3) +
  
  ggplot2::scale_x_continuous(name = "False negative rate") +
  ggplot2::scale_y_continuous(name = "False positive rate") +
  
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 12),
    axis.title.y = ggplot2::element_text(size = 12),
    axis.text.x = ggplot2::element_text(size = 10),
    axis.text.y = ggplot2::element_text(size = 10)
  )

p

}

