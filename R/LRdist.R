
#' Likelihood ratio distribution: a function for plotting expected log10(LR) distributions under relatedness and unrelatedness.
#'
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#' @param type Select between a density plot (type = 1, default) or a violin plot (type = 2).
#'
#' @export
#' @return A plot showing likelihood ratio distributions under relatedness and unrelatedness hypothesis.
#' @examples
#' library(forrel)
#' x = linearPed(2)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' datasim = simLRgen(x, missing = 5, 10, 123)
#' LRdist(datasim)
#' @import plotly
#' @import dplyr
#' @import highcharter
#' @import tidyr


LRdist = function(datasim, type = 1) {

TPED = log10(datasim$Related)
RPED = log10(datasim$Unrelated)

if(type == 1) {
	hc <- hchart(
  	stats::density(TPED), type = "area", 
  	color = "steelblue", name = "Related"
  	) %>%
  	hc_add_series(
    	stats::density(RPED), type = "area",
    	color = "#B71C1C", 
    	name = "Unrelated"
    	) %>%
    	hc_title(text = "LR distributions") %>%
    	hc_xAxis(title = list(text = "Log10(LR)")) %>%
    	hc_yAxis(title = list(text = "Density"))
	}
else if(type == 2) {
	datalog <- log10(datasim)
	datalog <- tidyr::gather(datalog)
	colnames(datalog) <- c("tipo", "LR")
	hc <- datalog %>%
	  plot_ly(
    	  	x = ~tipo,
    		y = ~LR,
    		type = 'violin',
    		box = list(
      		visible = T
    		),
    	        split = ~tipo,
    		meanline = list(
      		visible = T
   	 		)
  			 )

       hc <- hc %>%
  	plotly::layout(
    	xaxis = list(
      	title = "hypothesis"
    	),
    	yaxis = list(
      	title = "Log10(LR)",
      	zeroline = F
    	       )
  		)
       }
hc
}
