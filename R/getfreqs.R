#' Function for getting STR allele frequencies from different world populations.
#' 
#' @param region select the place of the allele frequency database. Possible values are listed: "Argentina", "Asia", "Europe", "USA", "Austria", "BosniaHerz", "China" and "Japan". Sources: https://doi.org/10.1016/j.fsigss.2009.08.178; https://doi.org/10.1016/j.fsigen.2016.06.008; https://doi.org/10.1016/j.fsigen.2018.07.013.
#'
#' @return An allele frequency database. 
#' @export
#' @import tidyverse
#' @examples
#' getfreqs("Argentina")

getfreqs <- function(region){
#Get alelos de leapdna:
indice.es.alelos <- data.frame(country = c("Argentina",
                                              "Asia",
                                              "Austria",
					      "BosniaHerz",
					      "China",
					      "Europe",
					      "Japan",
					      "USA"
                                              ),
                                    
                                    urls = c('data/Argentina.rda',
                                             "data/Asia.rda",
				    	     "data/Austria.rda",
				    	     "data/BosniaHerz.rda",
				             "data/China.rda",
				             "data/Europe.rda",
				             "data/Japan.rda",
				             "data/USA.rda")
                                    ) %>% 
  arrange(country)
  
  indice = indice.es.alelos
  leap.url <- indice[indice$country == region, "urls"]
  
    prueba <- get(load(file = leap.url))
    Freqs <- as.list(prueba)
    for(i in 2:length(Freqs)){
      names(Freqs[[i]]) <- Freqs[[1]]
    }
    Freqs$Allele <- NULL
    
    return(Freqs)
  
  }
