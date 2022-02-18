#' Function for getting STR allele frequencies from different world populations.
#' 
#' @param region select the place of the allele frequency database. Possible values are listed: "Argentina", "Asia", "Europe", "USA", "Austria", "BosniaHerz", "China" and "Japan".
#' @source https://doi.org/10.1016/j.fsigss.2009.08.178; https://doi.org/10.1016/j.fsigen.2016.06.008; https://doi.org/10.1016/j.fsigen.2018.07.013.
#' @return An allele frequency database adapted compatible with pedtools format. 
#' @export
#' @import tidyverse

getfreqs <- function(region){
    Freqs <- as.list(region)
    for(i in 2:length(Freqs)){
      names(Freqs[[i]]) <- Freqs[[1]]
    }
    Freqs$Allele <- NULL
    
    return(Freqs)
  
  }
