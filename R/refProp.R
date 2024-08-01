#' Generate a dataframe with hair colour, skin colour, eye colour and their specific combination frequencies
#'
#' This function creates a dataframe that lists every unique combination of hair colour,
#' skin colour, and eye colour in the provided dataset, along with the proportion of
#' occurrences of each combination.
#'
#' @param data A data.frame containing the characteristics of individuals.
#' @return A data.frame with columns for hair_colour, skin_colour, eye_colour, and f_h_s_y.
#' 
#' @examples
#' data <- simRef(1000)
#' refProp(data)
#' @export
refProp <- function(data) {
  temp <- as.data.frame(table(data$hair_colour, data$skin_colour, data$eye_colour))
  names(temp) <- c("hair_colour", "skin_colour", "eye_colour", "count")
  
  temp$f_h_s_y <- temp$count / nrow(data)
  
  result <- temp[, c("hair_colour", "skin_colour", "eye_colour", "f_h_s_y")]
  
  return(result)
}
