#' simLR2dataframe: A function for extracting LR distributions in a dataframe from simLRgen() output.
#'
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#'
#' @export
#' @return A dataframe with LR values obtained from simulations.

simLR2dataframe = function(datasim) {

  unrelated_values <- vector()
  related_values <- vector()

  list_length <- length(datasim[["Unrelated"]])

  for (i in 1:list_length) {
    unrelated_value <- datasim[["Unrelated"]][[i]][["LRtotal"]][["H1:H2"]]
    related_value <- datasim[["Related"]][[i]][["LRtotal"]][["H1:H2"]]

    unrelated_values <- c(unrelated_values, unrelated_value)
    related_values <- c(related_values, related_value)
  }

  results_df <- data.frame(Unrelated = unrelated_values, Related = related_values)
  datasim <- results_df
  }
