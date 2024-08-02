#' Compute Likelihood Ratios based con color characteristics
#'
#' This function calculates the Likelihood Ratios (LRs) for each combination of hair colour,
#' skin colour, and eye colour between two datasets. It assumes one dataset (`conditioned`)
#' contains numerators and the other (`unconditioned`) contains denominators.
#'
#' @param conditioned A dataframe with at least the columns 'hair_colour', 'skin_colour',
#' 'eye_colour', and 'numerators'.
#' @param unconditioned A dataframe with at least the columns 'hair_colour', 'skin_colour',
#' 'eye_colour', and 'f_h_s_y'.
#' @return A dataframe with the merged data and computed LRs.
#' @export
#' @examples
#' data <- simRef()
#' conditioned <- conditionedProp(data, 1, 1, 1, 0.01, 0.01, 0.01) 
#' unconditioned <- refProp(data)
#' compute_LRs_colors(conditioned, unconditioned)
compute_LRs_colors <- function(conditioned, unconditioned) {
    merged_data <- merge(unconditioned, conditioned, by = c("hair_colour", "skin_colour", "eye_colour"))
    
    merged_data$LR <- merged_data$numerators / merged_data$f_h_s_y
    
    return(merged_data)
}
