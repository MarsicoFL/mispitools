#' Compute likelihood ratios from pigmentation traits
#'
#' This function calculates likelihood ratios (LRs) for each combination of hair, skin, and eye
#' colour across two datasets. It assumes `conditioned` provides the numerators and
#' `unconditioned` provides the denominators for the same trait combinations.
#' Combinations not present in both inputs are dropped by the merge.
#'
#' @param conditioned A data.frame with at least the columns `hair_colour`, `skin_colour`,
#'   `eye_colour`, and `numerators`.
#' @param unconditioned A data.frame with at least the columns `hair_colour`, `skin_colour`,
#'   `eye_colour`, and `f_h_s_y`.
#' @return A data.frame with the merged data and the computed `LR` column.
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
