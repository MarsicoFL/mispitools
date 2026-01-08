#' Compute Reference Population Proportions for Pigmentation Traits
#'
#' @description
#' Computes the frequency of each unique combination of hair color, skin color,
#' and eye color in the reference population. These proportions serve as the
#' denominator (H2 probabilities) in LR calculations for pigmentation traits.
#'
#' @param data A data.frame with columns \code{hair_colour}, \code{skin_colour},
#'   and \code{eye_colour}, typically output from \code{\link{sim_reference_pop}}.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{hair_colour}: Hair color category
#'     \item \code{skin_colour}: Skin color category
#'     \item \code{eye_colour}: Eye color category
#'     \item \code{f_h_s_y}: Population frequency of this combination
#'   }
#'
#' @seealso
#' \code{\link{sim_reference_pop}} for generating the input data,
#' \code{\link{compute_conditioned_prop}} for conditioned proportions,
#' \code{\link{lr_compute_pigmentation}} for computing LRs.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Generate reference population
#' pop_data <- sim_reference_pop(n = 500, seed = 123)
#'
#' # Compute proportions
#' ref_prop <- compute_reference_prop(pop_data)
#' head(ref_prop)
#'
#' # Most common combinations
#' ref_prop[order(-ref_prop$f_h_s_y), ][1:5, ]

compute_reference_prop <- function(data) {

  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data.frame (typically from sim_reference_pop)")
  }

  required_cols <- c("hair_colour", "skin_colour", "eye_colour")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (nrow(data) == 0) {
    stop("data has no rows. Cannot compute proportions from empty data.")
  }

  temp <- as.data.frame(table(data$hair_colour, data$skin_colour, data$eye_colour))
  names(temp) <- c("hair_colour", "skin_colour", "eye_colour", "count")

  temp$f_h_s_y <- temp$count / nrow(data)

  result <- temp[, c("hair_colour", "skin_colour", "eye_colour", "f_h_s_y")]

  return(result)
}
