#' Bidirectional Kullback-Leibler Divergence for Genetic Markers
#'
#' @description
#' Computes the Kullback-Leibler (KL) divergence between allele frequency
#' distributions of two populations. Calculates divergence in both directions
#' to assess asymmetric differences between populations.
#'
#' @param data1 A data frame with allele frequencies for the first population.
#'   First column should be "Allele", subsequent columns are marker frequencies.
#' @param data2 A data frame with allele frequencies for the second population.
#'   Same format as \code{data1}.
#' @param minFreq Numeric. Minimum frequency value to replace zeros (to avoid
#'   undefined logarithms). Default: 1e-10.
#'
#' @return A named list with two elements:
#'   \itemize{
#'     \item \code{"KL from data1 to data2"}: KL(P || Q) - divergence from
#'       population 1 to population 2
#'     \item \code{"KL from data2 to data1"}: KL(Q || P) - divergence from
#'       population 2 to population 1
#'   }
#'   Values are computed using log base 10.
#'
#' @details
#' The KL divergence measures how one probability distribution differs from
#' another. It is asymmetric: KL(P || Q) != KL(Q || P).
#'
#' Higher values indicate greater divergence between populations, which may
#' affect the reliability of LR calculations when using frequency data from
#' one population to analyze individuals from another.
#'
#' The function:
#' \enumerate{
#'   \item Finds markers common to both populations
#'   \item Merges allele frequency tables
#'   \item Replaces missing/zero frequencies with \code{minFreq}
#'   \item Normalizes frequencies to sum to 1
#'   \item Computes KL divergence in both directions
#' }
#'
#' @seealso
#' \code{\link{kl_multi}} for comparing multiple populations,
#' \code{\link{kl_pie}} for matrix-based KL divergence.
#'
#' @references
#' Kullback S, Leibler RA (1951). "On Information and Sufficiency."
#' \emph{The Annals of Mathematical Statistics}, 22(1), 79-86.
#'
#' @export
#' @import dplyr
#' @examples
#' # Compare Argentina and Bosnia-Herzegovina populations
#' result <- kl_bidirectional(Argentina, BosniaHerz)
#' print(result)
#'
#' # Compare Argentina and Europe
#' kl_bidirectional(Argentina, Europe)

kl_bidirectional <- function(data1, data2, minFreq = 1e-10) {

  common_markers <- intersect(names(data1), names(data2))

  # Ensure only unique columns are selected
  common_cols1 <- c("Allele", common_markers[!common_markers %in% "Allele"])
  common_cols2 <- c("Allele", common_markers[!common_markers %in% "Allele"])

  merged1 <- dplyr::full_join(data1[, common_cols1, drop = FALSE],
                              data2[, common_cols2, drop = FALSE],
                              by = "Allele", suffix = c(".1", ".2"))
  merged2 <- dplyr::full_join(data2[, common_cols2, drop = FALSE],
                              data1[, common_cols1, drop = FALSE],
                              by = "Allele", suffix = c(".2", ".1"))

  # Replace NA with 0
  merged1[is.na(merged1)] <- 0
  merged2[is.na(merged2)] <- 0

  # Replace 0 with minimum frequency
  merged1[merged1 == 0] <- minFreq
  merged2[merged2 == 0] <- minFreq

  # Normalize
  merged1 <- merged1 %>%
    dplyr::mutate(dplyr::across(dplyr::matches("\\.1$"), ~ . / sum(.))) %>%
    dplyr::mutate(dplyr::across(dplyr::matches("\\.2$"), ~ . / sum(.)))
  merged2 <- merged2 %>%
    dplyr::mutate(dplyr::across(dplyr::matches("\\.2$"), ~ . / sum(.))) %>%
    dplyr::mutate(dplyr::across(dplyr::matches("\\.1$"), ~ . / sum(.)))

  # Calculate KL divergence
  kl1 <- sum(merged1[, grep("\\.1$", names(merged1))] *
               log10(merged1[, grep("\\.1$", names(merged1))] /
                       merged1[, grep("\\.2$", names(merged1))]))
  kl2 <- sum(merged2[, grep("\\.2$", names(merged2))] *
               log10(merged2[, grep("\\.2$", names(merged2))] /
                       merged2[, grep("\\.1$", names(merged2))]))

  return(list("KL from data1 to data2" = kl1, "KL from data2 to data1" = kl2))
}
