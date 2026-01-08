#' Multi-Population Kullback-Leibler Divergence Matrix
#'
#' @description
#' Computes pairwise Kullback-Leibler (KL) divergence between multiple
#' population allele frequency databases. Returns a matrix of divergence
#' values useful for assessing population differentiation.
#'
#' @param datasets A list of data frames, each containing allele frequencies
#'   for a different population. Each data frame should have "Allele" as
#'   the first column and marker frequencies in subsequent columns.
#' @param minFreq Numeric. Minimum frequency to replace zeros.
#'   Default: 1e-10.
#'
#' @return A square numeric matrix where element (i,j) contains the KL
#'   divergence from population i to population j. Diagonal elements are 0.
#'
#' @details
#' This function is useful for:
#' \itemize{
#'   \item Comparing multiple reference populations
#'   \item Selecting the most appropriate frequency database for a case
#'   \item Assessing potential bias from population mismatch
#' }
#'
#' Higher values indicate greater divergence between populations.
#'
#' @seealso
#' \code{\link{kl_bidirectional}} for pairwise comparison of two populations.
#'
#' @references
#' Kullback S, Leibler RA (1951). "On Information and Sufficiency."
#' \emph{The Annals of Mathematical Statistics}, 22(1), 79-86.
#'
#' @export
#' @import dplyr
#' @examples
#' # Compare three populations
#' kl_matrix <- kl_multi(list(Argentina, BosniaHerz, Europe))
#' print(kl_matrix)
#'
#' # Visualize as heatmap
#' # heatmap(kl_matrix, main = "KL Divergence Between Populations")

kl_multi <- function(datasets, minFreq = 1e-10) {

  num_datasets <- length(datasets)
  kl_matrix <- matrix(0, nrow = num_datasets, ncol = num_datasets,
                      dimnames = list(NULL, NULL))

  for (i in 1:num_datasets) {
    for (j in 1:num_datasets) {
      if (i != j) {
        data1 <- datasets[[i]]
        data2 <- datasets[[j]]
        common_markers <- intersect(names(data1), names(data2))

        # Check that 'Allele' is in both datasets
        if (!"Allele" %in% names(data1)) next

        common_cols1 <- c("Allele", common_markers[common_markers != "Allele"])
        common_cols2 <- c("Allele", common_markers[common_markers != "Allele"])

        merged <- dplyr::full_join(data1[, common_cols1, drop = FALSE],
                                   data2[, common_cols2, drop = FALSE],
                                   by = "Allele", suffix = c(".1", ".2"))
        merged[is.na(merged)] <- 0
        merged[merged == 0] <- minFreq

        # Normalization
        merged <- merged %>%
          dplyr::mutate(dplyr::across(dplyr::matches("\\.1$"), ~ . / sum(.))) %>%
          dplyr::mutate(dplyr::across(dplyr::matches("\\.2$"), ~ . / sum(.)))

        # KL divergence calculation
        kl_div <- sum(merged[, grep("\\.1$", names(merged))] *
                        log10(merged[, grep("\\.1$", names(merged))] /
                                merged[, grep("\\.2$", names(merged))]))

        kl_matrix[i, j] <- kl_div
      }
    }
  }

  return(kl_matrix)
}
