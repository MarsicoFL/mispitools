#' Multi-dataset Kullback-Leibler Divergence Calculation
#'
#' @description
#' This function calculates the Kullback-Leibler divergence for all pairs of provided datasets,
#' considering allele frequencies. It normalizes data, adjusts zero frequencies,
#' and computes KL divergence in both directions for each pair.
#'
#' @param datasets List of dataframes, each containing allele frequencies for different populations.
#' @param minFreq Minimum frequency to be considered for unobserved or poorly observed alleles.
#' @import dplyr
#' @export
#' @return A matrix containing the Kullback-Leibler divergence for each dataset pair.
#' @examples
#' kl_matrix <- multi_kl_divergence(list(Argentina, BosniaHerz, Europe))
multi_kl_divergence <- function(datasets, minFreq = 1e-10) {
  num_datasets <- length(datasets)
  kl_matrix <- matrix(0, nrow = num_datasets, ncol = num_datasets, dimnames = list(NULL, NULL))

  for (i in 1:num_datasets) {
    for (j in 1:num_datasets) {
      if (i != j) {
        data1 <- datasets[[i]]
        data2 <- datasets[[j]]
        common_markers <- intersect(names(data1), names(data2))

        if (!"Allele" %in% names(data1)) next  # Verificar que 'Allele' esté en ambos datasets

        common_cols1 <- c("Allele", common_markers[common_markers != "Allele"])
        common_cols2 <- c("Allele", common_markers[common_markers != "Allele"])

        merged <- full_join(data1[, common_cols1, drop = FALSE], data2[, common_cols2, drop = FALSE], by = "Allele", suffix = c(".1", ".2"))
        merged[is.na(merged)] <- 0
        merged[merged == 0] <- minFreq

        # Normalization
        merged <- merged %>%
          mutate(across(matches("\\.1$"), ~ ./sum(.))) %>%
          mutate(across(matches("\\.2$"), ~ ./sum(.)))

        # KL divergence calculation
        kl_div <- sum(merged[, grep("\\.1$", names(merged))] * log10(merged[, grep("\\.1$", names(merged))] / merged[, grep("\\.2$", names(merged))]))
        
        kl_matrix[i, j] <- kl_div
      }
    }
  }

  return(kl_matrix)
}
