#' Kullback-Leibler Divergence Calculation for Genetic Markers
#'
#' @description
#' This function calculates the Kullback-Leibler divergence for shared genetic markers
#' between two populations, considering allele frequencies. It normalizes data, adjusts zero frequencies,
#' and calculates divergence in both directions.
#'
#' @param data1 DataFrame with allele frequencies for the first population.
#' @param data2 DataFrame with allele frequencies for the second population.
#' @param minFreq Minimum frequency to be considered for unobserved or poorly observed alleles.
#' @import dplyr
#' @export
#' @return A list containing the Kullback-Leibler divergence from data1 to data2 and vice versa.
#' @examples
#' bidirectionalKL(Argentina, BosniaHerz)

bidirectionalKL <- function(data1, data2, minFreq = 1e-10) {

  common_markers <- intersect(names(data1), names(data2))

  # Ensure only unique columns are selected
  common_cols1 <- c("Allele", common_markers[!common_markers %in% "Allele"])
  common_cols2 <- c("Allele", common_markers[!common_markers %in% "Allele"])

  merged1 <- full_join(data1[, common_cols1, drop = FALSE], data2[, common_cols2, drop = FALSE], by = "Allele", suffix = c(".1", ".2"))
  merged2 <- full_join(data2[, common_cols2, drop = FALSE], data1[, common_cols1, drop = FALSE], by = "Allele", suffix = c(".2", ".1"))

  merged1[is.na(merged1)] <- 0
  merged2[is.na(merged2)] <- 0

  merged1[merged1 == 0] <- minFreq
  merged2[merged2 == 0] <- minFreq

  merged1 <- merged1 %>%
    mutate(across(matches("\\.1$"), ~ ./sum(.))) %>%
    mutate(across(matches("\\.2$"), ~ ./sum(.)))
  merged2 <- merged2 %>%
    mutate(across(matches("\\.2$"), ~ ./sum(.))) %>%
    mutate(across(matches("\\.1$"), ~ ./sum(.)))

  kl1 <- sum(merged1[, grep("\\.1$", names(merged1))] * log10(merged1[, grep("\\.1$", names(merged1))] / merged1[, grep("\\.2$", names(merged1))]))
  kl2 <- sum(merged2[, grep("\\.2$", names(merged2))] * log10(merged2[, grep("\\.2$", names(merged2))] / merged2[, grep("\\.1$", names(merged2))]))

  return(list("KL from data1 to data2" = kl1, "KL from data2 to data1" = kl2))
}
