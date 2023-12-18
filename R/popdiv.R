# FUNCTION: Normalize Allele Frequencies
#
# This function checks and normalizes allele frequencies in a dataset for each genetic marker.
# It ensures the sum of allele frequencies for each marker is approximately 1.
#
# @param data A data frame with genetic markers as columns and allele frequencies as values.
# @return A data frame with normalized allele frequencies, if necessary.

normalizeAlleleFrequencies <- function(data) {
  markers <- colnames(data)[-1]
  not_normalized_markers <- sapply(markers, function(marker) {
    sum_val <- sum(data[[marker]], na.rm = TRUE)
    return(abs(sum_val - 1) >= 1e-6)
  })

  for (marker in markers[not_normalized_markers]) {
    total_frequency <- sum(data[[marker]], na.rm = TRUE)
    if (total_frequency != 0) {
      data[[marker]] <- data[[marker]] / total_frequency
    }
  }

  return(data)
}

# FUNCTION: Compute bidirectional KL Divergence
#
# This function calculates the Kullback-Leibler Divergence for common genetic markers between two data sets of allele frequencies.
# It first identifies common markers between the two data sets, excluding specified columns.
# Returns NA if there are no common markers.
#
# @param data_1 First data set of allele frequencies.
# @param data_2 Second data set of allele frequencies.
# @param ignore_columns Columns to exclude from the comparison (default is "X").
# @return A numeric value representing the KL divergence or NA if there are no common markers.

bidirectionalKL <- function(data_1, data_2, ignore_columns = c("X")) {
  colnames_data_1 <- setdiff(colnames(data_1), ignore_columns)
  colnames_data_2 <- setdiff(colnames(data_2), ignore_columns)
  common_markers <- intersect(colnames_data_1, colnames_data_2)

  if (length(common_markers) == 0) {
    return(NA)
  }

  P <- data_1[, common_markers, drop = FALSE]
  Q <- data_2[, common_markers, drop = FALSE]
  P[is.na(P)] <- 1e-16
  Q[is.na(Q)] <- 1e-16

  kl_divergence <- sum(sapply(1:nrow(P), function(i) {
    row_P <- P[i, ]
    row_Q <- Q[i, ]
    sum(row_P * log(row_P / row_Q), na.rm = TRUE)
  }), na.rm = TRUE)

  return(kl_divergence)
}

# FUNCTION: Test KL Divergence for All Pairs in Dataset
#
# Computes the KL divergence for each pair of countries in the given dataset.
# It returns a matrix with KL divergence values or NA for each pair of countries,
# with NA indicating pairs with no common markers.
#
# @param datasets A list of data frames, each representing allele frequencies for a different country.
# @return A matrix containing KL divergence values or NA for each pair of countries.
allPairsKL <- function(datasets) {
  country_names <- names(datasets)
  results_matrix <- matrix(nrow = length(country_names), ncol = length(country_names),
                           dimnames = list(country_names, country_names))

  for (i in 1:length(country_names)) {
    for (j in 1:length(country_names)) {
      if (i != j) {
        kl_divergence <- bidirectionalKL(datasets[[country_names[i]]], datasets[[country_names[j]]])
        results_matrix[i, j] <- kl_divergence
      }
    }
  }

  return(results_matrix)
}
