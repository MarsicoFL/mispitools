#' Extract Marker Data from Pedigree Object
#'
#' @description
#' Internal function to extract genetic marker data for a specific individual
#' from a pedigree object. Converts the marker data to a data frame format
#' suitable for further analysis.
#'
#' @param pedigree A pedigree object with genetic markers (from pedtools/forrel).
#' @param target_id Character or numeric. The ID of the individual to extract
#'   marker data for.
#'
#' @return A data.frame with marker data for the target individual.
#'   Each row corresponds to a marker, columns contain allele information.
#'
#' @details
#' This is an internal helper function used by other mispitools functions
#' to process pedigree marker data. It removes the first 4 columns (typically
#' id, fid, mid, sex) and returns only the genetic marker columns.
#'
#' @keywords internal
#' @noRd

markers_from_ped <- function(pedigree, target_id) {
  dataframe_list <- list()

  for (marker in pedigree) {
    marker_df <- as.data.frame(marker)
    dataframe_list <- append(dataframe_list, list(marker_df[target_id, ]))
  }
  ped_mrkrs_appnd <- do.call(rbind, dataframe_list)
  ped_mrkrs_appnd <- ped_mrkrs_appnd[, -1:-4]

  return(ped_mrkrs_appnd)
}
