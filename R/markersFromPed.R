#library(pedtools)
#library(forrel)
#library(mispitools)


markers_from_ped <- function(pedigree, target_id){
  dataframe_list <- list()

  for (marker in pedigree) {
    marker_df <- as.data.frame(marker)
    dataframe_list <- append(dataframe_list, list(marker_df[target_id, ]))
  }
  ped_mrkrs_appnd <- do.call(rbind, dataframe_list)
  ped_mrkrs_appnd <- ped_mrkrs_appnd[, -1:-4]
  
  return(ped_mrkrs_appnd)
}

