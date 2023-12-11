#library(pedtools)
#library(forrel)
#library(mispitools)


split_df <- function(ped_mrkrs_appnd, is_target){

  split_column <- function(column) {
    split_data <- strsplit(column, "/")
    first_numbers <- sapply(split_data, function(x) as.numeric(x[1]))
    second_numbers <- sapply(split_data, function(x) as.numeric(x[2]))

    return(data.frame(First = first_numbers, Second = second_numbers))
  }

  ped_columns <- lapply(ped_mrkrs_appnd, split_column)
  ped_split <- do.call(cbind, ped_columns)
  ped_split$target <- is_target
  
  return(ped_split)
}
