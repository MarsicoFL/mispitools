#library(pedtools)
#library(forrel)
#library(mispitools)

# silent allele / drop out error (errors that change from heterozygous to homozygous)
zygosity_alter_pedigree <- function(ped, target_id, n_errors, verbose=FALSE){
  marker_names <- colnames(as.data.frame(ped))[-c(1,2,3,4)]
  markers_selected <- sample(marker_names,n_errors)
  if (verbose){
    print(markers_selected)
  }
  for (marker_name in markers_selected){
    genotype(ped, marker_name, target_id)[2] <- genotype(ped, marker_name, target_id)[1]
  }
  return(ped)
}