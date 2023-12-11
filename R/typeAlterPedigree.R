#library(pedtools)
#library(forrel)
#library(mispitools)

# typing error
type_alter_pedigree <- function(ped, target_id, n_errors, verbose=FALSE){
  marker_names <- colnames(as.data.frame(ped))[-c(1,2,3,4)]
  markers_selected <- sample(marker_names,n_errors)
  if (verbose){
    print(markers_selected)
  }
  for (marker_name in markers_selected){
    freqs <- afreq(getMarkers(ped, marker_name))
    possible_alleles <- alleles(getMarkers(ped, marker_name))[freqs > 0]
    altered_allele <- sample(possible_alleles,1)
    if (runif(1) < 0.5){
      genotype(ped, marker_name, target_id)[1] <- altered_allele
    } else {
      genotype(ped, marker_name, target_id)[2] <- altered_allele
    }  
  }
  return(ped)
}