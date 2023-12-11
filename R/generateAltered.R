#library(pedtools)
#library(forrel)
#library(mispitools)

generate_altered <- function(kinship, zygosity, n_errors, n_madres=10000){

  targetID <- c()
  targetID['madre'] <- 3
  targetID['tio'] <- 5
  targetID['abuelo'] <- 5
  targetID['bisabuelo'] <- 7

  for (seed in 1:n_madres){ 
    print(paste("generando",kinship,seed))
    ped <- create_pedigree(kinship, seed)
    target_id <- targetID[kinship]
    hijos <- efficientProfileSim(ped, n_hijos, target_id, seed = seed+n_madres, numCores = n_cores)
    #al seed de generacion de hijos le agrego n_madres (si madres fue de 1 a 10 esto va de 11 a 20) para que no sean los mismos hijos del df sin alterar
    if (zygosity){
      hijos <- lapply(hijos, zygosity_alter_pedigree, target_id=target_id, n_errors)
    } else {
      hijos <- lapply(hijos, type_alter_pedigree, target_id=target_id, n_errors)
    }
    
    hijos_appnd <- markers_from_ped(hijos,target_id)
    hijos_split <- split_df(hijos_appnd, is_target=TRUE)
    hijos_split$target <-  if (zygosity) 'zygosity' else 'type'
    error_name <- if (zygosity) '_zygosity_errors_' else '_type_errors_'
    out_fname <- paste("./data/",kinship,"_only/",kinship,n_errors,error_name,seed,".csv", sep="")
    write.csv(hijos_split, out_fname, row.names=FALSE)
  }
}