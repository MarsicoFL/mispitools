#library(pedtools)
#library(forrel)
#library(mispitools)

markers_from_ped <- function(pedigree,target_id){
  dataframe_list <- list()
  
  for (marker in pedigree) {
    marker_df <- as.data.frame(marker)
    dataframe_list <- append(dataframe_list, list(marker_df[target_id,]))
  }
  
  ped_mrkrs_appnd <- do.call(rbind, dataframe_list)
  ped_mrkrs_appnd <- ped_mrkrs_appnd[,-1:-4]
  
  return(ped_mrkrs_appnd)
}

split_column <- function(column) {
  split_data <- strsplit(column, "/")
  first_numbers <- sapply(split_data, function(x) as.numeric(x[1]))
  second_numbers <- sapply(split_data, function(x) as.numeric(x[2]))
  return(data.frame(First = first_numbers, Second = second_numbers))
}

split_df <- function(ped_mrkrs_appnd, is_target){
  ped_columns <- lapply(ped_mrkrs_appnd, split_column)
  ped_split <- do.call(cbind, ped_columns)
  ped_split$target <- is_target
  
  return(ped_split)
}

shuffle_and_save <- function(ped_split, sample_split, name, number){
  ped_sample <- rbind(ped_split, sample_split)
  random_indices <- sample(2*n_hijos)
  ped_sample <- ped_sample[random_indices,]
  write.csv(ped_sample, paste("./data/",name,"/", name,"_",number,".csv", sep=""), row.names=FALSE)
}

mix_and_save <- function(filename_true, filename_false){
  true <- read.csv(filename_true)
  false <- read.csv(filename_false)
  
  true$target <- TRUE
  false$target <- FALSE
  
  mixed <- rbind(true, false)
  random_indices <- sample(2*n_hijos)
  mixed <- mixed[random_indices,]
  
  out_suffix <- strsplit(filename_true, "/")[[1]][4]
  write.csv(mixed, paste("./data/LRs/mixed/",out_suffix, sep=""), row.names=FALSE)
}

# memory efficient profileSim
efficientProfileSim <- function(ped, n_hijos, target_id, seed, numCores, n_fractions=10, markers = NULL){
  seed <- seed*100
  if(n_hijos%%n_fractions != 0){
    print('ALERTA: LA CANTIDAD DE FRACCIONES DEBE SER DIVISOR DE LA CANTIDAD DE HIJOS')
  }
  
  peds <- list()
  
  for (i_frac in 1:n_fractions){
    frac_ped <- profileSim(ped, n_hijos/n_fractions, target_id, seed = seed+i_frac, numCores = numCores, markers=markers, verbose=FALSE)
    peds <- c(peds, frac_ped)
  }
  
  return(peds)
}


create_pedigree <- function(kinship, seed){
  if (kinship == 'madre'){
    ped <- linearPed(1)
    ped <- profileSim(ped, 1, 2,  markers = getfreqs(Argentina), seed = seed)
  } else if (kinship == 'tio'){
    ped <- linearPed(2) 
    ped <- addChildren(ped, 1, 2)
    ped <- profileSim(ped, 1, 6, markers = getfreqs(Argentina), seed = seed)
  } else if (kinship == 'abuelo'){
    ped <- linearPed(2)
    ped <- profileSim(ped, 1, 1, markers = getfreqs(Argentina), seed = seed)
  } else if (kinship == 'bisabuelo'){
    ped <- linearPed(3)
    ped <- profileSim(ped, 1, 1,  markers = getfreqs(Argentina), seed = seed)
  }
  return(ped)
}


targetID <- function(kinship){
  if (kinship == 'madre'){
    target_id <- 3
  } else if (kinship == 'tio'){
    target_id <- 5
  } else if (kinship == 'abuelo'){
    target_id <- 5
  } else if (kinship == 'bisabuelo'){
    target_id <- 7
  } 
  return(target_id)
}

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


generate_altered <-function(kinship, zygosity, n_errors){
  for (seed in 1:n_madres){ 
    print(paste("generando",kinship,seed))
    ped <- create_pedigree(kinship, seed)
    target_id <- targetID(kinship)
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
