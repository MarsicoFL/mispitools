#library(pedtools)
#library(forrel)
#library(mispitools)

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
