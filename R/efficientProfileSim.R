#library(pedtools)
#library(forrel)
#library(mispitools)

# memory efficient profileSim
efficientProfileSim <- function(ped, n_hijos, target_id, seed, numCores, n_fractions=10, markers = NULL){
  seed <- seed*100
  if(n_hijos%%n_fractions != 0){
    print('ALERTA: LA CANTIDAD DE FRACCIONES DEBE SER DIVISOR DE LA CANTIDAD DE HIJOS')
  }
  
  peds <- list()
  
  for (i_frac in 1:n_fractions){
    frac_ped <- profileSim(ped, n_hijos/n_fractions, 
                           target_id, seed = seed+i_frac, numCores = numCores,
                           markers=markers, verbose=FALSE)
    peds <- c(peds, frac_ped)
  }
  
  return(peds)
}
