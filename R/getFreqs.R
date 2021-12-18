#' Function for getting STR allele frequencies from 48 world populations.
#'
#' @param region select the place of the allele frequency database. Possible values are listed: "Argentina", "Asia","Algeria (M’zab) – Mozabite","Austria","Belgium","Bosnia and Herzegowina","Brazil – Karitiana", "Brazil – Suruí", "Bugainville – NAN Melanesian","Cambodia – Cambodian", "Central African Republic – Biaka Pygmies","China – Dai", "China – Lahu", "China – Yizu", "China – Miaozu", "China – Han", "China – Tu", "China – Oroqen", "China – Hezhen", "China – Uygur", "China – Daur", "China – She", "China – Mongola", "China – Xibo", "China – Naxi", "China – Tujia", "Colombia – Colombian", "Czech Republic", "Democratic Republic of the Congo – Mbuti Pygmies", "Denmark", "Dominican Republic", "Europe", "Finland","France","Germany","Greece","Hungary","Ireland","Israel (Negev) – Bedouin","Israel (Central) – Palestinian","Israel (Carmel) – Druze","Japan – Japanese","Kenya – Bantu (North East)","Mexico – Maya","Mexico – Pima","Montenegro","Namibia – San","New Guinea – Papuan","Nigeria – Yoruba", "Norway", "Pakistan – Burusho", "Pakistan – Sindhi", "Pakistan – Balochi","Pakistan – Kalash","Pakistan – Pathan", "Pakistan – Makrani","Pakistan – Hazara","Poland","Saudi Arabia","Senegal – Mandenka","Siberia - Yakut","Slovenia","Somalia","South Africa – Bantu","Spain","Spain (North-West)","Sweden","Switzerland","Thailand","U.S. All (length-based)","U.S. All (sequence-based)","U.S. African American (length-based)","U.S. African American (sequence-based)","U.S. Asian (sequence-based)", "U.S. Asian (length-based)","U.S. Caucasian (sequence-based)","U.S. Caucasian (length-based)","U.S. Hispanic (length-based)","U.S. Hispanic (sequence-based)"
#'
#' @return An allele frequency database. 
#' @export
#' @import rjson
#' @import httr
#'
#' @examples
#' library(rjson) 
#' getFreqs("Argentina")

getfreqs <- function(region){
#Get alelos de leapdna:
indice.es.alelos <- data.frame(country = c("Argentina",
                                              "Asia",
                                              "Algeria (M’zab) – Mozabite",
                                              "Austria",
                                              "Belgium",
                                              "Bosnia and Herzegowina",
                                              "Brazil – Karitiana",
                                              "Brazil – Suruí",
                                              "Bugainville – NAN Melanesian",
                                              "Cambodia – Cambodian",
                                              "Central African Republic – Biaka Pygmies",
                                              "China – Dai",
                                              "China – Lahu",
                                              "China – Yizu",
                                              "China – Miaozu",
                                              "China – Han",
                                              "China – Tu", 
                                              "China – Oroqen",
                                              "China – Hezhen",
                                              "China – Uygur",
                                              "China – Daur",
                                              "China – She",
                                              "China – Mongola",
                                              "China – Xibo",
                                              "China – Naxi",
                                              "China – Tujia",
                                              "Colombia – Colombian",
                                              "Czech Republic",
                                              "Democratic Republic of the Congo – Mbuti Pygmies",
                                              "Denmark",
                                              "Dominican Republic",
                                              "Europe",
                                              "Finland",
                                              "France",
                                              "Germany",
                                              "Greece",
                                              "Hungary",
                                              "Ireland",
                                              "Israel (Negev) – Bedouin",
                                              "Israel (Central) – Palestinian",
                                              "Israel (Carmel) – Druze",
                                              "Japan – Japanese",
                                              "Kenya – Bantu (North East)",
                                              "Mexico – Maya",
                                              "Mexico – Pima",
                                              "Montenegro",
                                              "Namibia – San",
                                              "New Guinea – Papuan",
                                              "Nigeria – Yoruba",
                                              "Norway",
                                              "Pakistan – Burusho",
                                              "Pakistan – Sindhi",
                                              "Pakistan – Balochi",
                                              "Pakistan – Kalash",
                                              "Pakistan – Pathan",
                                              "Pakistan – Makrani",
                                              "Pakistan – Hazara",
                                              "Poland",
                                              "Saudi Arabia",
                                              "Senegal – Mandenka",
                                              "Siberia - Yakut",
                                              "Slovenia",
                                              "Somalia",
                                              "South Africa – Bantu",
                                              "Spain",
                                              "Spain (North-West)",
                                              "Sweden",
                                              "Switzerland",
                                              "Thailand",
                                              "U.S. All (length-based)",
                                              "U.S. All (sequence-based)",
                                              "U.S. African American (length-based)",
                                              "U.S. African American (sequence-based)",
                                              "U.S. Asian (sequence-based)",
                                              "U.S. Asian (length-based)",
                                              "U.S. Caucasian (sequence-based)",
                                              "U.S. Caucasian (length-based)",
                                              "U.S. Hispanic (length-based)",
                                              "U.S. Hispanic (sequence-based)"
                                              ),
                                    
                                    urls = c('data/freq_Arg.rda',
                                             "https://api.leapdna.org/studies/strider_asia_all.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_algeria.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_austria.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_belgium.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_bosina_and_herzegowina.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_brazil_karitiana.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_brazil_surui.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_bougainville.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_cambodia.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_central_african_republic.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_dai.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_lahu.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_yizu.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_miaozu.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_han.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_tu.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_oroqen.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_hezhen.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_uygur.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_daur.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_she.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_mongola.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_xibo.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_naxi.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_china_tujia.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_colombia.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_czech_republic.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_dpcongo_mbuti.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_denmark.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_dominican_republic.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_europe_all.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_finland.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_france.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_germany.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_greece.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_hungary.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_ireland.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_israel_bedouin.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_israel_palestinian.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_israel_druze.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_japan.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_kenya_bantu.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_mexico_maya.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_mexico_pima.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_montenegro.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_namibia_san.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_new_guinea_papuan.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_nigeria_yoruba.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_norway.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_pakistan_burusho.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_pakistan_sindhi.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_pakistan_balochi.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_pakistan_kalash.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_pakistan_pathan.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_pakistan_makrani.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_pakistan_hazara.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_poland.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_saudi_arabia.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_senegal_mandenka.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_siberia_yakut.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_slovenia.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_somalia.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_south_africa_bantu.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_spain.leapdna.json",
                                             "https://api.leapdna.org/studies/popstr_spain_nw.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_sweden.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_switzerland.leapdna.json",
                                             "https://api.leapdna.org/studies/strider_thailand.leapdna.json",
                                             "https://api.leapdna.org/studies/nist1036_us.leapdna.json",
                                             "https://api.leapdna.org/studies/nist1036_seq_us.leapdna.json",
                                             "https://api.leapdna.org/studies/nist1036_us_african_american.leapdna.json",
                                             "https://api.leapdna.org/studies/nist1036_seq_us_african_american.leapdna.json",
                                             "https://api.leapdna.org/studies/nist1036_seq_us_asian.leapdna.json",
                                             "https://api.leapdna.org/studies/nist1036_us_asian.leapdna.json",
                                             "https://api.leapdna.org/studies/nist1036_seq_us_caucasian.leapdna.json",
                                             "https://api.leapdna.org/studies/nist1036_us_caucasian.leapdna.json",
                                             "https://api.leapdna.org/studies/nist1036_us_hispanic.leapdna.json",
                                             "https://api.leapdna.org/studies/nist1036_seq_us_hispanic.leapdna.json")
                                    ) %>% 
  arrange(country)
  
  indice = indice.es.alelos
  leap.url <- indice[indice$country == region, "urls"]
  
  if(region == "Argentina"){
    prueba <- read.csv(file = leap.url)
    prueba.lista <- as.list(prueba)
    for(i in 2:length(prueba.lista)){
      names(prueba.lista[[i]]) <- prueba.lista[[1]]
    }
    prueba.lista$Allele <- NULL
    
    return(prueba.lista)
  
  }else{
    leap <- httr::GET(leap.url)
    df.leap <- fromJSON(rawToChar(leap$content))$loci
    
    lista.alelos <- list()
    for(i in unique(df.leap$name)){
      df.alelo <- df.leap[df.leap$name == i, "alleles"][[1]]
      vector.alelo <- df.alelo$frequency
      names(vector.alelo) <- df.alelo$name
      lista.alelos[[i]] <- vector.alelo
    }
    
    return(lista.alelos)
   
  } 
  }
